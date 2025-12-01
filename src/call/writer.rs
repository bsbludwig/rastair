use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat},
    call::CallParams,
    metrics::PileupMetrics,
    sequence::ChunkRegion,
    utils::logging::ThisIsABug as _,
};
use color_eyre::{Result, eyre::Context as _};
use ordered_channel::Receiver;
use std::thread;
use tracing::info;

/// Build the VCF writer thread
pub fn writer_thread(
    params: &CallParams,
    regions: &[ChunkRegion],
    vcf_receiver: Receiver<Vec<PileupMetrics>>,
) -> Result<thread::JoinHandle<Result<()>>> {
    let vcf_output = params.vcf.vcf.clone();
    let vcf_filter = params.record_filters.clone();
    let metadata = [
        format!("rastairVersion={}", env!("CARGO_PKG_VERSION")),
        format!("rastairCommand={}", std::env::args().skip(1).collect::<Vec<_>>().join(" ")),
        format!(
            "rastairConfig={}",
            serde_json::to_string(params)
                .wrap_err("Failed to serialize config to JSON")
                .this_is_a_bug()?
        ),
        format!("reference={}", params.segments.fasta_file),
    ];
    let mut vcf_writer =
        params.vcf.writer(regions, &metadata).wrap_err("Failed to create VCF writer")?;

    let bed = params.bed.clone();
    let mut bed_writer = bed.writer().wrap_err("Failed to create BED writer")?;
    let bed_params =
        BedRecordsConvertParams { ml_threshold: params.ml.ml, filters: bed.filters.clone() };

    let ml_threshold = params.ml.threshold();

    // Spawn the actual VCF writer thread. Everything in here is driven by the
    // incoming records from the processing threads.
    //
    // The result returned from this thread is evaluated when the handle is joined.
    thread::Builder::new()
        .name("writer".to_string())
        .spawn(move || -> Result<()> {
            let span =
                tracing::debug_span!("writer", vcf=%vcf_output.is_some(), bed=%bed_writer.is_some());
            let _guard = span.enter();

            // Since we only have the region index to ensure order, each
            // processing thread will send a vector of VCF records when it's
            // done with a region.
            for records in vcf_receiver {
                let span = tracing::debug_span!("recv_records", records=%records.len());
                let _guard = span.enter();

                for record in records {
                    // Write BED record if requested
                    if let Some(bed_writer) = bed_writer.as_mut()
                        && let Some(bed_record) = Rastair1BedFormat::from_metrics(&record, &bed_params)
                            .wrap_err("Failed to convert record to BED format")
                            .this_is_a_bug()?
                        {
                            bed_writer
                                .write_record(&bed_record)
                                .wrap_err("Failed to write record to BED")?;
                        }

                    if let Some(vcf_writer) = vcf_writer.as_mut()
                        && vcf_filter.matches(&record, ml_threshold)
                    {
                        use crate::io::vcf_writer::Writer;
                        match vcf_writer {
                            Writer::Vcf(writer) => {
                                let mut records = record
                                    .to_vcf_records(ml_threshold)
                                    .wrap_err("Failed to convert metrics to VCF record")
                                    .this_is_a_bug()?;
                                records.filter(&vcf_filter);
                                for vcf_record in records.iter() {
                                    writer.add(&vcf_record).wrap_err("Failed to write VCF record")?;
                                }
                            }
                            Writer::MessagePack(writer) => {
                                writer
                                    .add(&record)
                                    .wrap_err("Failed to write MessagePack VCF record")?;
                            }
                        }
                    }
                }
            }

            if let Some(vcf_output) = vcf_output.as_ref() {
                drop(vcf_writer);
                info!(file = %vcf_output, "Wrote VCF output");
            }
            if let Some(bed_output) = bed.bed.as_ref()
                && let Some(bed_writer) = bed_writer
            {
                bed_writer.close().wrap_err("Failed to close BED writer")?;
                info!(file = %bed_output, "Wrote BED output");
            }
            Ok(())
        })
        .wrap_err("Failed to spawn VCF writer thread")
}
