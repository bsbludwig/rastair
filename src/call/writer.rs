use crate::{
    bed::rastair1::{BedRecordsConvertParams, Rastair1BedFormat},
    call::CallParams,
    metrics::PileupMetrics,
    sequence::ChunkRegion,
    utils::logging::ThisIsABug as _,
    vcf::{self, VcfOutputFilter},
};
use color_eyre::{Result, eyre::Context as _};
use ordered_channel::Receiver;
use smol_str::SmolStr;
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

    // Spawn the actual VCF writer thread. Everything in here is driven by the
    // incoming records from the processing threads.
    //
    // The result returned from this thread is evaluated when the handle is joined.
    thread::Builder::new()
        .name("writer".to_string())
        .spawn(move || -> Result<()> {
            let filters = VcfOutputFilter { reject_low_quality_variants: true };

            let mut last_seen = LastSeen::default();

            // Since we only have the region index to ensure order, each
            // processing thread will send a vector of VCF records when it's
            // done with a region.
            for records in vcf_receiver {
                'current_batch: for mut record in records {
                    if !last_seen.is_new(record.contig(), record.pos()) {
                        continue 'current_batch;
                    }

                    if filters.reject_low_quality_variants {
                        record.alts.retain(|alt| alt.filters.filters.is_empty());
                    }

                    let mut vcf_record: Option<vcf::Record> = None;

                    if let Some(vcf_writer) = vcf_writer.as_mut()
                        && vcf_filter.matches(&record)
                    {
                        use crate::io::vcf_writer::Writer;
                        match vcf_writer {
                            Writer::Vcf(writer) => {
                                let vcf = record
                                    .to_vcf_record()
                                    .wrap_err("Failed to convert metrics to VCF record")
                                    .this_is_a_bug()?;
                                writer.add(&vcf).wrap_err("Failed to write VCF record")?;
                                vcf_record = Some(vcf);
                            }
                            Writer::MessagePack(writer) => {
                                writer
                                    .add(&record)
                                    .wrap_err("Failed to write MessagePack VCF record")?;
                            }
                        }
                    }

                    if let Some(bed_writer) = bed_writer.as_mut()
                        && (*record.pos_metrics.cpg || record.forms_denovo())
                    {
                        let vcf = if let Some(vcf) = vcf_record {
                            vcf
                        } else {
                            record
                                .to_vcf_record()
                                .wrap_err("Failed to convert metrics to VCF record")
                                .this_is_a_bug()?
                        };
                        if let Some(bed_record) = Rastair1BedFormat::from_record(&vcf, &bed_params)
                            .wrap_err("Failed to convert VCF record to BED format")
                            .this_is_a_bug()?
                        {
                            bed_writer
                                .write_record(&bed_record)
                                .wrap_err("Failed to write record to BED")?;
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

/// The segments we get have some overlap between them, so we need
/// to ensure that we don't write the same record multiple times.
#[derive(Default)]
struct LastSeen {
    contig: Option<SmolStr>,
    pos: Option<u32>,
}

impl LastSeen {
    /// If this is new, returns true and updates the last seen record
    fn is_new(&mut self, contig: SmolStr, pos: u32) -> bool {
        if self.contig.as_ref() == Some(&contig) && self.pos >= Some(pos) {
            false
        } else {
            self.contig = Some(contig);
            self.pos = Some(pos);
            true
        }
    }
}
