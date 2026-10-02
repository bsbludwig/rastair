use crate::{
    bed::{rastair1::Rastair1BedFormat, writer::BedWriter},
    call::{CallParams, RecordFilters, variant_calling::ErrorModel},
    io::vcf_writer::Writer,
    metrics::PileupMetrics,
    runtime::{
        fault_injection, partial_output::PartialOutput, progress::ProgressTracker, segments,
    },
    sequence::ChunkRegion,
    utils::{Probability, logging::ThisIsABug as _},
};
use color_eyre::{Result, eyre::Context as _};
use tracing::{instrument, trace};

/// Writes finished segments, in order, to the VCF and BED outputs.
pub struct SegmentWriter {
    outputs: Vec<PartialOutput>,
    vcf_writer: Option<Writer>,
    bed_writer: Option<BedWriter<Rastair1BedFormat>>,
    record_filter: RecordFilters,
    ml_threshold: Option<Probability>,
    error_model: ErrorModel,
    progress: ProgressTracker,
}

impl SegmentWriter {
    pub fn new(params: &CallParams, regions: &[ChunkRegion]) -> Result<Self> {
        let vcf_output = params.vcf.vcf.as_ref().map(|path| PartialOutput::new("VCF output", path));
        let bed_output = params.bed.bed.as_ref().map(|path| PartialOutput::new("BED output", path));
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
        let vcf_writer = vcf_output
            .as_ref()
            .map(|output| params.vcf.writer(&output.write_path()?, regions, &metadata))
            .transpose()
            .wrap_err("Failed to create VCF writer")?;
        let bed_writer = bed_output
            .as_ref()
            .map(|output| params.bed.writer(&output.write_path()?))
            .transpose()
            .wrap_err("Failed to create BED writer")?;

        Ok(Self {
            outputs: vcf_output.into_iter().chain(bed_output).collect(),
            vcf_writer,
            bed_writer,
            record_filter: params.record_filters.clone(),
            ml_threshold: params.ml.threshold(),
            error_model: params.variant_calling.error_model,
            progress: ProgressTracker::new(regions.len()),
        })
    }

    #[instrument(level = "debug", skip_all, fields(records = records.len()))]
    pub fn write(&mut self, records: Vec<PileupMetrics>) -> Result<()> {
        for record in records {
            if !self.record_filter.matches(&record) {
                trace!(pos=%record.contig_pos(), "Record did not pass filters, skipping");
                continue;
            }

            if let Some(bed_writer) = self.bed_writer.as_mut()
                && let Some(bed_record) = Rastair1BedFormat::from_metrics(&record)
                    .wrap_err("Failed to convert record to BED format")
                    .this_is_a_bug()?
            {
                bed_writer.write_record(&bed_record).wrap_err("Failed to write record to BED")?;
            }

            match self.vcf_writer.as_mut() {
                Some(Writer::Vcf(writer)) => {
                    writer
                        .emit(&record, self.ml_threshold, &self.error_model, &self.record_filter)
                        .wrap_err("Failed to write VCF record")?;
                }
                Some(Writer::MessagePack(writer)) => {
                    writer.add(&record).wrap_err("Failed to write MessagePack VCF record")?;
                }
                None => {}
            }
        }

        self.progress.segment_done();
        fault_injection::fault_point(fault_injection::FaultPoint::Writer)
    }

    /// Closes every output, even if one fails to close, so each holds a valid
    /// prefix. Also returns the outputs, to commit once the run succeeded.
    pub fn close(self) -> (Result<()>, Vec<PartialOutput>) {
        let vcf_closed = match self.vcf_writer {
            Some(Writer::Vcf(mut writer)) => {
                writer.finish().wrap_err("Failed to finish VCF output")
            }
            Some(Writer::MessagePack(writer)) => writer.finish(),
            None => Ok(()),
        };
        let bed_closed = self
            .bed_writer
            .map_or(Ok(()), |writer| writer.close().wrap_err("Failed to close BED writer"));
        (segments::first_error(vcf_closed, bed_closed), self.outputs)
    }
}
