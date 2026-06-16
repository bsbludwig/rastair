use crate::{
    call::{RecordFilters, variant_calling::ErrorModel},
    io::{
        formats::FromFileExtension,
        mpk::{MessagePackWriter, format::MpkVcfHeader},
    },
    metrics::PileupMetrics,
    sequence::ChunkRegion,
    utils::{cli, logging::ThisIsABug as _},
    vcf::{Contig, FieldConfig, Schema, emit_pileup, register},
};
use better_default::Default;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use clio::ClioPath;
use color_eyre::eyre::{ContextCompat, Result, WrapErr};
use rustc_hash::FxHashSet;
use seqair::vcf::{OutputFormat, Ready, Writer as SeqWriter};
use seqair_types::{Probability, SmolStr};
use std::{ffi::OsStr, io::Write, num::NonZeroUsize};
use tracing::{debug, warn};

#[derive(Debug, Clone, Default, clap::Parser)]
pub struct VcfParams {
    /// VCF/BCF output file path (use - to write to stdout)
    ///
    /// Format is guessed based on the file extension:
    /// `.vcf` for VCF (uncompressed),
    /// `.vcf.gz` for VCF (compressed),
    /// `.bcf` for BCF (compressed)
    /// `.mpk.lz4` for internal format (Message Pack, LZ4-compressed)
    #[arg(short = 'o', long, required = false, default_missing_value = "-", num_args = 0..=1)]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub vcf: Option<ClioPath>,

    /// Number of threads to use for writing (and compressing) VCF files
    ///
    /// This is subtracted from `--threads` but never below 1. Adjust this if
    /// you think that VCF writing is a bottleneck, e.g. when the output files
    /// contain a lot of positions.
    // Default value chosen after profiling on a machine with 14 cores.
    #[arg(long, default_value = "1")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[default(NonZeroUsize::new(1).expect("1 > 0"))]
    pub vcf_threads: NonZeroUsize,

    /// Additional INFO fields to include in VCF output (comma-separated VCF field IDs)
    ///
    /// By default, only a minimal set is included.
    #[arg(
        long,
        value_delimiter = ',',
        value_parser = PossibleValuesParser::new(crate::vcf::InfoFieldId::ALL_IDS)
            .map(|s| s.parse::<crate::vcf::InfoFieldId>().unwrap())
    )]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub vcf_info_fields: Vec<crate::vcf::InfoFieldId>,

    /// Additional FORMAT fields to include in VCF output (comma-separated VCF field IDs)
    ///
    /// By default, only a minimal set is included.
    #[arg(
        long,
        value_delimiter = ',',
        value_parser = PossibleValuesParser::new(crate::vcf::FormatFieldId::ALL_IDS)
            .map(|s| s.parse::<crate::vcf::FormatFieldId>().unwrap())
    )]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub vcf_format_fields: Vec<crate::vcf::FormatFieldId>,

    // Include all possible fields in VCF
    #[arg(long, default_value_t = false, conflicts_with_all = &["vcf_info_fields", "vcf_format_fields"])]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub vcf_all_fields: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Vcf(VcfFormat),
    /// Rastair-internal `MessagePack` format, always LZ4-compressed
    MessagePack,
}

impl clap::ValueEnum for Format {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            Format::Vcf(VcfFormat::Vcf),
            Format::Vcf(VcfFormat::VcfCompressed),
            Format::Vcf(VcfFormat::Bcf),
            Format::MessagePack,
        ]
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        match self {
            Format::Vcf(format) => format.to_possible_value(),
            Format::MessagePack => Some(clap::builder::PossibleValue::new("mpk.lz4")),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum VcfFormat {
    /// Text-based VCF format (.vcf)
    Vcf,
    /// Compressed text-based VCF format (.vcf.gz)
    VcfCompressed,
    /// Binary VCF format (.bcf)
    Bcf,
}

impl From<VcfFormat> for OutputFormat {
    fn from(format: VcfFormat) -> Self {
        match format {
            VcfFormat::Vcf => OutputFormat::Vcf,
            VcfFormat::VcfCompressed => OutputFormat::VcfGz,
            VcfFormat::Bcf => OutputFormat::Bcf,
        }
    }
}

impl VcfParams {
    /// Create a new instance of `Params` with the specified VCF output path.
    pub fn guess_format(&self) -> Format {
        let Some(vcf_output) = &self.vcf else {
            // No VCF output, so the format doesn't matter.
            return Format::Vcf(VcfFormat::Vcf);
        };

        if vcf_output.is_std() {
            return Format::Vcf(VcfFormat::Vcf);
        }

        let Some(name) = vcf_output.file_name().and_then(OsStr::to_str) else {
            warn!(
                filename=%vcf_output,
                "No file name found in VCF output path, defaulting to VCF format without compression."
            );
            return Format::Vcf(VcfFormat::Vcf);
        };

        let Some(format) = Format::from_file_extension(name) else {
            warn!(
                filename=%vcf_output,
                "Could not determine format from file extension, defaulting to VCF format without compression."
            );
            return Format::Vcf(VcfFormat::Vcf);
        };

        format
    }

    /// Build the field configuration from the CLI flags.
    pub fn field_config(&self) -> FieldConfig {
        let config = FieldConfig::default();
        if self.vcf_all_fields {
            config.with_all_fields()
        } else {
            config.with_field_ids(&self.vcf_info_fields, &self.vcf_format_fields)
        }
    }

    pub fn writer(&self, regions: &[ChunkRegion], metadata: &[String]) -> Result<Option<Writer>> {
        let Some(_) = &self.vcf else {
            return Ok(None);
        };

        let contigs = header_contigs(regions);
        let samples = vec![SmolStr::new("sample")]; // Note: we only deal with one sample for now

        let format = match self.guess_format() {
            Format::MessagePack => {
                let writer = self
                    .create_mpk_writer(contigs, samples, metadata)
                    .wrap_err("Failed to create MessagePack writer")?;
                return Ok(Some(Writer::MessagePack(
                    writer.wrap_err("No VCF output path present").this_is_a_bug()?,
                )));
            }
            Format::Vcf(f) => f.into(),
        };

        Ok(Some(Writer::Vcf(
            self.seqair_writer(&contigs, &samples, metadata, format)
                .wrap_err("Failed to create VCF writer")?
                .wrap_err("No VCF output path present")
                .this_is_a_bug()?,
        )))
    }

    /// Build a seqair-backed VCF/BCF writer for the configured output path.
    pub fn seqair_writer(
        &self,
        contigs: &[Contig],
        samples: &[SmolStr],
        metadata: &[String],
        format: OutputFormat,
    ) -> Result<Option<SeqairVcfWriter>> {
        let Some(vcf_output) = &self.vcf else {
            return Ok(None);
        };

        debug!(target=?vcf_output.display(), ?format, "Creating VCF writer");

        let (header, schema) =
            register(contigs, samples, metadata).wrap_err("Failed to build VCF header")?;

        let inner: Box<dyn Write + Send> = Box::new(
            vcf_output
                .clone()
                .create()
                .wrap_err_with(|| format!("Failed to create output {vcf_output}"))?,
        );

        let writer = SeqWriter::new(inner, format)
            .write_header(&header)
            .wrap_err("Failed to write VCF header")?;

        Ok(Some(SeqairVcfWriter {
            writer: Some(writer),
            schema,
            config: self.field_config(),
            last_contig: None,
        }))
    }

    pub fn create_mpk_writer(
        &self,
        contigs: Vec<Contig>,
        samples: Vec<SmolStr>,
        metadata: &[String],
    ) -> Result<Option<MessagePackWriter>> {
        let Some(path) = &self.vcf else {
            return Ok(None);
        };
        warn!(
            %path,
            "MessagePack format only for internal use, no stability guarantees",
        );
        let mut w = MessagePackWriter::new(path)
            .wrap_err_with(|| format!("Failed to create MessagePack writer for {path}"))?;

        w.add_metadata(MpkVcfHeader { contigs, samples, metadata: metadata.to_owned() })?;

        Ok(Some(w))
    }
}

/// A seqair-backed VCF/BCF writer plus the resolved schema and field selection.
pub struct SeqairVcfWriter {
    writer: Option<SeqWriter<Box<dyn Write + Send>, Ready>>,
    schema: Schema,
    config: FieldConfig,
    last_contig: Option<seqair::vcf::ContigId>,
}

impl SeqairVcfWriter {
    /// Encode every VCF record this pileup produces.
    pub fn emit(
        &mut self,
        pileup: &PileupMetrics,
        ml_threshold: Option<Probability>,
        error_model: &ErrorModel,
        record_filter: &RecordFilters,
    ) -> Result<()> {
        let name = pileup.contig_name();
        if self.last_contig.as_ref().is_none_or(|c| c.name() != name) {
            self.last_contig = Some(
                self.schema
                    .contig(name)
                    .wrap_err_with(|| format!("Contig {name} not in header"))?
                    .clone(),
            );
        }
        let contig = self.last_contig.as_ref().wrap_err("contig cache unset").this_is_a_bug()?;

        let writer =
            self.writer.as_mut().wrap_err("VCF writer already finished").this_is_a_bug()?;
        emit_pileup(
            pileup,
            &self.schema,
            contig,
            &self.config,
            ml_threshold,
            error_model,
            record_filter,
            writer,
        )
    }

    /// Flush the BGZF EOF block / finalize the stream. Must be called once.
    pub fn finish(&mut self) -> Result<()> {
        if let Some(writer) = self.writer.take() {
            writer.finish().wrap_err("Failed to finish VCF output")?;
        }
        Ok(())
    }
}

pub enum Writer {
    Vcf(SeqairVcfWriter),
    MessagePack(MessagePackWriter),
}

/// Build the VCF header contig list: one entry per distinct contig touched by
/// `regions`, carrying the *true* contig length, in first-appearance order.
///
/// Only the processed contigs are emitted (a `--region chr7` run lists just the
/// contigs it covers), so for a whole-file run this is the full reference
/// dictionary and for a subset run it is the subset's contigs.
///
/// Two invariants the callers uphold for this to be correct:
/// * `regions` arrive in coordinate (tid) order — segmentation sorts them — so
///   first-appearance order equals tid order. The header tid assignment must
///   match the order records are emitted, or seqair's single-pass index builder
///   rejects the stream as unsorted. (Lexical order, e.g. via `BTreeSet`, would
///   put `chr10` before `chr2` and break this on GRCh38-style references.)
/// * Every `ChunkRegion` carries the true contig length in `last_position`
///   (the contig end for a whole contig, the BAM header length for a
///   user-defined subset), so we take it directly rather than summing chunk
///   lengths — chunks overlap, and a subset chunk is shorter than its contig.
fn header_contigs(regions: &[ChunkRegion]) -> Vec<Contig> {
    let mut contigs: Vec<Contig> = Vec::new();
    let mut seen: FxHashSet<SmolStr> = FxHashSet::default();
    for region in regions {
        if seen.insert(region.contig.clone()) {
            contigs.push(Contig { name: region.contig.clone(), length: region.last_position });
        }
    }
    contigs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::{ChunkRegion, Region};

    /// `contig_len` is the true contig length carried by every chunk in its
    /// `last_position`; `start`/`end` are the (possibly partial) chunk bounds.
    fn chunk(contig: &str, start: u64, end: u64, contig_len: u64) -> ChunkRegion {
        ChunkRegion {
            region: Region { contig: SmolStr::from(contig), start, end },
            last_position: contig_len,
            overlap_start: 0,
            overlap_end: 0,
        }
    }

    /// The header keeps the order contigs appear in `regions` (which segmentation
    /// has sorted into tid order) rather than lexical order: `chr10` appears
    /// before `chr2`, and emitting against a lexically-sorted header is exactly
    /// what tripped seqair's "input not sorted" index check.
    #[test]
    fn header_contigs_preserve_reference_order_not_lexical() {
        let regions = [
            chunk("chr1", 1, 100, 100),
            chunk("chr2", 1, 100, 100),
            chunk("chr10", 1, 100, 100),
            chunk("chrUn_decoy", 1, 100, 100),
        ];

        let contigs = header_contigs(&regions);
        let names: Vec<&str> = contigs.iter().map(|c| c.name.as_str()).collect();

        assert_eq!(names, ["chr1", "chr2", "chr10", "chrUn_decoy"]);
    }

    /// Multiple (overlapping) chunks of one contig collapse to a single header
    /// entry whose length is the true contig length from `last_position`, not the
    /// sum of chunk lengths — summing would overcount because chunks overlap and a
    /// subset chunk is shorter than its contig.
    #[test]
    fn header_contigs_use_true_contig_length_not_chunk_sum() {
        let regions =
            [chunk("chr1", 1, 100, 250), chunk("chr1", 90, 250, 250), chunk("chr2", 1, 50, 50)];

        let contigs = header_contigs(&regions);

        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[0], Contig { name: SmolStr::from("chr1"), length: 250 });
        assert_eq!(contigs[1], Contig { name: SmolStr::from("chr2"), length: 50 });
    }
}
