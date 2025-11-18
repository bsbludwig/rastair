use crate::{
    io::{
        formats::FromFileExtension,
        mpk::{MessagePackWriter, format::MpkVcfHeader},
    },
    sequence::ChunkRegion,
    utils::{cli, logging::ThisIsABug as _},
    vcf::Record,
};
use better_default::Default;
use clio::ClioPath;
use color_eyre::eyre::{ContextCompat, Result, WrapErr};
use rastair_types::SmolStr;
use rastair_vcf::{Compression, Contig, VcfBuilder, VcfFile, VcfFormat as HtsVcfFormat};
use std::{collections::BTreeSet, ffi::OsStr, num::NonZeroUsize};
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
    #[arg(long, default_value = "2")]
    #[arg(help_heading = cli::sections::PROCESSING)]
    #[default(NonZeroUsize::new(1).expect("3 > 0"))]
    pub vcf_threads: NonZeroUsize,
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

impl From<VcfFormat> for (HtsVcfFormat, Compression) {
    fn from(format: VcfFormat) -> Self {
        match format {
            VcfFormat::Vcf => (HtsVcfFormat::Vcf, Compression::Off),
            VcfFormat::VcfCompressed => (HtsVcfFormat::Vcf, Compression::On),
            VcfFormat::Bcf => (HtsVcfFormat::Bcf, Compression::On),
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

    pub fn writer(&self, regions: &[ChunkRegion], metadata: &[String]) -> Result<Option<Writer>> {
        let Some(_) = &self.vcf else {
            return Ok(None);
        };

        let contigs: BTreeSet<Contig> =
            regions.iter().map(|r| Contig { name: r.contig.clone(), length: r.len() }).collect();
        let contigs: Vec<Contig> = contigs.into_iter().collect();
        let samples = vec![SmolStr::new("sample")]; // Note: we only deal with one sample for now

        let (format, compression) = match self.guess_format() {
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
            self.vcf_writer(&contigs, &samples, metadata, format, compression)
                .wrap_err("Failed to create VCF writer")?
                .wrap_err("No VCF output path present")
                .this_is_a_bug()?,
        )))
    }

    pub fn vcf_writer(
        &self,
        contigs: &[Contig],
        samples: &[SmolStr],
        metadata: &[String],
        format: HtsVcfFormat,
        compression: Compression,
    ) -> Result<Option<VcfFile<Record>>> {
        let Some(vcf_output) = &self.vcf else {
            return Ok(None);
        };

        debug!(
            target=?vcf_output.display(), ?format, ?compression,
            "Creating VCF writer",
        );
        let mut writer = VcfBuilder::new(vcf_output, format, compression, self.vcf_threads.get())
            .wrap_err("Failed to create VCF writer")?;

        for line in metadata {
            writer.add_header_line(format!("##{line}"));
        }

        Some(writer.build(contigs, samples).wrap_err("Failed to build VCF writer")).transpose()
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

pub enum Writer {
    Vcf(VcfFile<Record>),
    MessagePack(MessagePackWriter),
}
