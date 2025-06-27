use clio::ClioPath;
use color_eyre::eyre::{Result, WrapErr};
use rastair2_vcf::{Compression, VcfBuilder, VcfFile, VcfFormat as HtsVcfFormat};
use rustc_hash::FxHashSet;
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::{ffi::OsStr, num::NonZeroUsize};
use tracing::{debug, warn};

use crate::{
    io::{
        formats::FromFileExtension,
        mpk::{MessagePackWriter, MpkVcfHeader},
    },
    sequence::ChunkRegion,
    vcf::Record,
};

#[derive(Debug, Clone, clap::Parser)]
pub struct Params {
    /// VCF/BCF output file path (use - to write to stdout)
    ///
    /// Format is guessed based on the file extension:
    /// `.vcf` for VCF (uncompressed),
    /// `.vcf.gz` for VCF (compressed),
    /// `.bcf` for BCF (compressed)
    /// `.mpk.lz4` for internal format (Message Pack, LZ4-compressed)
    #[arg(short = 'o', long, default_value = "-")]
    pub vcf_output: ClioPath,

    /// Number of threads to use for writing (and compressing) VCF files
    ///
    /// This is subtracted from `--threads` but never below 1
    // Default value chosen after profiling on a machine with 14 cores.
    #[arg(long, default_value = "2")]
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
    /// Text-based VCF format
    Vcf,
    /// Compressed text-based VCF format
    VcfCompressed,
    /// Binary VCF format (BCF)
    Bcf,
}

impl Params {
    /// Create a new instance of `Params` with the specified VCF output path.
    pub fn guess_format(&self) -> Format {
        if self.vcf_output.is_std() {
            return Format::Vcf(VcfFormat::Vcf);
        }

        let Some(name) = self.vcf_output.file_name().and_then(OsStr::to_str) else {
            warn!(
                filename=%self.vcf_output,
                "No file name found in VCF output path, defaulting to VCF format without compression."
            );
            return Format::Vcf(VcfFormat::Vcf);
        };

        let Some(format) = Format::from_file_extension(name) else {
            warn!(
                filename=%self.vcf_output,
                "Could not determine format from file extension, defaulting to VCF format without compression."
            );
            return Format::Vcf(VcfFormat::Vcf);
        };

        format
    }

    pub fn writer(&self, regions: &[ChunkRegion], metadata: &[String]) -> Result<Writer> {
        let (format, compression) = match self.guess_format() {
            Format::MessagePack => {
                let writer = self
                    .create_mpk_writer(regions, metadata)
                    .wrap_err("Failed to create MessagePack writer")?;
                return Ok(Writer::MessagePack(writer));
            }
            Format::Vcf(VcfFormat::Vcf) => (HtsVcfFormat::Vcf, Compression::Off),
            Format::Vcf(VcfFormat::VcfCompressed) => (HtsVcfFormat::Vcf, Compression::On),
            Format::Vcf(VcfFormat::Bcf) => (HtsVcfFormat::Bcf, Compression::On),
        };

        Ok(Writer::Vcf(
            self.vcf_writer(regions, metadata, format, compression)
                .wrap_err("Failed to create VCF writer")?,
        ))
    }

    pub fn vcf_writer(
        &self,
        regions: &[ChunkRegion],
        metadata: &[String],
        format: HtsVcfFormat,
        compression: Compression,
    ) -> Result<VcfFile<Record>> {
        debug!(
            target=?self.vcf_output.display(), ?format, ?compression,
            "Creating VCF writer",
        );
        let mut writer =
            VcfBuilder::new(&self.vcf_output, format, compression, self.vcf_threads.get())
                .wrap_err("Failed to create VCF writer")?;

        for line in metadata {
            writer.add_header_line(format!("##{line}"));
        }

        // List all chromosomes in the regions for VCF header
        let contigs = regions.iter().map(|r| r.chromosome.clone()).fold(
            SmallVec::<SmolStr, 20>::new(),
            |mut acc, chrom| {
                if !acc.contains(&chrom) {
                    acc.push(chrom);
                }
                acc
            },
        );
        let samples = [SmolStr::new("sample")]; // TODO: we have one sample for now

        writer.build(&contigs, &samples).wrap_err("Failed to build VCF writer")
    }

    pub fn create_mpk_writer(
        &self,
        regions: &[ChunkRegion],
        metadata: &[String],
    ) -> Result<MessagePackWriter> {
        let path = &self.vcf_output;
        warn!(
            %path,
            "MessagePack format only for internal use, no stability guarantees",
        );
        let mut w = MessagePackWriter::new(path)
            .wrap_err_with(|| format!("Failed to create MessagePack writer for {path}"))?;

        // List all chromosomes in the regions for VCF header
        let contigs = regions
            .iter()
            .map(|r| r.chromosome.clone())
            .collect::<FxHashSet<_>>()
            .into_iter()
            .collect();
        let samples = vec![SmolStr::new("sample")]; // TODO: we have one sample for now

        w.add_metadata(MpkVcfHeader { contigs, samples, metadata: metadata.to_owned() })?;

        Ok(w)
    }
}

pub enum Writer {
    Vcf(VcfFile<Record>),
    MessagePack(MessagePackWriter),
}

impl Writer {
    pub fn add(&mut self, record: &Record) -> Result<()> {
        match self {
            Writer::Vcf(writer) => writer.add(record).wrap_err("Failed to write record to VCF"),
            Writer::MessagePack(writer) => writer.add(record),
        }
    }
}
