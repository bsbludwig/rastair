use clio::ClioPath;
use color_eyre::eyre::{Result, WrapErr};
use rastair2_vcf::{Compression, VcfBuilder, VcfFile, VcfFormat};
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::{ffi::OsStr, fs::File, io::BufWriter, num::NonZeroUsize};
use tracing::{debug, warn};

use crate::{sequence::ChunkRegion, vcf::Record};

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

#[derive(Debug, Clone, Copy)]
pub enum Format {
    /// Text-based VCF format
    Vcf,
    /// Binary VCF format (BCF)
    Bcf,
    /// Rastair-internal `MessagePack` format, always LZ4-compressed
    MessagePack,
}

impl Params {
    /// Create a new instance of `Params` with the specified VCF output path.
    pub fn guess_format(&self) -> (Format, Compression) {
        if self.vcf_output == ClioPath::new("-").expect("static path is valid") {
            return (Format::Vcf, Compression::Off);
        }

        let Some(name) = self.vcf_output.file_name().and_then(OsStr::to_str) else {
            warn!(
                "No file name found in VCF output path, defaulting to VCF format without compression."
            );
            return (Format::Vcf, Compression::Off);
        };

        if name.ends_with("bcf") {
            (Format::Bcf, Compression::On)
        } else if name.ends_with("vcf") {
            (Format::Vcf, Compression::Off)
        } else if name.ends_with("vcf.gz") {
            (Format::Vcf, Compression::On)
        } else if name.ends_with("mpk.lz4") {
            (Format::MessagePack, Compression::On)
        } else {
            warn!("Unexpected file extension, defaulting to VCF format without compression.");
            (Format::Vcf, Compression::Off)
        }
    }

    pub fn writer(&self, regions: &[ChunkRegion], metadata: &[String]) -> Result<Writer> {
        let (format, compression) = self.guess_format();

        let format = match format {
            Format::MessagePack => {
                let writer =
                    self.create_direct_writer().wrap_err("Failed to create MessagePack writer")?;
                return Ok(Writer::MessagePack(writer));
            }
            Format::Vcf => VcfFormat::Vcf,
            Format::Bcf => VcfFormat::Bcf,
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
        format: VcfFormat,
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

    pub fn create_direct_writer(&self) -> Result<DirectWriter> {
        let path = &self.vcf_output;
        warn!(
            %path,
            "MessagePack format only for internal use, no stability guarantees",
        );
        let file = BufWriter::new(
            File::create(path.path()).wrap_err_with(|| format!("Failed to create file {path}"))?,
        );

        lz4::EncoderBuilder::new().level(0).build(file).wrap_err("Failed to create LZ4 encoder")
    }
}

type DirectWriter = lz4::Encoder<BufWriter<File>>;

pub enum Writer {
    Vcf(VcfFile<Record>),
    MessagePack(DirectWriter),
}

impl Writer {
    pub fn add(&mut self, record: &Record) -> Result<()> {
        match self {
            Writer::Vcf(writer) => writer.add(record).wrap_err("Failed to write record to VCF"),
            Writer::MessagePack(writer) => rmp_serde::encode::write_named(writer, record)
                .wrap_err("Failed to write record to MessagePack file"),
        }
    }
}
