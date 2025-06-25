use clio::ClioPath;
use color_eyre::eyre::{Result, WrapErr};
use rastair2_vcf::{Compression, Vcf, VcfBuilder, VcfFormat};
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::{ffi::OsStr, num::NonZeroUsize};
use tracing::{debug, warn};

use crate::{sequence::ChunkRegion, vcf::Record};

#[derive(Debug, clap::Parser)]
pub struct Params {
    /// VCF/BCF output file path (use - to write to stdout)
    ///
    /// Format is guessed based on the file extension:
    /// `.vcf` for VCF (uncompressed),
    /// `.vcf.gz` for VCF (compressed),
    /// `.bcf` for BCF (compressed)
    #[arg(short = 'o', long, default_value = "-")]
    pub vcf_output: ClioPath,

    /// Number of threads to use for writing (and compressing) VCF files
    ///
    /// This is subtracted from `--threads` but never below 1
    // Default value chosen after profiling on a machine with 14 cores.
    #[arg(long, default_value = "2")]
    pub vcf_threads: NonZeroUsize,
}

impl Params {
    /// Create a new instance of `Params` with the specified VCF output path.
    pub fn guess_format(&self) -> (VcfFormat, Compression) {
        if self.vcf_output == ClioPath::new("-").expect("static path is valid") {
            return (VcfFormat::Vcf, Compression::Off);
        }

        let Some(name) = self.vcf_output.file_name().and_then(OsStr::to_str) else {
            warn!(
                "No file name found in VCF output path, defaulting to VCF format without compression."
            );
            return (VcfFormat::Vcf, Compression::Off);
        };

        if name.ends_with("bcf") {
            (VcfFormat::Bcf, Compression::On)
        } else if name.ends_with("vcf") {
            (VcfFormat::Vcf, Compression::Off)
        } else if name.ends_with("vcf.gz") {
            (VcfFormat::Vcf, Compression::On)
        } else {
            warn!("Unexpected file extension, defaulting to VCF format without compression.");
            (VcfFormat::Vcf, Compression::Off)
        }
    }

    pub fn vcf_writer(&self, regions: &[ChunkRegion]) -> Result<Vcf<Record>> {
        let (format, compression) = self.guess_format();
        debug!(
            target=?self.vcf_output.display(), ?format, ?compression,
            "Creating VCF writer",
        );
        let writer = VcfBuilder::new(&self.vcf_output, format, compression, self.vcf_threads.get())
            .wrap_err("Failed to create VCF writer")?;

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
}
