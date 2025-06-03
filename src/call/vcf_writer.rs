use std::ffi::OsStr;

use clio::ClioPath;
use color_eyre::eyre::{Result, WrapErr};
use rastair2_vcf::{Vcf, VcfBuilder, VcfFormat};
use smallvec::SmallVec;
use smol_str::SmolStr;
use tracing::{debug, warn};

use crate::{call::vcf::Record, sequence::ChunkRegion};

#[derive(Debug, clap::Parser)]
pub struct Params {
    /// VCF/BCF output file path (use - to write to stdout)
    ///
    /// Format is guessed based on the file extension:
    /// `.vcf` for VCF (uncompressed),
    /// `.bcf` for BCF (compressed)
    #[arg(short = 'o', long)]
    pub vcf_output: ClioPath,
}

impl Params {
    /// Create a new instance of `Params` with the specified VCF output path.
    pub fn guess_format(&self) -> VcfFormat {
        match self.vcf_output.extension().map(|s| s.to_ascii_lowercase()) {
            Some(ext) if ext == OsStr::new("bcf") => VcfFormat::Bcf,
            Some(ext) if ext == OsStr::new("vcf") => VcfFormat::Vcf,
            None => VcfFormat::Vcf,
            _ => {
                warn!("Unexpected file extension, defaulting to VCF format without compression.");
                VcfFormat::Vcf
            }
        }
    }

    pub fn vcf_writer(&self, regions: &[ChunkRegion]) -> Result<Vcf<Record>> {
        let format = self.guess_format();
        let compression = match format {
            VcfFormat::Bcf => rastair2_vcf::Compression::On,
            VcfFormat::Vcf => rastair2_vcf::Compression::Off,
        };
        debug!(
            target=?self.vcf_output.display(), ?format, ?compression,
            "Creating VCF writer",
        );
        let writer = VcfBuilder::new(&self.vcf_output, format, compression)
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
