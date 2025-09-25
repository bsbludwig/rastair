use crate::{
    bed::{BedFormat, BedRecord, writer::BedWriter},
    io::formats::FromFileExtension as _,
    utils::cli,
    vcf::{GenotypeConfidence, GenotypeLikelihood},
};
use clio::ClioPath;
use color_eyre::{Result, eyre::Context as _};
use rastair_types::Probability;
use rastair_vcf::standard_fields::{Genotype, GenotypeAllele};
use smol_str::{SmolStr, format_smolstr};
use std::io::Write;
use tracing::{debug, instrument};

#[derive(Debug, Clone, clap::Args)]
pub struct BedParams {
    /// Output BED file with the called methylated positions
    #[arg(long = "bed", required = false, default_missing_value = "-", num_args = 0..=1)]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub bed: Option<ClioPath>,

    /// Format of the output BED file
    ///
    /// If not specified, the format is guessed based on the file extension.
    #[arg(long, requires = "bed")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub bed_format: Option<BedFormat>,

    #[command(flatten)]
    pub filters: BedRecordsFilterParams,
}

impl BedParams {
    pub fn bed_format(&self) -> BedFormat {
        if let Some(format) = self.bed_format {
            format
        } else if let Some(path) = &self.bed
            && let Some(path) = path.path().to_str()
            && let Some(format) = BedFormat::from_file_extension(path)
        {
            format
        } else {
            debug!(
                "Could not determine BED output format from file extension, defaulting to uncompressed"
            );
            BedFormat::Bed
        }
    }

    #[instrument(level = "debug")]
    pub fn writer(&self) -> Result<Option<BedWriter<Rastair1BedFormat>>> {
        let Some(path) = &self.bed else {
            return Ok(None);
        };

        let format = self.bed_format();
        let writer = BedWriter::new(path, format)
            .wrap_err_with(|| format!("Failed to create BED writer for {path}"))?;
        Ok(Some(writer))
    }
}

#[derive(Debug, Clone)]
pub struct Rastair1BedFormat {
    pub contig: SmolStr,
    pub pos: usize,
    pub r#ref: SmolStr,
    pub beta: Option<Probability>,
    pub unmod: u32,
    pub r#mod: u32,
    pub no_snp: u32,
    pub snp: u32,
    pub coverage: usize,
    pub genotype: Genotype,
    pub genotype_likelihood: GenotypeLikelihood,
    pub genotype_confidence: GenotypeConfidence,
    pub de_novo: bool,
}

/// Parameters for filtering BED records
#[derive(Debug, Clone, clap::Args)]
pub struct BedRecordsFilterParams {
    /// Include CpG positions with zero coverage
    ///
    /// This can be useful to get a complete list of CpG positions in the output BED file.
    /// Note that this requires the input data to contain a complete list of CpG positions,
    /// e.g. by using the `--cpgs-only` option when calling methylation.
    #[arg(long = "bed-include-empty")]
    #[arg(help_heading = cli::sections::FILTER)]
    pub include_empty: bool,
}

// Constructed from command line arguments by the `convert` subcommand, or from
// other existing parameters when used by `call`
#[derive(Debug, Clone, clap::Args)]
pub struct BedRecordsConvertParams {
    /// Minimum ML score to consider a position as variant
    ///
    /// This does nothing if the input data does not contain ML scores.
    #[arg(long = "bed-ml", default_value_t = Probability::new(0.8).expect("valid default probability"))]
    #[arg(help_heading = cli::sections::FILTER)]
    pub ml_threshold: Probability,
    #[command(flatten)]
    pub filters: BedRecordsFilterParams,
}

mod internal_to_bed;
mod vcf_to_bed;

impl BedRecord for Rastair1BedFormat {
    const HEADER: &'static str = "#chr\tstart\tend\tname\tbeta_est\tstrand\tunmod\tmod\tno_snp\tsnp\tcoverage\tgenotype\tgt_p_score\tgt_conf_score\tcpg";

    fn write<W: Write>(&self, f: &mut W) -> Result<()> {
        let Rastair1BedFormat {
            contig,
            pos: start,
            r#ref,
            beta,
            unmod,
            r#mod,
            no_snp,
            snp,
            coverage,
            genotype,
            genotype_likelihood,
            genotype_confidence,
            de_novo,
        } = self;
        let end = start + 1;
        let name = ".";
        let strand = match r#ref.as_str() {
            "C" => "+",
            "G" => "-",
            _ => ".",
        };
        let beta = if let Some(beta) = beta {
            format_smolstr!("{beta:.2}")
        } else {
            SmolStr::new_static("")
        };

        write!(
            f,
            "{contig}\t{start}\t{end}\t{name}\t{beta}\t{strand}\t{unmod}\t{mod}\t{no_snp}\t{snp}\t{coverage}"
        )?;

        let genotype = genotype_to_rastair1_string(genotype, r#ref);
        let likelihood =
            genotype_likelihood.first().and_then(|x| *x).map(|x| x.as_int()).unwrap_or_default();
        let confidence =
            genotype_confidence.first().and_then(|x| *x).map(|x| x.as_int()).unwrap_or_default();
        write!(f, "\t{genotype}\t{likelihood:.2}\t{confidence:.2}")?;
        write!(f, "\t{}", if *de_novo { "NEW" } else { "REF" })?;
        writeln!(f)?;

        Ok(())
    }

    fn chr(&self) -> &str {
        &self.contig
    }

    fn start(&self) -> usize {
        self.pos
    }

    fn end(&self) -> usize {
        self.pos + 1
    }
}

pub fn genotype_to_rastair1_string(genotype: &Genotype, ref_base: &str) -> SmolStr {
    match genotype.0.as_slice() {
        [GenotypeAllele::Phased(0)]
        | [GenotypeAllele::Unphased(0)]
        | [GenotypeAllele::Phased(0), GenotypeAllele::Phased(0)]
        | [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)] => {
            if ref_base == "C" {
                SmolStr::new_static("C/C")
            } else {
                SmolStr::new_static("G/G")
            }
        }
        [GenotypeAllele::Phased(0), GenotypeAllele::Phased(1)]
        | [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)] => {
            if ref_base == "C" {
                SmolStr::new_static("C/T")
            } else {
                SmolStr::new_static("G/A")
            }
        }
        [GenotypeAllele::Phased(1)]
        | [GenotypeAllele::Unphased(1)]
        | [GenotypeAllele::Phased(1), GenotypeAllele::Phased(1)]
        | [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)] => {
            if ref_base == "C" {
                SmolStr::new_static("T/T")
            } else {
                SmolStr::new_static("A/A")
            }
        }
        _ => {
            // would want to return `SmolStr::new_static(".")` but let's be compatible with rastair1 for now
            if ref_base == "C" { SmolStr::new_static("C/C") } else { SmolStr::new_static("G/G") }
        }
    }
}
