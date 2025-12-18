use crate::bed::BedRecord;
use color_eyre::Result;
use rastair_types::{Base, Phred, Probability, SmolStr, Strand, smol_str::format_smolstr};
use rastair_vcf::standard_fields::{Genotype, GenotypeAllele};
use std::io::Write;
use tracing::debug;

#[derive(Debug, Clone)]
pub struct Rastair1BedFormat {
    pub contig: SmolStr,
    pub pos: usize,
    pub r#ref: SmolStr,
    pub beta: Option<Probability>,
    pub strand: Strand,
    pub unmod: u32,
    pub r#mod: u32,
    pub no_snp: u32,
    pub snp: u32,
    pub coverage: usize,
    pub genotype: GenotypeString,
    pub genotype_likelihood: Phred,
    pub genotype_confidence: Phred,
    pub de_novo: bool,
}

impl BedRecord for Rastair1BedFormat {
    const HEADER: &'static str = "#chr\tstart\tend\tname\tbeta_est\tstrand\tunmod\tmod\tno_snp\tsnp\tcoverage\tgenotype\tgt_p_score\tgt_conf_score\tcpg";

    fn write<W: Write>(&self, f: &mut W) -> Result<()> {
        let Rastair1BedFormat {
            contig,
            pos: start,
            r#ref: _,
            beta,
            strand,
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
        let strand = strand.as_symbol();
        let beta = if let Some(beta) = beta {
            format_smolstr!("{beta:.2}")
        } else {
            debug!("position {contig}:{start} has no beta value");
            SmolStr::new_inline(".")
        };

        write!(
            f,
            "{contig}\t{start}\t{end}\t{name}\t{beta}\t{strand}\t{unmod}\t{mod}\t{no_snp}\t{snp}\t{coverage}"
        )?;

        // let genotype = genotype_to_rastair1_string(genotype, r#ref);
        let genotype_likelihood = genotype_likelihood.as_int();
        let genotype_confidence = genotype_confidence.as_int();
        write!(f, "\t{genotype}\t{genotype_likelihood}\t{genotype_confidence}")?;
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenotypeString {
    CC,
    CT,
    TT,
    GG,
    GA,
    AA,
    Unknown,
}

impl GenotypeString {
    pub fn as_str(&self) -> &str {
        match self {
            GenotypeString::CC => "C/C",
            GenotypeString::CT => "C/T",
            GenotypeString::TT => "T/T",
            GenotypeString::GG => "G/G",
            GenotypeString::GA => "G/A",
            GenotypeString::AA => "A/A",
            GenotypeString::Unknown => ".",
        }
    }

    pub fn from_genotype(genotype: &Genotype, ref_base: Base) -> GenotypeString {
        use rastair_types::Base::*;

        match genotype.0.as_slice() {
            [GenotypeAllele::Phased(0)]
            | [GenotypeAllele::Unphased(0)]
            | [GenotypeAllele::Phased(0), GenotypeAllele::Phased(0)]
            | [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)] => {
                if ref_base == C {
                    GenotypeString::CC
                } else {
                    GenotypeString::GG
                }
            }
            [GenotypeAllele::Phased(0), GenotypeAllele::Phased(1)]
            | [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)] => {
                if ref_base == C {
                    GenotypeString::CT
                } else {
                    GenotypeString::GA
                }
            }
            [GenotypeAllele::Phased(1)]
            | [GenotypeAllele::Unphased(1)]
            | [GenotypeAllele::Phased(1), GenotypeAllele::Phased(1)]
            | [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)] => {
                if ref_base == C {
                    GenotypeString::TT
                } else {
                    GenotypeString::AA
                }
            }
            _ => {
                // would want to return `SmolStr::new_static(".")` but let's be compatible with rastair1 for now
                if ref_base == C { GenotypeString::CC } else { GenotypeString::GG }
            }
        }
    }
}

impl std::fmt::Display for GenotypeString {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}
