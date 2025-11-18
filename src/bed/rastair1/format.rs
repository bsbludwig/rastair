use crate::bed::BedRecord;
use color_eyre::Result;
use rastair_types::{Phred, Probability, SmolStr, Strand, smol_str::format_smolstr};
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
    pub genotype: Genotype,
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
            r#ref,
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

        let genotype = genotype_to_rastair1_string(genotype, r#ref);
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
