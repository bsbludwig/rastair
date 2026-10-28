use crate::bed::BedRecord;
use crate::call::variant_calling::GenotypeTag;
use color_eyre::Result;
use seqair_types::{Base, Phred, Probability, SmolStr, Strand, smol_str::format_smolstr};
use std::io::Write;

/// The genotype quality column. The format has no missing marker for it, so a
/// record without a genotype (`./.`, no `GQ`) writes 0 on both the `call --bed`
/// and the `convert` path.
pub fn gq_column(gq: Option<Phred>) -> Phred {
    gq.unwrap_or_else(|| Phred::from_phred(0_u8))
}

#[derive(Debug, Clone)]
pub struct Rastair1BedFormat {
    pub contig: SmolStr,
    pub pos: usize,
    pub r#ref: SmolStr,
    pub beta: Option<Probability>,
    /// `None` where the position is neither a CpG nor next to a de-novo one.
    pub strand: Option<Strand>,
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
        let strand = strand.as_ref().map_or(".", Strand::as_symbol);
        let beta = beta.map_or_else(SmolStr::default, |beta| format_smolstr!("{beta:.2}"));

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
pub struct GenotypeString(pub Base, pub Base);

impl GenotypeString {
    /// Creates a [`GenotypeString`] from a [`GenotypeTag`] and the actual bases.
    ///
    /// `alt_bases` should contain the alternative alleles in order (index 0 = first alt, etc.)
    pub fn from_genotype_tag(
        genotype: GenotypeTag,
        ref_base: Base,
        alt_bases: &[Base],
    ) -> GenotypeString {
        match genotype {
            GenotypeTag::HomRef => GenotypeString(ref_base, ref_base),
            GenotypeTag::RefHet(alt_idx) => {
                let alt_base =
                    alt_bases.get(usize::from(alt_idx.get()) - 1).copied().unwrap_or(ref_base);
                GenotypeString(ref_base, alt_base)
            }
            GenotypeTag::HomAlt(alt_idx) => {
                let alt_base =
                    alt_bases.get(usize::from(alt_idx.get()) - 1).copied().unwrap_or(ref_base);
                GenotypeString(alt_base, alt_base)
            }
            GenotypeTag::AltHet(alt1_idx, alt2_idx) => {
                let alt1 =
                    alt_bases.get(usize::from(alt1_idx.get()) - 1).copied().unwrap_or(ref_base);
                let alt2 =
                    alt_bases.get(usize::from(alt2_idx.get()) - 1).copied().unwrap_or(ref_base);
                GenotypeString(alt1, alt2)
            }
        }
    }
}

impl std::fmt::Display for GenotypeString {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}/{}", self.0, self.1)
    }
}

#[cfg(test)]
mod gq_column_tests {
    use super::*;

    #[test]
    fn a_missing_gq_writes_zero_and_a_real_one_is_kept() {
        assert_eq!(gq_column(None).as_int(), 0);
        assert_eq!(gq_column(Some(Phred::from_phred(37_u8))).as_int(), 37);
    }
}
