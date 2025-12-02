use crate::metrics::PileupMetrics;
use color_eyre::eyre::{Context, Result};
use cstr8::{CStr8, cstr8};
use rastair_types::Base::*;
use rastair_vcf::{HeaderField, InfoField, InfoFieldNumber, VcfField};
use rust_htslib::bcf;

#[derive(Debug, Default, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct MethylationEvidenceStrandInfo {
    /// Number of unmethylated reads
    ///
    /// - for C: C OT reads
    /// - for G: G OB reads
    pub unmod: u32,
    /// Number of methylated reads
    ///
    /// - for C: T OT reads
    /// - for G: A OB reads
    pub modified: u32,
    /// Number of reads with no SNP evidence
    ///
    /// - for C: C OB reads
    /// - for G: G OT reads
    pub no_snp: u32,
    /// Number of reads with SNP evidence
    ///
    /// - for C: T OB reads
    /// - for G: A OT reads
    pub snp: u32,
}

impl MethylationEvidenceStrandInfo {
    pub fn from_c(pileup: &PileupMetrics) -> Self {
        let c = pileup.allele(C).map(|a| a.strand_count).unwrap_or_default();
        let t = pileup.allele(T).map(|a| a.strand_count).unwrap_or_default();

        Self { unmod: c.ot, modified: t.ot, no_snp: c.ob, snp: t.ob }
    }

    pub fn from_g(pileup: &PileupMetrics) -> Self {
        let g = pileup.allele(G).map(|a| a.strand_count).unwrap_or_default();
        let a = pileup.allele(A).map(|a| a.strand_count).unwrap_or_default();

        Self { unmod: g.ob, modified: a.ob, no_snp: g.ot, snp: a.ot }
    }
}

impl VcfField for MethylationEvidenceStrandInfo {
    const ID: &'static CStr8 = cstr8!("M5cM_Strands");
}

impl HeaderField for MethylationEvidenceStrandInfo {
    const DESCRIPTION: &'static str = "Number of methylated and unmethylated reads supporting each strand, as well as reads with and without SNP evidence";
}

impl InfoField for MethylationEvidenceStrandInfo {
    type Type = u32;
    const NUMBER: InfoFieldNumber = InfoFieldNumber::Num(4);

    #[expect(clippy::cast_possible_wrap, reason = "vcf integer fields with number below i32::MAX")]
    fn write(&self, record: &mut bcf::Record) -> Result<()> {
        record
            .push_info_integer(
                Self::ID,
                &[self.modified as i32, self.unmod as i32, self.no_snp as i32, self.snp as i32],
            )
            .wrap_err("Failed to set M5cM_Strands field")
    }
}
