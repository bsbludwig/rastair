use crate::metrics::{AltCall, DenovoAdjecent, FormsDenovo, PileupMetrics};
use seqair_types::Base::*;

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
    pub fn from_pileup(pileup: &PileupMetrics) -> Self {
        let ref_base = pileup.ref_base();
        let cpg = *pileup.pos_metrics.cpg;

        if cpg && ref_base == C {
            MethylationEvidenceStrandInfo::from_c(pileup)
        } else if cpg && ref_base == G {
            MethylationEvidenceStrandInfo::from_g(pileup)
        } else if pileup.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingC {
            MethylationEvidenceStrandInfo::from_c(pileup)
        } else if pileup.pos_metrics.denovo_adj == DenovoAdjecent::ThisIsTheMatchingG {
            MethylationEvidenceStrandInfo::from_g(pileup)
        } else if let Some(denovo) =
            pileup.alts.iter().find(|a| *a.metrics.denovo && a.call == AltCall::RealVariant)
        {
            if denovo.metrics.denovo == FormsDenovo::ThisBecomesC {
                MethylationEvidenceStrandInfo::from_c(pileup)
            } else if denovo.metrics.denovo == FormsDenovo::ThisBecomesG {
                MethylationEvidenceStrandInfo::from_g(pileup)
            } else {
                // Should be never happen since we filtered above
                MethylationEvidenceStrandInfo::default()
            }
        } else {
            MethylationEvidenceStrandInfo::default()
        }
    }

    pub fn from_pileup_with_methylation(pileup: &PileupMetrics) -> Self {
        let m = &pileup.pos_metrics.extended.methylated;
        let denovo_only = m.denovo().is_some_and(|b| b.has_evidence())
            && !m.original().is_some_and(|b| b.has_evidence());
        if denovo_only {
            Self::from_denovo_context(pileup).unwrap_or_else(|| Self::from_pileup(pileup))
        } else {
            Self::from_pileup(pileup)
        }
    }

    fn from_c(pileup: &PileupMetrics) -> Self {
        let c = pileup.allele(C).map(|a| a.strand_count).unwrap_or_default();
        let t = pileup.allele(T).map(|a| a.strand_count).unwrap_or_default();

        Self { unmod: c.ot, modified: t.ot, no_snp: c.ob, snp: t.ob }
    }

    fn from_g(pileup: &PileupMetrics) -> Self {
        let g = pileup.allele(G).map(|a| a.strand_count).unwrap_or_default();
        let a = pileup.allele(A).map(|a| a.strand_count).unwrap_or_default();

        Self { unmod: g.ob, modified: a.ob, no_snp: g.ot, snp: a.ot }
    }

    fn from_denovo_context(pileup: &PileupMetrics) -> Option<Self> {
        let context = &pileup.context;
        let has_denovo_c =
            pileup.alts.iter().any(|alt| alt.base == C && alt.call == AltCall::RealVariant);
        if has_denovo_c && context.after_1 == Some(G) {
            return Some(Self::from_c(pileup));
        }

        let has_denovo_g =
            pileup.alts.iter().any(|alt| alt.base == G && alt.call == AltCall::RealVariant);
        if has_denovo_g && context.before_1 == Some(C) {
            return Some(Self::from_g(pileup));
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics;
    use crate::pileups;
    use color_eyre::Result;

    #[test]
    fn denovo_inside_ref_cpg_uses_denovo_context() -> Result<()> {
        let (_segment, pileups) = pileups!(
            [ C C G ] Ref,
            [ C G G ] OB,
            [ C A G ] OB,
            [ C G G ] OB,
        );

        let mut metrics = PileupMetrics::new(pileups[1].clone())?;
        if let Some(alt) = metrics.alts.iter_mut().find(|alt| alt.base == G) {
            alt.call = AltCall::RealVariant;
        }

        metrics.pos_metrics.extended.methylated =
            metrics::methylation::call(&metrics)?.unwrap_or_default();

        assert!(metrics.pos_metrics.extended.methylated.denovo().is_some_and(|b| b.has_evidence()));
        assert!(
            !metrics.pos_metrics.extended.methylated.original().is_some_and(|b| b.has_evidence())
        );

        let info = MethylationEvidenceStrandInfo::from_pileup_with_methylation(&metrics);
        assert_eq!(info.unmod, 2);
        assert_eq!(info.modified, 1);
        assert_eq!(info.no_snp, 0);
        assert_eq!(info.snp, 0);

        Ok(())
    }

    #[test]
    fn matching_c_with_denovo_candidate_prefers_matching_context() -> Result<()> {
        let (_segment, pileups) = pileups!(
            [ C C A ] Ref,
            [ C C A ] OT,
            [ C C A ] OT,
            [ C T A ] OT,
            [ C C A ] OB,
            [ C C A ] OB,
            [ C G A ] OT,
        );

        let mut metrics = PileupMetrics::new(pileups[1].clone())?;
        metrics.pos_metrics.extended.denovo_adj = DenovoAdjecent::ThisIsTheMatchingC;

        let info = MethylationEvidenceStrandInfo::from_pileup(&metrics);
        assert_eq!(info.unmod, 2);
        assert_eq!(info.modified, 1);
        assert_eq!(info.no_snp, 2);
        assert_eq!(info.snp, 0);

        Ok(())
    }
}
