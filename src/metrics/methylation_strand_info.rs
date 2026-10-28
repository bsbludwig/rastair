use crate::metrics::{
    PileupMetrics,
    methylation::{self, CpgReads},
};
use crate::vcf::CpgSide;

#[derive(Debug, Default, Clone, Copy)]
pub struct MethylationEvidenceStrandInfo {
    /// Number of unmethylated reads, the ones β counts
    ///
    /// - for C: C OT reads followed by a G
    /// - for G: G OB reads preceded by a C
    pub unmod: u32,
    /// Number of methylated reads, the ones β counts
    ///
    /// - for C: T OT reads followed by a G
    /// - for G: A OB reads preceded by a C
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
    /// The reads of [`methylation::evidence`]; none outside a CpG.
    pub fn from_pileup(pileup: &PileupMetrics) -> Self {
        methylation::evidence(pileup)
            .map(|(side, _)| Self::from_side(pileup, side))
            .unwrap_or_default()
    }

    /// TAPS converts `side` on one strand, so there its unmodified and
    /// modified bases are methylation, counted as β counts them; the other
    /// strand reads the reference base or a variant.
    fn from_side(pileup: &PileupMetrics, side: CpgSide) -> Self {
        let CpgReads { modified, unmodified } = methylation::read_counts(pileup, side);
        let plain = |base| {
            let reads = pileup.allele(base).map(|a| a.strand_count).unwrap_or_default();
            side.converting_and_plain(reads).plain
        };
        Self {
            unmod: unmodified,
            modified,
            no_snp: plain(side.unmod_base()),
            snp: plain(side.mod_base()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics;
    use crate::pileups;
    use color_eyre::Result;
    use seqair_types::Base::*;

    #[test]
    fn denovo_inside_ref_cpg_uses_denovo_context() -> Result<()> {
        let (segment, pileups) = pileups!(
            [ C C G ] Ref,
            [ C G G ] OB,
            [ C A G ] OB,
            [ C G G ] OB,
        );

        let mut metrics = pileups[1].metrics(&segment)?;
        if let Some(alt) = metrics.alts.iter_mut().find(|alt| alt.base == G) {
            alt.verdict = crate::metrics::AltVerdict::Called(seqair_types::Probability::ONE);
        }

        metrics.pos_metrics.methylated = metrics::methylation::call(&metrics)?;

        assert!(metrics.pos_metrics.methylated.denovo.is_some());
        assert!(metrics.pos_metrics.methylated.original.is_none());

        let info = MethylationEvidenceStrandInfo::from_pileup(&metrics);
        assert_eq!(info.unmod, 2);
        assert_eq!(info.modified, 1);
        assert_eq!(info.no_snp, 0);
        assert_eq!(info.snp, 0);

        Ok(())
    }

    #[test]
    fn matching_c_with_denovo_candidate_prefers_matching_context() -> Result<()> {
        let (segment, pileups) = pileups!(
            [ C C A ] Ref,
            [ C C G ] OT,
            [ C C G ] OT,
            [ C T G ] OT,
            [ C C A ] OB,
            [ C C A ] OB,
            [ C G A ] OT,
        );

        let mut metrics = pileups[1].metrics(&segment)?;
        assert_eq!(
            MethylationEvidenceStrandInfo::from_pileup(&metrics).unmod,
            0,
            "a G candidate after it whose alt is not called makes no CpG"
        );
        metrics.beside_called_alts = enumset::EnumSet::only(CpgSide::C);

        let info = MethylationEvidenceStrandInfo::from_pileup(&metrics);
        assert_eq!(info.unmod, 2);
        assert_eq!(info.modified, 1);
        assert_eq!(info.no_snp, 2);
        assert_eq!(info.snp, 0);

        Ok(())
    }
}
