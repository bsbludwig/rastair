use rastair_types::Probability;

use crate::{metrics::PileupMetrics, utils::cli};

#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct RecordFilters {
    /// Output all positions, even if they do not pass filters.
    ///
    /// If combined with `--cpgs-only`, only CpG positions will be reported,
    /// including non-passing ones.
    #[arg(long = "all")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub vcf_all: bool,

    /// Report CpGs only and default to BED output
    ///
    /// Only report positions that are CpGs in the reference or variants that
    /// would result in a de-novo CpG.
    ///
    /// Only if combined with `--all`, non-passing CpG positions will also be
    /// reported.
    #[arg(short = 'c', long, default_value_t = false)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub cpgs_only: bool,
}

impl RecordFilters {
    /// Check if a VCF record matches the filter criteria
    pub fn matches(&self, record: &PileupMetrics, ml_threshold: Option<Probability>) -> bool {
        // filter for CpGs, this takes precedence over "all"
        if self.cpgs_only && !(*record.pos_metrics.cpg || record.forms_denovo()) {
            return false;
        }

        // filter for passing records if desired
        if self.vcf_all {
            return true;
        }

        // reject records without alts
        if record.alts.is_empty() {
            return false;
        }

        // okay and now only those that pass
        record.pass(ml_threshold)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        call::pileup::Pileup,
        metrics::Alt,
        sequence::{ChunkRegion, Region},
        vcf::{InCpG, lowDp},
    };
    use rastair_types::Base::*;
    use rastair_vcf::standard_fields::PASS;
    use smallvec::SmallVec;

    use super::*;

    fn default<T: Default>() -> T {
        T::default()
    }

    fn default_record() -> PileupMetrics {
        PileupMetrics {
            pileup: Pileup {
                region: ChunkRegion {
                    region: Region { contig: "chr6".into(), start: 0, end: 1000 },
                    last_position: 10_000,
                },
                context: crate::vcf::SequenceContext {
                    before_2: Some(A),
                    before_1: Some(C),
                    me: C,
                    after_1: Some(G),
                    after_2: Some(A),
                },
                pos: 123,
                reads: crate::call::pileup::SimpleReads(SmallVec::new()),
                reference_base: C,
            },
            pos_metrics: default(),
            pos_filters: default(),
            ref_metrics: default(),
            alts: [Alt { base: T, filters: default(), metrics: default() }].into(),
        }
    }

    #[test]
    fn test_pass() {
        let filters = RecordFilters { vcf_all: false, cpgs_only: false };
        let ml_threshold = Some(Probability::new(0.9).unwrap());

        // nothing going on, filters are empty
        assert!(
            filters.matches(&default_record(), ml_threshold),
            "should match record with not filters"
        );

        // record fails filters
        let mut r = default_record();
        r.pos_filters.add(lowDp, || true);
        assert!(!filters.matches(&r, ml_threshold), "should not match record with failing filters");
    }

    #[test]
    fn test_cpgs() {
        let filters = RecordFilters { vcf_all: false, cpgs_only: true };
        let ml_threshold = Some(Probability::new(0.9).unwrap());

        // empty record is not in CpG
        let mut r = default_record();
        assert!(!filters.matches(&r, ml_threshold), "should not match non-CpG");

        r.pos_filters.add(lowDp, || true);
        assert!(!filters.matches(&r, ml_threshold), "should not match non-CpG failing filters");

        // explicit non-CpG record
        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::No;
        assert!(!filters.matches(&r, ml_threshold), "should not match non-CpG record");

        // CpG record
        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::C;
        assert!(filters.matches(&r, ml_threshold), "should match CpG record");

        r.pos_filters.add(lowDp, || true);
        assert!(!filters.matches(&r, ml_threshold), "should not match CpG record failing filters");

        // denovo CpG candidate
        let mut r = default_record();
        r.alts[0].metrics.denovo = crate::metrics::FormsDenovo::ThisBecomesG;
        assert!(filters.matches(&r, ml_threshold), "should match de-novo CpG candidate record");

        r.pos_filters.add(lowDp, || true);
        assert!(
            !filters.matches(&r, ml_threshold),
            "should not match de-novo CpG candidate record failing filters"
        );
    }

    #[test]
    fn test_all_cpgs() {
        let filters = RecordFilters { vcf_all: true, cpgs_only: true };
        let ml_threshold = Some(Probability::new(0.9).unwrap());

        // empty record is not in CpG
        let mut r = default_record();
        assert!(!filters.matches(&r, ml_threshold), "should not match non-CpG");

        r.pos_filters.add(lowDp, || true);
        assert!(!filters.matches(&r, ml_threshold), "should not match non-CpG failing filters");

        // explicit non-CpG record
        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::No;
        assert!(!filters.matches(&r, ml_threshold), "should not match non-CpG record");

        // CpG record
        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::C;
        assert!(filters.matches(&r, ml_threshold), "should match CpG record");

        r.pos_filters.add(lowDp, || true);
        assert!(filters.matches(&r, ml_threshold), "should match CpG record failing filters");

        // denovo CpG candidate
        let mut r = default_record();
        r.alts[0].metrics.denovo = crate::metrics::FormsDenovo::ThisBecomesG;
        assert!(filters.matches(&r, ml_threshold), "should match de-novo CpG candidate record");

        r.pos_filters.add(lowDp, || true);
        assert!(
            filters.matches(&r, ml_threshold),
            "should match de-novo CpG candidate record failing filters"
        );
    }

    #[test]
    fn test_all() {
        let filters = RecordFilters { vcf_all: true, cpgs_only: false };
        let ml_threshold = Some(Probability::new(0.9).unwrap());

        // empty record
        let mut r = default_record();
        assert!(filters.matches(&r, ml_threshold), "should match record");

        r.pos_filters.add(lowDp, || true);
        assert!(filters.matches(&r, ml_threshold), "should match record failing filters");

        r.pos_filters.add(PASS, || true);
        assert!(filters.matches(&r, ml_threshold), "should match passing record");
    }
}
