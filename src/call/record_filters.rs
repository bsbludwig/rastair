use crate::{metrics::PileupMetrics, utils::cli};
use rastair_types::Probability;

/// Filters to apply when deciding whether to output a VCF record
///
/// Here's a table with all the combinations:
///
/// | `vcf_all` | `cpgs_only` | Output behavior                                  |
/// | --------- | ----------- | ------------------------------------------------ |
/// | ``        | ``          | All CpG sites and variants that PASS             |
/// | ``        | `-c`        | All CpG sites and PASSing de-novo CpG candidates |
/// | `--all`   | ``          | All CpG sites and all variants                   |
/// | `--all`   | `-c`        | All CpG and de-novo CpG candidates               |
///
/// Note: We alwasy report both positions of a CpG.
#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize)]
pub struct RecordFilters {
    /// Output all positions, even if they do not pass filters.
    ///
    /// If combined with `--cpgs-only`, only CpG positions will be reported,
    /// including non-passing de-novo CpGs.
    #[arg(long = "all")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub vcf_all: bool,

    /// Report CpGs only and default to BED output
    ///
    /// Only report positions that are CpGs in the reference or variants that
    /// would result in a de-novo CpG.
    ///
    /// Only if combined with `--all`, non-passing de-novo CpG positions will
    /// also be reported.
    #[arg(short = 'c', long, default_value_t = false)]
    #[arg(help_heading = cli::sections::FILTER)]
    pub cpgs_only: bool,
}

impl RecordFilters {
    /// Check if a VCF record matches the filter criteria
    pub fn matches(&self, record: &PileupMetrics, ml_threshold: Option<Probability>) -> bool {
        match (self.vcf_all, self.cpgs_only) {
            (false, false) => {
                // default behavior: only passing records with alts
                *record.pos_metrics.cpg || record.pass(ml_threshold)
            }
            (false, true) => {
                // passing CpGs and passing de-novo CpG candidates
                if *record.pos_metrics.cpg {
                    // - we're a CpG
                    return true;
                }
                if record.forms_denovo() && record.pass(ml_threshold) {
                    // - we're a passing de-novo CpG candidate
                    return true;
                }
                if record.alts.iter().any(|alt| {
                    *alt.metrics.denovo && alt.filters.filters.other_pos_in_denovo_passes
                }) {
                    // - other position passes
                    return true;
                }
                false
            }
            (true, false) => {
                // `--all` means all
                true
            }
            (true, true) => {
                // `--all -c` means all CpGs and de-novo CpG candidates
                *record.pos_metrics.cpg || record.forms_denovo()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        call::pileup::Pileup,
        metrics::{Alt, AltCall, PositionMetrics},
        sequence::{ChunkRegion, Region},
        utils::default,
        vcf::{InCpG, lowDp},
    };
    use rastair_types::Base::*;
    use rastair_types::SmallVec;
    use rastair_vcf::standard_fields::PASS;

    const ML_THRESHOLD: Option<Probability> = Some(Probability::new_panicky(0.9));

    fn default_record() -> PileupMetrics {
        PileupMetrics {
            pileup: Pileup {
                region: ChunkRegion {
                    region: Region { contig: "chr6".into(), start: 0, end: 1000 },
                    last_position: 10_000,
                    overlap_start: 0,
                    overlap_end: 0,
                },
                context: crate::vcf::SequenceContext {
                    before_2: Some(A),
                    before_1: Some(C),
                    me: A,
                    after_1: Some(G),
                    after_2: Some(A),
                },
                pos: 123,
                reads: crate::call::pileup::SimpleReads(SmallVec::new()),
                reference_base: A,
            },
            pos_metrics: PositionMetrics::default(),
            pos_filters: default(),
            ref_metrics: default(),
            alts: [Alt {
                base: C,
                filters: default(),
                metrics: crate::metrics::AlleleMetrics {
                    denovo: crate::metrics::FormsDenovo::No,
                    ..default()
                },
                call: default(),
            }]
            .into(),
        }
    }

    #[test]
    fn test_pass() {
        let filters = RecordFilters { vcf_all: false, cpgs_only: false };

        // nothing going on, filters are empty
        assert!(
            filters.matches(&default_record(), ML_THRESHOLD),
            "should match record with not filters"
        );

        // record fails filters
        let mut r = default_record();
        r.pos_filters.add(lowDp, || true);
        assert!(!filters.matches(&r, ML_THRESHOLD), "should not match record with failing filters");
    }

    #[test]
    fn test_cpgs() {
        let filters = RecordFilters { vcf_all: false, cpgs_only: true };

        // explicit non-CpG record
        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::No;
        assert!(!filters.matches(&r, ML_THRESHOLD), "should not match non-CpG");

        r.pos_filters.add(lowDp, || true);
        assert!(!filters.matches(&r, ML_THRESHOLD), "should not match non-CpG failing filters");

        // CpG record
        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::C;
        assert!(filters.matches(&r, ML_THRESHOLD), "should match CpG record");

        r.pos_filters.add(lowDp, || true);
        assert!(filters.matches(&r, ML_THRESHOLD), "should match even failing CpG record");

        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::C;
        r.pos_filters.add(lowDp, || true);
        r.pos_filters.other_pos_in_denovo_passes = true; // doesn't really change anything
        assert!(
            filters.matches(&r, ML_THRESHOLD),
            "should match CpG record if other position passes"
        );

        // denovo CpG candidate
        let mut r = default_record();
        r.alts[0].metrics.denovo = crate::metrics::FormsDenovo::ThisBecomesG;
        assert!(filters.matches(&r, ML_THRESHOLD), "should match de-novo CpG candidate");

        r.pos_filters.add(lowDp, || true);
        assert!(
            !filters.matches(&r, ML_THRESHOLD),
            "should not match de-novo CpG candidate failing pos filters"
        );

        let mut r = default_record();
        r.alts[0].metrics.denovo = crate::metrics::FormsDenovo::ThisBecomesG;
        r.alts[0].filters.filters.add(lowDp, || true);
        assert!(
            !filters.matches(&r, ML_THRESHOLD),
            "should not match de-novo CpG candidate failing alt filters"
        );

        let mut r = default_record();
        r.alts[0].metrics.denovo = crate::metrics::FormsDenovo::ThisBecomesG;
        r.alts[0].filters.filters.add(lowDp, || true);
        r.alts[0].filters.filters.other_pos_in_denovo_passes = true;
        assert!(
            filters.matches(&r, ML_THRESHOLD),
            "should match de-novo CpG candidate failing filters if other position passes"
        );
    }

    #[test]
    fn test_all_cpgs() {
        let filters = RecordFilters { vcf_all: true, cpgs_only: true };

        // empty record is not in CpG
        let mut r = default_record();
        assert!(!filters.matches(&r, ML_THRESHOLD), "should not match non-CpG");

        r.pos_filters.add(lowDp, || true);
        assert!(!filters.matches(&r, ML_THRESHOLD), "should not match non-CpG failing filters");

        // explicit non-CpG record
        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::No;
        assert!(!filters.matches(&r, ML_THRESHOLD), "should not match non-CpG record");

        // CpG record
        let mut r = default_record();
        r.pos_metrics.cpg = InCpG::C;
        assert!(filters.matches(&r, ML_THRESHOLD), "should match CpG record");

        r.pos_filters.add(lowDp, || true);
        assert!(filters.matches(&r, ML_THRESHOLD), "should match CpG record failing filters");

        // denovo CpG candidate
        let mut r = default_record();
        r.alts[0].metrics.denovo = crate::metrics::FormsDenovo::ThisBecomesG;
        assert!(filters.matches(&r, ML_THRESHOLD), "should match de-novo CpG candidate");

        r.pos_filters.add(lowDp, || true);
        assert!(
            filters.matches(&r, ML_THRESHOLD),
            "should match de-novo CpG candidate failing filters"
        );
    }

    #[test]
    fn test_all() {
        let filters = RecordFilters { vcf_all: true, cpgs_only: false };

        // empty record
        let mut r = default_record();
        assert!(filters.matches(&r, ML_THRESHOLD), "should match record");

        r.pos_filters.add(lowDp, || true);
        assert!(filters.matches(&r, ML_THRESHOLD), "should match record failing filters");

        r.pos_filters.add(PASS, || true);
        assert!(filters.matches(&r, ML_THRESHOLD), "should match passing record");
    }
}
