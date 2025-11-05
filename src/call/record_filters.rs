use crate::{utils::cli, vcf::Record};

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
    pub fn matches(&self, record: &Record) -> bool {
        // filter for CpGs, this takes precedence over "all"
        if self.cpgs_only && !(*record.info.in_cp_g || *record.info.de_novo_cp_g_candidate) {
            return false;
        }

        // filter for passing records if desired
        if self.vcf_all {
            return true;
        }

        // reject records without alts
        if record.main.alt.is_empty() {
            return false;
        }

        // okay and now only those that pass
        record.filters.pass()
    }
}

#[cfg(test)]
mod tests {
    use crate::vcf::{DeNovoCpGCandidate, InCpG, lowDp};
    use rastair_types::Base::*;
    use rastair_vcf::standard_fields::PASS;
    use smallvec::SmallVec;

    use super::*;

    fn default<T: Default>() -> T {
        T::default()
    }

    fn default_record() -> Record {
        Record {
            main: rastair_vcf::VcfFixedFields {
                chrom: "1".into(),
                pos: 7,
                id: default(),
                r#ref: "A".into(),
                alt: SmallVec::from(["C".into(), "G".into()]),
                qual: Some(50.0),
            },
            filters: default(),
            info: default(),
            samples: default(),
        }
    }

    #[test]
    fn test_pass() {
        let filters = RecordFilters { vcf_all: false, cpgs_only: false };

        // nothing going on, filters are empty
        assert!(filters.matches(&default_record()), "should match record with not filters");

        // record fails filters
        let mut r = default_record();
        r.filters.add_all(lowDp);
        assert!(!filters.matches(&r), "should not match record with failing filters");
    }

    #[test]
    fn test_cpgs() {
        let filters = RecordFilters { vcf_all: false, cpgs_only: true };

        // empty record is not in CpG
        let mut r = default_record();
        assert!(!filters.matches(&r), "should not match non-CpG");

        r.filters.add_all(lowDp);
        assert!(!filters.matches(&r), "should not match non-CpG failing filters");

        // explicit non-CpG record
        let mut r = default_record();
        r.info.in_cp_g = InCpG::No;
        assert!(!filters.matches(&r), "should not match non-CpG record");

        // CpG record
        let mut r = default_record();
        r.info.in_cp_g = InCpG::C;
        assert!(filters.matches(&r), "should match CpG record");

        r.filters.add_all(lowDp);
        assert!(!filters.matches(&r), "should not match CpG record failing filters");

        // denovo CpG candidate
        let mut r = default_record();
        r.info.de_novo_cp_g_candidate =
            DeNovoCpGCandidate::Candidate { ref_base: G, alt_base: A, alt_index: 1 };
        assert!(filters.matches(&r), "should match de-novo CpG candidate record");

        r.filters.add_all(lowDp);
        assert!(
            !filters.matches(&r),
            "should not match de-novo CpG candidate record failing filters"
        );
    }

    #[test]
    fn test_all_cpgs() {
        let filters = RecordFilters { vcf_all: true, cpgs_only: true };

        // empty record is not in CpG
        let mut r = default_record();
        assert!(!filters.matches(&r), "should not match non-CpG");

        r.filters.add_all(lowDp);
        assert!(!filters.matches(&r), "should not match non-CpG failing filters");

        // explicit non-CpG record
        let mut r = default_record();
        r.info.in_cp_g = InCpG::No;
        assert!(!filters.matches(&r), "should not match non-CpG record");

        // CpG record
        let mut r = default_record();
        r.info.in_cp_g = InCpG::C;
        assert!(filters.matches(&r), "should match CpG record");

        r.filters.add_all(lowDp);
        assert!(filters.matches(&r), "should match CpG record failing filters");

        // denovo CpG candidate
        let mut r = default_record();
        r.info.de_novo_cp_g_candidate =
            DeNovoCpGCandidate::Candidate { ref_base: G, alt_base: A, alt_index: 1 };
        assert!(filters.matches(&r), "should match de-novo CpG candidate record");

        r.filters.add_all(lowDp);
        assert!(filters.matches(&r), "should match de-novo CpG candidate record failing filters");
    }

    #[test]
    fn test_all() {
        let filters = RecordFilters { vcf_all: true, cpgs_only: false };

        // empty record
        let mut r = default_record();
        assert!(filters.matches(&r), "should match record");

        r.filters.add_all(lowDp);
        assert!(filters.matches(&r), "should match record failing filters");

        r.filters.add_all(PASS);
        assert!(filters.matches(&r), "should match passing record");
    }
}
