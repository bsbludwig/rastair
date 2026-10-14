use super::{FeatureCalculator, FeatureNames, FeatureNum, StandardFeatures, define_features};
use crate::{
    call::haplotypes::features::ColumnEvidence,
    metrics::{MetricsForAlt, MetricsForIndel, PileupMetrics},
};
use color_eyre::Result;
use ndarray::Array2;

define_features! {
    pub struct HaplotypeEvidence {
        scalar haplotype_available;
        scalar haplotype_alt_probability;
        scalar haplotype_heterozygous_probability;
        scalar haplotype_homozygous_probability;
        scalar haplotype_expected_fraction;
        scalar haplotype_log_bayes_factor;
        scalar haplotype_genotype_entropy;
        scalar haplotype_context_sites;
        scalar haplotype_linking_fragments;
        scalar haplotype_fragments;
        scalar haplotype_realigned_fraction;
        scalar haplotype_pruned_mass;
        scalar haplotype_boundary_resolved;
    }
}

impl Default for HaplotypeEvidence {
    fn default() -> Self {
        bytemuck::Zeroable::zeroed()
    }
}

impl std::fmt::Debug for HaplotypeEvidence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.as_row()).finish()
    }
}

#[derive(Debug)]
pub struct HaplotypeFeatures;

fn append(base: Result<Array2<f32>>, evidence: HaplotypeEvidence) -> Result<Array2<f32>> {
    let mut values = base?.into_raw_vec_and_offset().0;
    values.extend_from_slice(evidence.as_row());
    super::row_to_array(&values)
}

impl FeatureCalculator for HaplotypeFeatures {
    fn feature_num(&self) -> FeatureNum {
        let mut n = StandardFeatures.feature_num();
        n.cpg += HaplotypeEvidence::FEATURES;
        n.denovo_cpg += HaplotypeEvidence::FEATURES;
        n.others += HaplotypeEvidence::FEATURES;
        n.insertion += HaplotypeEvidence::FEATURES;
        n.deletion += HaplotypeEvidence::FEATURES;
        n
    }

    fn feature_names(&self) -> FeatureNames {
        let mut names = StandardFeatures.feature_names();
        for row in [
            &mut names.cpg,
            &mut names.denovo_cpg,
            &mut names.others,
            &mut names.insertion,
            &mut names.deletion,
        ] {
            row.extend(HaplotypeEvidence::names());
        }
        names
    }

    fn calculate_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f32>> {
        append(StandardFeatures.calculate_cpg(current, before, after), ColumnEvidence::snp(current))
    }
    fn calculate_denovo_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f32>> {
        append(
            StandardFeatures.calculate_denovo_cpg(current, before, after),
            ColumnEvidence::snp(current),
        )
    }
    fn calculate_others(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f32>> {
        append(
            StandardFeatures.calculate_others(current, before, after),
            ColumnEvidence::snp(current),
        )
    }
    fn calculate_insertion(&self, current: &MetricsForIndel) -> Result<Array2<f32>> {
        append(StandardFeatures.calculate_insertion(current), ColumnEvidence::indel(current))
    }
    fn calculate_deletion(&self, current: &MetricsForIndel) -> Result<Array2<f32>> {
        append(StandardFeatures.calculate_deletion(current), ColumnEvidence::indel(current))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn molecular_features_append_to_each_frozen_legacy_layout() {
        let old = StandardFeatures.feature_names();
        let new = HaplotypeFeatures.feature_names();
        let counts = HaplotypeFeatures.feature_num();
        for (old, new, expected) in [
            (old.cpg, new.cpg, counts.cpg),
            (old.denovo_cpg, new.denovo_cpg, counts.denovo_cpg),
            (old.others, new.others, counts.others),
            (old.insertion, new.insertion, counts.insertion),
            (old.deletion, new.deletion, counts.deletion),
        ] {
            assert_eq!(new.get(..old.len()), Some(old.as_slice()));
            assert_eq!(new.get(old.len()..), Some(HaplotypeEvidence::names().as_slice()));
            assert_eq!(new.len(), expected);
            assert_eq!(new.len(), old.len() + 13);
        }
        assert_eq!(HaplotypeEvidence::default().as_row(), &[0.0; 13]);
    }
}
