//! Feature calculation implementations for ML models
//!
//! Trait-based abstraction for calculating features from variant metrics.
//! Different implementations can be swapped to support various feature sets and
//! model variants.

use super::types::MlFeatureSet;
use crate::metrics::{MetricsForAlt, PileupMetrics};
use color_eyre::{Result, eyre::Context as _};
use ndarray::Array2;
use std::fmt;

pub mod shared;
pub mod standard;
pub mod utils;

#[derive(Debug, Clone, Copy)]
pub struct FeatureNum {
    pub cpg: usize,
    pub denovo_cpg: usize,
    pub others: usize,
}

/// Calculate ML features from variant metrics
pub trait FeatureCalculator: fmt::Debug + Send + Sync {
    fn feature_num(&self) -> FeatureNum;

    /// Calculate features for a CpG methylation candidate
    fn calculate_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>>;

    /// Calculate features for a denovo CpG candidate
    fn calculate_denovo_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>>;

    /// Calculate features for other variants
    fn calculate_others(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>>;
}

impl MlFeatureSet {
    pub fn get_calculator(&self) -> Box<dyn FeatureCalculator> {
        match self {
            MlFeatureSet::Standard => Box::new(StandardFeatures),
            MlFeatureSet::Simple => Box::new(SimpleFeatures),
        }
    }
}

/// Standard implementation of feature calculation using all features
#[derive(Debug, Clone, Copy)]
pub struct StandardFeatures;

impl FeatureCalculator for StandardFeatures {
    fn feature_num(&self) -> FeatureNum {
        FeatureNum { cpg: 56, denovo_cpg: 56, others: 48 }
    }

    fn calculate_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        standard::cpg(current, before, after)
    }

    fn calculate_denovo_cpg(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        standard::denovo_cpg(current, before, after)
    }

    fn calculate_others(
        &self,
        current: &MetricsForAlt,
        before: Option<&PileupMetrics>,
        after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        standard::others(current, before, after)
    }
}

/// Very basic feature calculation using small subset of features
#[derive(Debug, Clone, Copy)]
pub struct SimpleFeatures;

impl SimpleFeatures {
    fn calculate_basic(&self, current: &MetricsForAlt) -> Result<Array2<f64>> {
        let common = shared::CommonFeatures::extract(current);
        let mut features = Vec::with_capacity(48);
        features.extend_from_slice(&common.base_encoding);
        features.extend_from_slice(&common.position_metrics);
        features.extend_from_slice(&common.context_encoding);
        features.extend_from_slice(&common.depth_ratios);
        features.extend_from_slice(&common.base_quality_metrics);
        features.extend_from_slice(&common.mapping_quality_metrics);
        features.extend_from_slice(&common.read_metrics);
        Array2::from_shape_vec((1, features.len()), features)
            .wrap_err("Failed to create basic feature array")
    }
}

impl FeatureCalculator for SimpleFeatures {
    fn feature_num(&self) -> FeatureNum {
        FeatureNum { cpg: 48, denovo_cpg: 48, others: 48 }
    }

    fn calculate_cpg(
        &self,
        current: &MetricsForAlt,
        _before: Option<&PileupMetrics>,
        _after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        self.calculate_basic(current)
    }

    fn calculate_denovo_cpg(
        &self,
        current: &MetricsForAlt,
        _before: Option<&PileupMetrics>,
        _after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        self.calculate_basic(current)
    }

    fn calculate_others(
        &self,
        current: &MetricsForAlt,
        _before: Option<&PileupMetrics>,
        _after: Option<&PileupMetrics>,
    ) -> Result<Array2<f64>> {
        self.calculate_basic(current)
    }
}
