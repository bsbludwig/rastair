//! Random Forest based classification of positions
//!
//! This module implements feature extraction for machine learning classification
//! of methylation sites in CpG contexts. The feature extraction logic replicates
//! Ben's R notebook analysis for training random forest models on VCF data.
//!
//! ## CpG Feature Extraction
//!
//! The `cpg` module contains functions to extract features from VCF records
//! including:
//! - Basic variant information (ref, alt, mapping quality, etc.)
//! - Sequence context one-hot encoding
//! - Normalized allele depths and strand bias counts
//! - Base and mapping quality metrics by strand
//! - Adjacent position features for methylation evidence
//!
//! The extracted features are returned as an ndarray Array1<f64> suitable
//! for input to a random forest classifier.

#![allow(unused)]

pub use cpg::params_from_record;
pub use models::{MachineLearning, MlResult};

mod models;
mod utils;

mod cpg;
mod denovo_cpg;
