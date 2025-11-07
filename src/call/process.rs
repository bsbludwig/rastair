//! Module for processing reads from BAM files to generate pileups and calculate metrics.

// First off, let's build Rastair-specific pileups from the BAM file
mod pileups;
pub use pileups::{PileupMappingParams, get_pileups};

// Next, we calculate various metrics based on the generated pileups
mod pileup_metrics;
pub use pileup_metrics::{PileupMetricsParams, calculate_pileup_metrics};

// Finally, we add machine learning-based metrics to the pileups
mod calc_ml;
pub use calc_ml::add_ml_metrics;
