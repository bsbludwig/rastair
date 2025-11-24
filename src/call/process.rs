//! Module for processing reads from BAM files to generate pileups and calculate metrics.

// First off, let's build Rastair-specific pileups from the BAM file
mod pileups;
pub use pileups::{PileupMappingParams, get_pileups};

// Next, we calculate various metrics based on the generated pileups
mod pileup_metrics;
pub use pileup_metrics::calculate_pileup_metrics;

// Then, we apply filters to the calculated pileup metrics
mod threshold_filters;
pub use threshold_filters::{ThresholdFilterParams, apply_threshold_filters};

// And machine learning-based metrics are just more filters
mod calc_ml;
pub use calc_ml::add_ml_metrics;

mod cpg_sites;
pub use cpg_sites::{propagate_cpg_pass_flags, set_denovo_adj};
