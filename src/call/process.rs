//! Module for processing reads from BAM files to generate pileups and calculate metrics.

mod calc_ml;
mod denovo;
mod pileup_metrics;
mod pileups;
mod position_tags;
mod set_alt_calls;
mod threshold_filters;

// First off, let's build Rastair-specific pileups from the BAM file
pub use pileups::{PileupMappingParams, get_pileups};

// Next, we calculate various metrics based on the generated pileups
pub use pileup_metrics::calculate_pileup_metrics;
// and set de-novo adjacency information now that we know where de-novo candidates are
pub use denovo::set_denovo_adj;

// Then, we apply filters to the calculated pileup metrics
pub use threshold_filters::{ThresholdFilterParams, apply_threshold_filters};

// And machine learning-based metrics, which are just more filters
pub use calc_ml::add_ml_metrics;

// Now, let's mark both positions of de-novo CpG sites as pass if one passes
pub use denovo::propagate_denovo_pass_flags;

// Finally, we set the actual variant calls based on all the metrics and filters
pub use set_alt_calls::set_alt_calls;
// and add tags for positions based on the alt calls
pub use position_tags::add_position_tags;
