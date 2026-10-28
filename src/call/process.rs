//! Module for processing reads from BAM files to generate pileups and calculate metrics.

mod calc_ml;
mod denovo;
mod inference;
#[cfg(test)]
mod pileup_metrics;
mod pileups;
mod set_alt_calls;

#[cfg(test)]
pub use pileups::get_pileups;
pub use pileups::{ColumnKeep, Engine, PileupMappingParams, UncoveredCpgs};
pub use pileups::{HeldReads, get_pileups_held};

#[cfg(test)]
pub use pileup_metrics::calculate_pileup_metrics;

pub use calc_ml::{GPU_BATCH_BUFFER_SIZE, Unscored, ml_row, score_on_cpu};
// run on one dedicated thread when there is a GPU, so workers never own a forest
pub use inference::{GpuScoring, InferenceStage, score_on_gpu};

pub use denovo::mark_called_cpgs;
pub use set_alt_calls::set_alt_calls;
