#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

// CLI modules
pub(crate) mod call;
pub use call::{CallParams, call};
mod call_reads;
pub use call_reads::{PerReadParams, call_reads};
mod convert;
pub use convert::{ConvertParams, convert};
mod mbias;
pub use bam::{BamRewriteArgs, rewrite as rewrite_bam};
pub use io::mpk::viewer::{MpkViewParams, view as mpk_view};
pub use mbias::{MBiasParams, mbias};

// utils
pub use utils::logging::setup_logging;
pub use vcf::Record as VcfRecord;

pub(crate) mod bam;
pub(crate) mod bed;
pub(crate) mod metrics;
pub(crate) mod train;
pub use train::{TrainModelParams, train_model};
pub(crate) mod io {
    pub mod formats;
    pub mod mpk;
    pub mod vcf_writer;
}
pub(crate) mod vcf;
pub(crate) mod utils {
    pub use rastair_types::*;

    mod base_counter;
    pub mod file_helpers;
    pub use base_counter::Counter;

    pub mod logging;

    mod surrounding;
    pub use surrounding::PileupMetricsIteratorExt;

    pub mod cli;

    mod conversion;
    pub use conversion::{IntoF64, default};

    mod grouping;
    pub use grouping::{ByAllele, ByStrand};

    mod sequence_context;
    pub use sequence_context::SequenceContext;
}
pub(crate) mod sequence;
