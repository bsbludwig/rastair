#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

// CLI modules
mod call;
pub use call::{CallParams, SegmentationParams, call};
mod call_reads;
pub use call_reads::{PerReadParams, call_reads};
mod convert;
pub use convert::{ConvertParams, convert};
mod mbias;
pub use mbias::{MBiasParams, mbias};

pub mod bam;
pub mod bed;
pub mod metrics;
pub mod io {
    pub mod formats;
    pub mod mpk;
    pub mod vcf_writer;
}
pub mod vcf;
pub mod utils {
    pub use rastair_types::*;

    mod base_counter;
    pub mod file_helpers;
    pub use base_counter::Counter;

    pub mod logging;

    mod surrounding;
    pub use surrounding::{Surrounding, surrounding_pileups, surrounding_records};

    mod dedupe_reads;
    pub use dedupe_reads::ReadDeduplicator;

    pub mod cli;

    mod conversion;
    pub use conversion::{IntoF64, default};

    mod grouping;
    pub use grouping::{ByAllele, ByStrand};

    mod sequence_context;
    pub use sequence_context::SequenceContext;
}
pub mod sequence;
