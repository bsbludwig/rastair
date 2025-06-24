#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

pub mod call;
pub mod vcf;
pub mod vcf_writer;
pub mod utils {
    mod base;
    pub use base::Base;
    pub mod file_helpers;
    mod region_string;
    pub use region_string::RegionString;
    mod rms;
    pub use rms::RootMeanSquare;
    mod base_modification;
    pub use base_modification::MethylatedPositions;
    mod phred;
    pub use phred::Phred;
    mod base_counter;
    pub use base_counter::Counter;
    mod strand;
    pub use strand::{Strand, StrandFromRecord};
}
pub mod sequence;
