#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

pub mod base;
pub mod phred;
pub mod region_string;
pub mod rms;
pub mod strand;

pub use {
    base::Base,
    phred::Phred,
    region_string::RegionString,
    rms::RootMeanSquare,
    strand::{Strand, StrandFromRecord},
};
