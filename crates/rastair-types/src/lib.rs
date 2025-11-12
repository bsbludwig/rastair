//! This crate provides common types used throughout rastair2.
#![deny(missing_docs)]
#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

mod base;
mod phred;
mod probability;
mod region_string;
mod rms;
mod strand;

pub use {
    base::Base,
    phred::Phred,
    probability::Probability,
    region_string::RegionString,
    rms::{RootMeanSquare, RootMeanSquareExt},
    strand::{Strand, StrandFromRecord},
};
