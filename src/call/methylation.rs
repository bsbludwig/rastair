//! Methylation calling
//!
//! Right now, Rastair collects all possible variant positions from aligned reads relative to a given reference genome.
//! These are represented as the same `Record`s that can also be converted to `VCF` lines.
//! The methylation calling code plugs in at the point where the `Record`s are being processed,
//! i.e. it is called with an already constructed `Record` object.
//! If a methylation call is made, the `Record` is updated.

mod filters;
pub use filters::add_filters;

pub mod params;
pub use params::ThresholdParams;
