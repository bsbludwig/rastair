//! Segment processing for genomic analysis.
//!
//! This module provides functionality for:
//! - Dividing genomic regions into manageable segments
//! - Reading reference sequences from FASTA files
//! - Accessing aligned reads from BAM files
//! - Processing segments with configurable overlap between chunks

mod chunked;
mod readers;
mod regions;
mod segementation;
#[cfg(test)]
mod tests;

pub use readers::{ReaderParams, Readers};
pub use regions::{ChunkRegion, Region, SelectedRegion};
pub use segementation::{Segment, SegmentationParams};
