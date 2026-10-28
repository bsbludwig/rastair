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
mod segmentation;
#[cfg(test)]
mod tests;

#[cfg(test)]
pub(crate) use readers::RastairRecordFilter;
pub use readers::{
    IndexDepth, RastairReadExtras, ReaderParams, Readers, ReferenceWindow, SeqairReaders,
    StrandSource,
};
pub use regions::{ChunkRegion, Region, SelectedRegion};
pub use segmentation::{ByteBudget, Segment, SegmentationParams, segment_length};
