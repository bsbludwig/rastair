#![expect(clippy::cast_possible_truncation, reason = "Test code")]

use super::*;
use crate::utils::Base;
use clio::ClioPath;
use color_eyre::Result;
use proptest::prelude::*;
use rust_htslib::bam::FetchDefinition;
use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("data")
}

fn get_test_bam() -> ClioPath {
    ClioPath::new(test_data_dir().join("test.bam")).expect("test bam path should be valid")
}

fn get_test_fasta() -> ClioPath {
    ClioPath::new(test_data_dir().join("test.fasta.gz")).expect("test fasta path should be valid")
}

#[test]
fn test_segment_reading() -> Result<()> {
    // Create a test region and params
    let region = ChunkRegion {
        region: Region { contig: "chr19".into(), start: 6105700, end: 6105800 },
        last_position: 6105900,
    };

    let params =
        ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), region: None };

    // Initialize readers
    let mut readers = params.readers()?;

    // Test reading a single segment
    let segment = readers.segment(&region, 2)?;

    // Verify segment properties
    assert_eq!(segment.range, region);

    Ok(())
}

#[test]
fn test_calculate_segments() -> Result<()> {
    let params = ReaderParams {
        bam_file: get_test_bam(),
        fasta_file: get_test_fasta(),
        region: Some("chr19:6105700-6105800".parse().unwrap()),
    };

    let readers = params.readers()?;
    // Collect all segments into a Vec since we need to verify properties across all segments
    let segments = readers.segments(1000, 100)?.collect::<Vec<_>>();

    assert!(!segments.is_empty(), "Should have at least one segment");

    // Check segment properties
    for segment in segments {
        assert_eq!(segment.region.contig, "chr19");
        assert!(segment.region.start >= 6105700);
        assert!(segment.region.end <= 6105800);
        // last position is the end of the contig, the same for all segments
        assert_eq!(segment.last_position, 61431566);
    }

    Ok(())
}

#[test]
fn test_overlapping_segments() -> Result<()> {
    let params = ReaderParams {
        bam_file: get_test_bam(),
        fasta_file: get_test_fasta(),
        region: Some("chr19:6105700-6105900".parse().unwrap()),
    };

    let mut readers = params.readers()?;
    let segment_max_length = 100; // Small max length to force multiple segments
    let segment_overlap = 20;
    let segments = readers.segments(segment_max_length, segment_overlap)?.collect::<Vec<_>>();

    assert!(segments.len() > 1, "Should have multiple segments");

    // Check overlaps between adjacent segments
    for pair in segments.windows(2) {
        let first = readers.segment(&pair[0], 2)?;
        let second = readers.segment(&pair[1], 2)?;

        // Verify overlap amount
        let overlap = first.range.region.end - second.range.region.start;
        assert_eq!(overlap, segment_overlap, "Overlap amount should match configured value");

        // Calculate overlap regions accounting for 0-based sequence indexing
        let first_start_idx =
            (first.range.region.end - overlap - first.range.region.start) as usize;
        let first_end_idx = (first.range.region.end - first.range.region.start) as usize;
        let second_start_idx = 0;
        let second_end_idx = overlap as usize;

        // Verify overlapping sequence content matches
        let overlap_first = &first.sequence[first_start_idx..first_end_idx];
        let overlap_second = &second.sequence[second_start_idx..second_end_idx];
        assert_eq!(overlap_first, overlap_second, "Overlapping sequences should match");
    }

    Ok(())
}

#[test]
fn test_segment_to_fetch_definition() -> Result<()> {
    // Create a test segment
    let region = ChunkRegion {
        region: Region { contig: "chr19".into(), start: 6105700, end: 6105800 },
        last_position: 6105900,
    };
    let segment = Segment { range: region.clone(), sequence: vec![65, 66, 67] }; // "ABC"

    // Convert to FetchDefinition
    let fetch_def = FetchDefinition::try_from(&segment.region)?;

    // Verify conversion worked correctly
    if let FetchDefinition::RegionString(chr, start, end) = fetch_def {
        assert_eq!(chr, b"chr19");
        assert_eq!(start, 6105700);
        assert_eq!(end, 6105801); // this is exclusive end
    } else {
        panic!("Unexpected FetchDefinition variant");
    }

    Ok(())
}

#[test]
fn test_segment_error_handling() -> Result<()> {
    let params =
        ReaderParams { bam_file: get_test_bam(), fasta_file: get_test_fasta(), region: None };

    let mut readers = params.readers()?;

    // Test with an invalid region (non-existent chromosome)
    let invalid_region = ChunkRegion {
        region: Region { contig: "nonexistent".into(), start: 100, end: 200 },
        last_position: 300,
    };

    let result = readers.segment(&invalid_region, 2);
    assert!(result.is_err());

    Ok(())
}

#[test]
fn test_non_overlapping_segments() -> Result<()> {
    let params = ReaderParams {
        bam_file: get_test_bam(),
        fasta_file: get_test_fasta(),
        region: Some("chr19:6105700-6105900".parse().unwrap()),
    };

    let readers = params.readers()?;

    let segment_max_length = 50; // Small max length to force multiple segments
    let segment_overlap = 0; // No overlap
    let segments = readers.segments(segment_max_length, segment_overlap)?.collect::<Vec<_>>();

    assert!(segments.len() > 1, "Should have multiple segments");

    // Check no overlaps between adjacent segments
    for pair in segments.windows(2) {
        // The first segment should end exactly where the next one starts
        assert_eq!(pair[0].region.end, pair[1].region.start);
    }

    Ok(())
}

#[test]
fn test_sequence_slice() -> Result<()> {
    let segment = Segment {
        range: ChunkRegion {
            region: Region { contig: "chr19".into(), start: 6105700, end: 6105800 },
            last_position: 6105900,
        },
        sequence: b"ATCGG".into(),
    };

    // Test valid slice
    let slice = segment.sequence_slice::<5>(1, 4)?;
    assert_eq!(slice.as_slice(), &[Base::T, Base::C, Base::G]);

    // Test out-of-bounds slice
    let slice = segment.sequence_slice::<5>(0, 10)?;
    assert_eq!(slice.as_slice(), &[Base::A, Base::T, Base::C, Base::G, Base::G]);

    Ok(())
}

proptest!(
    #[test]
    fn proptest_sequence_slice(start in 0usize..100, end in 0usize..100, seq in "[ATCG]{0,10}") {
        prop_assume!(start <= end, "Start must be less than or equal to end");
        let segment = Segment {
            range: ChunkRegion {
                region: Region { contig: "chr19".into(), start: 6105700, end: 6105800 },
                last_position: 6105900,
            },
            sequence: seq.into_bytes(),
        };

        // Ensure we don't panic on out-of-bounds slices
        segment.sequence_slice::<5>(start, end).unwrap();
    }
);
