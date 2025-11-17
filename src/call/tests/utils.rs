//! Test utilities for writing concise and readable tests for rastair call
#![allow(clippy::cast_possible_truncation, reason = "Test code")]

use crate::{
    CallParams,
    call::{
        pileup::{Pileup, SimpleRead, SimpleReads},
        process_region,
    },
    metrics::PileupMetrics,
    sequence::{ChunkRegion, Region, Segment},
    utils::SequenceContext,
    vcf::Record as VcfRecord,
};
pub use crate::{call::record_filters::RecordFilters, utils::default};
use clio::ClioPath;
pub(crate) use color_eyre::{Result, eyre::bail};
use rastair_types::{Base, Probability, SmallVec, Strand};
use std::{rc::Rc, str::FromStr};

#[macro_export]
macro_rules! pileups {
    (
        [ $($ref_base:ident)+ ] Ref,
        $([ $($base:ident)+ ] $strand:ident ),+ $(,)?
    ) => {{
        $crate::call::tests::utils::create_pileups(
            vec![$( rastair_types::Base::$ref_base ),+],
            vec![$(
                $crate::call::tests::utils::ReadLine {
                    bases: vec![$( rastair_types::Base::$base ),+],
                    strand: $crate::call::tests::utils::parse_strand(stringify!($strand))
                }
            ),+]
        )
    }};
}

pub(crate) struct ReadLine {
    pub bases: Vec<Base>,
    pub strand: Strand,
}

pub(crate) fn create_pileups(
    ref_bases: Vec<Base>,
    read_lines: Vec<ReadLine>,
) -> (Segment, Vec<Pileup>) {
    let num_positions = ref_bases.len();

    // Validate all reads have the same number of bases
    for read in &read_lines {
        assert_eq!(
            read.bases.len(),
            num_positions,
            "All reads must have the same number of bases as the reference"
        );
    }

    let start = 1000u64;
    let end = start + num_positions as u64;

    let segment = Segment {
        range: ChunkRegion {
            region: Region { contig: "chr_test".into(), start, end },
            last_position: end,
        },
        sequence: ref_bases.iter().map(|b| *b as u8).collect(),
    };

    let pileups = ref_bases
        .into_iter()
        .enumerate()
        .take(num_positions)
        .map(|(pos_idx, reference_base)| {
            let pos = start + pos_idx as u64;

            let reads = read_lines
                .iter()
                .enumerate()
                .map(|(read_idx, read_line)| SimpleRead {
                    base: read_line.bases[pos_idx],
                    strand: read_line.strand,
                    qname: SmallVec::from(format!("read{}", read_idx).into_bytes()),
                    ..default()
                })
                .collect();

            let context = SequenceContext::new(pos_idx, &segment).expect("valid context");

            Pileup {
                region: segment.range.clone(),
                context,
                pos: pos as u32,
                reads: SimpleReads(reads),
                reference_base,
            }
        })
        .collect();

    (segment, pileups)
}

pub(crate) fn parse_strand(s: &str) -> Strand {
    match s {
        "OT" => Strand::OT,
        "OB" => Strand::OB,
        _ => panic!("Unknown strand: {}", s),
    }
}

#[macro_export]
macro_rules! vcf {
    () => {{
        $crate::call::tests::utils::VcfMatcher { expected: vec![] }
    }};
    // Pattern with explicit PASS/FAIL status
    (
        $(( $ref:ident $($alt:tt),+ $(,)? ) $pass_status:ident),+ $(,)?
    ) => {{
        $crate::call::tests::utils::VcfMatcher {
            expected: vec![$(
                $crate::call::tests::utils::ExpectedVcfRecord {
                    ref_base: rastair_types::Base::$ref,
                    alt_bases: vec![$( $crate::call::tests::utils::parse_alt_token(stringify!($alt)) ),+],
                    pass_status: Some($crate::call::tests::utils::parse_pass_status(stringify!($pass_status))),
                }
            ),+]
        }
    }};
    // Pattern without pass/fail status (backward compatibility)
    (
        $(( $ref:ident $($alt:tt),+ $(,)? )),+ $(,)?
    ) => {{
        $crate::call::tests::utils::VcfMatcher {
            expected: vec![$(
                $crate::call::tests::utils::ExpectedVcfRecord {
                    ref_base: rastair_types::Base::$ref,
                    alt_bases: vec![$( $crate::call::tests::utils::parse_alt_token(stringify!($alt)) ),+],
                    pass_status: None,
                }
            ),+]
        }
    }};
}

pub(crate) fn parse_alt_token(s: &str) -> Option<Base> {
    if s == "." { None } else { rastair_types::Base::from_str(s).ok() }
}

pub(crate) fn parse_pass_status(s: &str) -> bool {
    match s {
        "PASS" => true,
        "FAIL" => false,
        _ => panic!("Expected PASS or FAIL, got: {}", s),
    }
}

#[derive(Debug)]
pub(crate) struct VcfMatcher {
    pub expected: Vec<ExpectedVcfRecord>,
}

#[derive(Debug)]
pub(crate) struct ExpectedVcfRecord {
    pub ref_base: Base,
    pub alt_bases: Vec<Option<Base>>,
    pub pass_status: Option<bool>,
}

impl VcfMatcher {
    /// Check if actual VCF records match expected records
    #[track_caller]
    pub(crate) fn matches(&self, actual: Vec<VcfRecord>) -> Result<()> {
        let mut errors = Vec::new();

        // Check count
        if self.expected.len() != actual.len() {
            errors.push(format!(
                "Expected {} VCF records, got {}",
                self.expected.len(),
                actual.len()
            ));
        }

        // Check each record
        for (idx, (expected, actual)) in self.expected.iter().zip(actual.iter()).enumerate() {
            // Check REF
            let actual_ref = actual.main.r#ref.as_str();
            let expected_ref: &str = expected.ref_base.into();
            if actual_ref != expected_ref {
                errors.push(format!(
                    "Record {}: Expected REF={}, got REF={}",
                    idx, expected_ref, actual_ref
                ));
            }

            // Check ALT
            let actual_alts = &actual.main.alt;
            let expected_alts = expected.alt_bases.iter().filter(|x| x.is_some()).count();
            if expected_alts != actual_alts.len() {
                errors.push(format!(
                    "Record {}: Expected {} ALT alleles, got {}",
                    idx,
                    expected.alt_bases.len(),
                    actual_alts.len()
                ));
            } else {
                for (alt_idx, (expected_alt, actual_alt)) in
                    expected.alt_bases.iter().zip(actual_alts.iter()).enumerate()
                {
                    let actual_alt_str = actual_alt.as_str();
                    match expected_alt {
                        None => {
                            if actual_alt_str != "." {
                                errors.push(format!(
                                    "Record {}, ALT {}: Expected '.', got '{}'",
                                    idx, alt_idx, actual_alt_str
                                ));
                            }
                        }
                        Some(base) => {
                            let expected_str: &str = (*base).into();
                            if actual_alt_str != expected_str {
                                errors.push(format!(
                                    "Record {}, ALT {}: Expected '{}', got '{}'",
                                    idx, alt_idx, expected_str, actual_alt_str
                                ));
                            }
                        }
                    }
                }
            }

            // Check FILTER status if expected
            if let Some(expected_pass) = expected.pass_status {
                let actual_filters = &actual.filters;
                let actual_passes = actual_filters.pass();

                if expected_pass {
                    // Expect PASS
                    if !actual_passes {
                        errors.push(format!(
                            "Record {}: Expected PASS, but record did not pass filters (len={})",
                            idx,
                            actual_filters.len()
                        ));
                    }
                } else {
                    // Expect FAIL (has filters)
                    if actual_passes {
                        errors.push(format!("Record {}: Expected FAIL, but got PASS", idx));
                    }
                }
            }

            // Todo: Check INFO fields
        }

        if !errors.is_empty() {
            if !actual.is_empty() {
                errors.push(format!(
                    "Got these VCF lines: {}",
                    actual
                        .iter()
                        .map(|r| format!("ref={} alt={:?}", r.main.r#ref, r.main.alt))
                        .collect::<Vec<_>>()
                        .join("\n")
                ));
            }
            bail!("VCF matching failed:\n{}", errors.join("\n"));
        }

        Ok(())
    }
}

/// Helper function to run `process_region` with fake pileup data
pub(crate) fn test_call(
    segment: Segment,
    pileups: Vec<Pileup>,
    params: RecordFilters,
) -> Result<Vec<PileupMetrics>> {
    let segment = Rc::new(segment);
    let params = CallParams {
        segments: crate::sequence::ReaderParams {
            bam_file: ClioPath::new("tests/data/test.bam").unwrap(),
            fasta_file: ClioPath::new("tests/data/test.fasta.gz").unwrap(),
            region: Some("chr_test".parse()?),
        },
        record_filters: params,
        segmentation: default(),
        variant_calling: default(),
        denovo_cpg: default(),
        methylation: default(),
        ml: default(),
        vcf: default(),
        bed: default(),
        total_threads: 2,
    };

    let ml = params.ml.init()?;

    process_region(segment, pileups.into_iter(), &params, &ml)
}

impl RecordFilters {
    pub(crate) fn variants() -> Self {
        Self { vcf_all: false, cpgs_only: false }
    }

    pub(crate) fn all() -> Self {
        Self { vcf_all: true, cpgs_only: false }
    }

    pub(crate) fn cpgs() -> Self {
        Self { vcf_all: false, cpgs_only: true }
    }

    pub(crate) fn all_cpgs() -> Self {
        Self { vcf_all: true, cpgs_only: true }
    }
}

pub(crate) fn metrics_to_vcf(metrics: &[PileupMetrics]) -> Result<Vec<VcfRecord>> {
    let mut vcf_records = Vec::new();
    for metric in metrics {
        let mut records = metric.to_vcf_records(None)?;
        vcf_records.append(&mut records);
    }
    Ok(vcf_records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use Base::*;

    #[test]
    fn test_parse_strand() {
        assert_eq!(parse_strand("OT"), Strand::OT);
        assert_eq!(parse_strand("OB"), Strand::OB);
    }

    #[test]
    fn test_create_simple_pileups() -> Result<()> {
        let (segment, pileups) = create_pileups(
            vec![A, C],
            vec![
                ReadLine { bases: vec![A, T], strand: Strand::OT },
                ReadLine { bases: vec![A, C], strand: Strand::OB },
            ],
        );

        assert_eq!(segment.sequence.len(), 2);
        assert_eq!(pileups.len(), 2);
        assert_eq!(pileups[0].reference_base, A);
        assert_eq!(pileups[1].reference_base, C);
        assert_eq!(pileups[0].reads.len(), 2);
        assert_eq!(pileups[1].reads.len(), 2);

        Ok(())
    }
}
