//! Test utilities for writing concise and readable tests for rastair call
#![allow(clippy::cast_possible_truncation, reason = "Test code")]

use crate::call::variant_calling::ErrorModel;
use crate::utils::IntoF64 as _;
use crate::{
    CallParams,
    call::{
        pileup::{Pileup, SimpleRead, SimpleReads},
        process, process_region,
    },
    metrics::{PileupMetrics, ml::types::MachineLearning},
    sequence::{ChunkRegion, Region, Segment},
    utils::{PileupMetricsIteratorExt as _, SequenceContext},
    vcf::Record as VcfRecord,
};
pub use crate::{call::record_filters::RecordFilters, utils::default};
use clio::ClioPath;
use color_eyre::eyre::ContextCompat as _;
pub(crate) use color_eyre::{Result, eyre::bail};
use seqair_types::{Base, Probability, Strand};
use std::{rc::Rc, str::FromStr, sync::OnceLock};

pub const ML_THRESHOLD: Probability = Probability::new_panicky(0.5);

#[macro_export]
macro_rules! pileups {
    (
        [ $($ref_base:ident)+ ] Ref,
        $([ $($base:ident)+ ] $strand:ident ),+ $(,)?
    ) => {{
        $crate::call::tests::utils::create_pileups(
            vec![$( seqair_types::Base::$ref_base ),+],
            vec![$(
                $crate::call::tests::utils::ReadLine {
                    bases: vec![$( seqair_types::Base::$base ),+],
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
            overlap_start: 0,
            overlap_end: 0,
        },
        sequence: ref_bases.iter().map(|b| *b as u8).collect(),
        overlap_start: 0,
        overlap_end: 0,
    };

    let pileups = ref_bases
        .into_iter()
        .enumerate()
        .take(num_positions)
        .map(|(pos_idx, reference_base)| {
            let pos = start + pos_idx as u64;

            let reads: Vec<SimpleRead> = read_lines
                .iter()
                .map(|read_line| SimpleRead {
                    base: read_line.bases[pos_idx],
                    strand: read_line.strand,
                    ..default()
                })
                .collect();

            let context = SequenceContext::new(pos_idx, &segment).expect("valid context");

            Pileup {
                region: segment.range.clone(),
                context,
                pos: pos as u32,
                reads: SimpleReads(reads.into()),
                reference_base,
                indel_observations: Default::default(),
                depth_offset: 0,
                homopolymer_run: 0,
                dinucleotide_run: 0,
                soft_clip_count: 0,
                indel_ref_window: Default::default(),
                indel_ref_anchor: 0,
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
macro_rules! vcf_assert {
    () => {{
        $crate::call::tests::utils::VcfMatcher { expected: vec![] }
    }};
    // Pattern with explicit PASS/FAIL status and field assertions
    (
        $(( $ref:ident $($alt:tt),+ $(,)? ) $pass_status:ident $($field:ident = $value:expr)*),+ $(,)?
    ) => {{
        $crate::call::tests::utils::VcfMatcher {
            expected: vec![$(
                $crate::call::tests::utils::ExpectedVcfRecord {
                    ref_base: seqair_types::Base::$ref,
                    alt_bases: vec![$( $crate::call::tests::utils::parse_alt_token(stringify!($alt)) ),+],
                    pass_status: Some($crate::call::tests::utils::parse_pass_status(stringify!($pass_status))),
                    fields: vec![$( (stringify!($field), $value.to_field_value()) ),*],
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
                    ref_base: seqair_types::Base::$ref,
                    alt_bases: vec![$( $crate::call::tests::utils::parse_alt_token(stringify!($alt)) ),+],
                    pass_status: None,
                    fields: vec![],
                }
            ),+]
        }
    }};
}

pub(crate) fn parse_alt_token(s: &str) -> Option<Base> {
    if s == "." { None } else { seqair_types::Base::from_str(s).ok() }
}

pub(crate) fn parse_pass_status(s: &str) -> bool {
    match s {
        "PASS" => true,
        "FAIL" => false,
        _ => panic!("Expected PASS or FAIL, got: {}", s),
    }
}

/// Value types for field assertions
#[derive(Debug, Clone)]
pub(crate) enum FieldValue {
    F64(f64),
    OptF64(Option<f64>),
    VecF64(Vec<f64>),
    String(String),
}

impl FieldValue {
    /// Compare with another FieldValue, using epsilon for float comparisons
    fn matches(&self, other: &Self, epsilon: f64) -> bool {
        match (self, other) {
            (Self::F64(a), Self::F64(b)) => (a - b).abs() < epsilon,
            // Allow F64 to match OptF64(Some(value)) for convenience
            (Self::F64(a), Self::OptF64(Some(b))) | (Self::OptF64(Some(b)), Self::F64(a)) => {
                (a - b).abs() < epsilon
            }
            (Self::OptF64(a), Self::OptF64(b)) => match (a, b) {
                (Some(a), Some(b)) => (a - b).abs() < epsilon,
                (None, None) => true,
                _ => false,
            },
            (Self::VecF64(a), Self::VecF64(b)) => {
                if a.len() != b.len() {
                    return false;
                }
                a.iter().zip(b.iter()).all(|(a_val, b_val)| (a_val - b_val).abs() < epsilon)
            }
            (Self::String(a), Self::String(b)) => a == b,
            _ => false,
        }
    }

    /// Compare field values with special handling for GT field
    fn matches_field(&self, other: &Self, field_name: &str, epsilon: f64) -> bool {
        if field_name == "GT" {
            // For GT field, convert VCF-style strings like "0/1" to expected format
            match (self, other) {
                (Self::String(expected), Self::String(actual)) => {
                    let normalized_expected = Self::normalize_gt_string(expected);
                    normalized_expected == *actual
                }
                _ => self.matches(other, epsilon),
            }
        } else {
            self.matches(other, epsilon)
        }
    }

    /// Normalize a genotype string - converts "0/1" to "Genotype([Unphased(0), Unphased(1)])"
    fn normalize_gt_string(s: &str) -> String {
        // If it's already in the expected format, return as-is
        if s.starts_with("Genotype(") {
            return s.to_string();
        }

        // Parse VCF-style format like "0/1" or "1/1"
        let parts: Vec<&str> = s.split('/').collect();
        if parts.len() != 2 {
            // Not a genotype format, return as-is
            return s.to_string();
        }

        if let (Ok(allele1), Ok(allele2)) = (parts[0].parse::<i32>(), parts[1].parse::<i32>()) {
            format!("Genotype([Unphased({}), Unphased({})])", allele1, allele2)
        } else {
            // Parse failed, return as-is
            s.to_string()
        }
    }
}

/// Trait for converting values into FieldValue for test assertions
pub(crate) trait ToFieldValue {
    fn to_field_value(self) -> FieldValue;
}

impl ToFieldValue for f64 {
    fn to_field_value(self) -> FieldValue {
        FieldValue::F64(self)
    }
}

impl ToFieldValue for Option<f64> {
    fn to_field_value(self) -> FieldValue {
        FieldValue::OptF64(self)
    }
}

impl ToFieldValue for &str {
    fn to_field_value(self) -> FieldValue {
        FieldValue::String(self.to_string())
    }
}

impl ToFieldValue for String {
    fn to_field_value(self) -> FieldValue {
        FieldValue::String(self)
    }
}

impl ToFieldValue for Vec<f64> {
    fn to_field_value(self) -> FieldValue {
        FieldValue::VecF64(self)
    }
}

/// Get a field value from a VCF record by field ID
fn get_field_value(record: &VcfRecord, field_id: &str) -> Result<FieldValue> {
    // FORMAT fields (sample 0)
    if let Some(sample) = record.samples.first() {
        match field_id {
            "M5mC" => {
                use crate::vcf::Methylated;
                return Ok(match &sample.methylated {
                    // Unknown: no data available
                    Methylated::Unknown => FieldValue::OptF64(None),
                    // NoEvidence: checked and found no methylation (beta = 0.0)
                    Methylated::NoEvidence => FieldValue::OptF64(Some(0.0)),
                    // Single context: one beta value
                    Methylated::OriginalCpG { beta } | Methylated::DeNovoCpG { beta } => {
                        FieldValue::OptF64(Some(beta.f()))
                    }
                    // Dual context: both beta values
                    Methylated::Both { original_beta, denovo_beta } => {
                        FieldValue::VecF64(vec![original_beta.f(), denovo_beta.f()])
                    }
                });
            }
            "ML" => {
                // ML is OnePerAlt, so we get the first value
                let ml_value = sample
                    .machine_learning_prediction
                    .first()
                    .ok_or_else(|| color_eyre::eyre::eyre!("ML field is empty"))?;
                return Ok(FieldValue::F64(*ml_value));
            }
            "GT" => {
                return Ok(FieldValue::String(format!("{:?}", sample.genotype)));
            }
            _ => {}
        }
    }

    bail!("Unknown or unsupported field: {}", field_id)
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
    pub fields: Vec<(&'static str, FieldValue)>,
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
            let expected_alts_count = expected.alt_bases.len();
            if expected_alts_count != actual_alts.len() {
                errors.push(format!(
                    "Record {}: Expected {} ALT alleles, got {}",
                    idx,
                    expected_alts_count,
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
                            "Record {}: Expected PASS, but record did not pass {:?}",
                            idx, actual_filters
                        ));
                    }
                } else {
                    // Expect FAIL (has filters)
                    if actual_passes {
                        errors.push(format!("Record {}: Expected FAIL, but got PASS", idx));
                    }
                }
            }

            // Check field assertions
            for (field_id, expected_value) in &expected.fields {
                match get_field_value(actual, field_id) {
                    Ok(actual_value) => {
                        if !expected_value.matches_field(&actual_value, field_id, 1e-3) {
                            errors.push(format!(
                                "Record {}: Field {} expected {:?}, got {:?}",
                                idx, field_id, expected_value, actual_value
                            ));
                        }
                    }
                    Err(e) => {
                        errors.push(format!(
                            "Record {}: Failed to get field {}: {}",
                            idx, field_id, e
                        ));
                    }
                }
            }
        }

        if !errors.is_empty() {
            let records = if actual.is_empty() {
                "[No records]".to_string()
            } else {
                format!(
                    "Got these VCF lines:\n    {}",
                    actual
                        .iter()
                        .enumerate()
                        .map(|(i, r)| format!(
                            "({i}) ref={} alt={:?} {:?}",
                            r.main.r#ref, r.main.alt, r.filters
                        ))
                        .collect::<Vec<_>>()
                        .join("\n    ")
                )
            };
            bail!(
                "VCF matching failed.\n    {records}\n    Errors:\n    - {}",
                errors.join("\n    - ")
            );
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
            regions: Some("chr_test".parse().unwrap()),
        },
        record_filters: params,
        segmentation: default(),
        require_tags: default(),
        variant_calling: default(),
        indel: default(),
        denovo_cpg: default(),
        methylation: default(),
        ml: default(),
        vcf: default(),
        bed: default(),
        total_threads: 2,
    };

    static ML: OnceLock<MachineLearning> = OnceLock::new();
    let ml = ML.get_or_init(|| params.ml.init().unwrap());

    process_region(segment, pileups.into_iter(), &params, ml)
}

#[track_caller]
pub(crate) fn set_pass(m: &mut PileupMetrics, base: Base) {
    assert!(
        m.ref_base() != base,
        "Cannot set reference base {base} as passing - only alt bases can be marked as passing/failing"
    );
    let alt = m.alt_filters_mut(base).wrap_err_with(|| format!("no {base} alt")).unwrap();
    alt.ml = Some(Probability::ONE);
}

rastair_vcf::filter!(MANUAL, "manual");

#[track_caller]
pub(crate) fn set_fail(m: &mut PileupMetrics, base: Base) {
    assert!(
        m.ref_base() != base,
        "Cannot set reference base {base} as failing - only alt bases can be marked as passing/failing"
    );
    m.pos_filters.other_pos_in_denovo_passes = false;

    let alt = m.alt_filters_mut(base).wrap_err_with(|| format!("no {base} alt")).unwrap();
    alt.filters.other_pos_in_denovo_passes = false;
    alt.ml = Some(Probability::ZERO);

    alt.filters.add(MANUAL, || true);
}

pub(crate) fn reprocess(records: Vec<PileupMetrics>) -> Result<Vec<PileupMetrics>> {
    use crate::{
        call::variant_calling::ErrorModel, metrics, metrics::MethylationEvidenceStrandInfo,
    };

    records
        .into_iter()
        .map_surrounding(|before, current, after| {
            process::propagate_denovo_pass_flags(before, current, after, Some(ML_THRESHOLD))
        })
        .map(|current| {
            let mut current = current?;
            process::set_alt_calls(&mut current, Some(ML_THRESHOLD))?;
            process::add_position_tags(&mut current);

            // Recalculate methylation strand info (needed for genotype estimation)
            current.pos_metrics.extended.methylation_strand_info =
                MethylationEvidenceStrandInfo::from_pileup(&current);

            // Recalculate genotype and methylation after changing alt calls
            current.pos_metrics.extended.genotype =
                current.estimate_genotype(Some(ML_THRESHOLD), ErrorModel::default());
            current.pos_metrics.extended.methylated =
                metrics::methylation::call(&current)?.unwrap_or_default();
            current.pos_metrics.extended.methylation_strand_info =
                MethylationEvidenceStrandInfo::from_pileup_with_methylation(&current);

            Ok(current)
        })
        .collect()
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

    #[expect(unused, reason = "for completeness")]
    pub(crate) fn all_cpgs() -> Self {
        Self { vcf_all: true, cpgs_only: true }
    }
}

pub(crate) fn metrics_to_vcf(
    metrics: &[PileupMetrics],
    filters: RecordFilters,
) -> Result<Vec<VcfRecord>> {
    let mut vcf_records = Vec::new();
    for metric in metrics {
        let records = metric.to_vcf_records(Some(ML_THRESHOLD), &ErrorModel::default())?;
        vcf_records.extend(records.to_vec(&filters).into_iter().cloned());
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
