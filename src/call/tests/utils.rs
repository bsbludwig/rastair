//! Test utilities for writing concise and readable tests for rastair call
#![allow(clippy::cast_possible_truncation, reason = "Test code")]
use crate::call::process_region;
use crate::call::variant_calling::ErrorModel;
use crate::{
    CallParams,
    call::{
        pileup::{Pileup, SimpleRead, SimpleReads},
        process,
    },
    metrics::{PileupMetrics, ml::types::MachineLearning},
    sequence::{ChunkRegion, Region, Segment},
    utils::{SequenceContext, map_surrounding},
    vcf::{Contig, FieldConfig, emit_pileup, register},
};
pub use crate::{call::record_filters::RecordFilters, utils::default};
use clio::ClioPath;
use color_eyre::eyre::{ContextCompat as _, WrapErr as _};
pub(crate) use color_eyre::{Result, eyre::bail};
use rustc_hash::FxHashMap;
use seqair::vcf::{OutputFormat, Writer as SeqWriter};
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
                .map(|read_line| {
                    let before_base = pos_idx.checked_sub(1).map(|i| read_line.bases[i]);
                    let after_base = read_line.bases.get(pos_idx + 1).copied();
                    SimpleRead {
                        base: read_line.bases[pos_idx],
                        strand: read_line.strand,
                        before_base,
                        after_base,
                        ..default()
                    }
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
                noisy_ref_count: 0,
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
    /// Compare with another `FieldValue`, using epsilon for float comparisons
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
            // Both expected and actual are now VCF-style genotype strings ("0/1").
            // Normalize separator order is not needed; compare directly.
            match (self, other) {
                (Self::String(expected), Self::String(actual)) => expected == actual,
                _ => self.matches(other, epsilon),
            }
        } else {
            self.matches(other, epsilon)
        }
    }
}

/// Trait for converting values into `FieldValue` for test assertions
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

/// A parsed VCF data line, used to assert against emitted output.
#[derive(Debug, Clone)]
pub(crate) struct VcfRecord {
    pub r#ref: String,
    /// ALT alleles; `["."]` for a reference-only site.
    pub alt: Vec<String>,
    /// FILTER codes; empty or `["PASS"]` means the record passes.
    filters: Vec<String>,
    /// FORMAT key → first sample's value.
    format: FxHashMap<String, String>,
}

impl VcfRecord {
    fn passes(&self) -> bool {
        self.filters.is_empty() || self.filters == ["."] || self.filters == ["PASS"]
    }
}

/// Parse the data lines of a (plain text) VCF into [`VcfRecord`]s.
fn parse_vcf(text: &str) -> Vec<VcfRecord> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            let r#ref = cols.get(3).copied().unwrap_or(".").to_string();
            let alt = cols.get(4).copied().unwrap_or(".").split(',').map(String::from).collect();
            let filters = match cols.get(6).copied().unwrap_or(".") {
                "." => Vec::new(),
                f => f.split(';').map(String::from).collect(),
            };
            let mut format = FxHashMap::default();
            if let (Some(keys), Some(sample)) = (cols.get(8), cols.get(9)) {
                for (k, v) in keys.split(':').zip(sample.split(':')) {
                    format.insert(k.to_string(), v.to_string());
                }
            }
            VcfRecord { r#ref, alt, filters, format }
        })
        .collect()
}

/// Get a field value from a parsed VCF record by field ID (FORMAT fields only,
/// matching the previous behaviour).
fn get_field_value(record: &VcfRecord, field_id: &str) -> Result<FieldValue> {
    match field_id {
        "M5mC" | "DPM5mC" | "ADM5mC" => Ok(match record.format.get(field_id) {
            None => FieldValue::OptF64(None),
            Some(v) if v == "." => FieldValue::OptF64(None),
            Some(v) => {
                let values: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
                match values.as_slice() {
                    [single] => FieldValue::OptF64(Some(*single)),
                    _ => FieldValue::VecF64(values),
                }
            }
        }),
        "ML" => {
            let v = record
                .format
                .get("ML")
                .filter(|v| *v != ".")
                .ok_or_else(|| color_eyre::eyre::eyre!("ML field is empty"))?;
            let first = v
                .split(',')
                .next()
                .and_then(|x| x.parse::<f64>().ok())
                .ok_or_else(|| color_eyre::eyre::eyre!("ML field is empty"))?;
            Ok(FieldValue::F64(first))
        }
        "GT" => {
            let gt = record.format.get("GT").cloned().unwrap_or_else(|| ".".to_string());
            Ok(FieldValue::String(gt))
        }
        other => bail!("Unknown or unsupported field: {}", other),
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
            let actual_ref = actual.r#ref.as_str();
            let expected_ref: &str = expected.ref_base.into();
            if actual_ref != expected_ref {
                errors.push(format!(
                    "Record {}: Expected REF={}, got REF={}",
                    idx, expected_ref, actual_ref
                ));
            }

            // Check ALT
            let actual_alts = &actual.alt;
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
                    let actual_alt_str: &str = actual_alt.as_str();
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
                let actual_passes = actual.passes();

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
                            r.r#ref, r.alt, r.filters
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

    alt.filters.add(crate::vcf::RastairFilter::LowMlScore, || true);
}

pub(crate) fn reprocess(records: Vec<PileupMetrics>) -> Result<Vec<PileupMetrics>> {
    use crate::{
        call::variant_calling::ErrorModel, metrics, metrics::MethylationEvidenceStrandInfo,
    };

    let mut records = records;
    map_surrounding(
        &mut records,
        |before, current, after| {
            process::propagate_denovo_pass_flags(before, current, after, Some(ML_THRESHOLD))
        },
        "failed to propagate CpG pass flags",
    );

    records
        .into_iter()
        .map(|mut current| {
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

/// Emit the metrics to an in-memory plain-text VCF (all fields enabled) and
/// parse it back, so assertions run against the real encoded output.
pub(crate) fn metrics_to_vcf(
    metrics: &[PileupMetrics],
    filters: RecordFilters,
) -> Result<Vec<VcfRecord>> {
    let contigs = [Contig { name: "chr_test".into(), length: 100_000 }];
    let samples = [seqair_types::SmolStr::new("sample")];
    let (header, schema) = register(&contigs, &samples, &[])?;
    let config = FieldConfig::default().with_all_fields();
    let error_model = ErrorModel::default();

    let mut buf = Vec::new();
    {
        let mut writer = SeqWriter::new(&mut buf, OutputFormat::Vcf).write_header(&header)?;
        for metric in metrics {
            let contig = schema
                .contig(metric.contig_name())
                .wrap_err_with(|| format!("Contig {} not in header", metric.contig_name()))?;
            emit_pileup(
                metric,
                &schema,
                contig,
                &config,
                Some(ML_THRESHOLD),
                &error_model,
                &filters,
                &mut writer,
            )?;
        }
        writer.finish()?;
    }

    let text = String::from_utf8(buf).wrap_err("VCF output was not valid UTF-8")?;
    Ok(parse_vcf(&text))
}

/// Emit the given metrics as binary BCF bytes. Unlike [`metrics_to_vcf`], this
/// exercises the BCF encoding path, where the distinction between a present
/// FORMAT field with zero values per sample and a proper missing value matters
/// (the former makes htslib-based float readers panic on `chunks(0)`).
pub(crate) fn metrics_to_bcf(metrics: &[PileupMetrics], filters: RecordFilters) -> Result<Vec<u8>> {
    let contigs = [Contig { name: "chr_test".into(), length: 100_000 }];
    let samples = [seqair_types::SmolStr::new("sample")];
    let (header, schema) = register(&contigs, &samples, &[])?;
    let config = FieldConfig::default().with_all_fields();
    let error_model = ErrorModel::default();

    let mut buf = Vec::new();
    {
        let mut writer = SeqWriter::new(&mut buf, OutputFormat::Bcf).write_header(&header)?;
        for metric in metrics {
            let contig = schema
                .contig(metric.contig_name())
                .wrap_err_with(|| format!("Contig {} not in header", metric.contig_name()))?;
            emit_pileup(
                metric,
                &schema,
                contig,
                &config,
                Some(ML_THRESHOLD),
                &error_model,
                &filters,
                &mut writer,
            )?;
        }
        writer.finish()?;
    }

    Ok(buf)
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
