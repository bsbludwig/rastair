#![allow(clippy::print_stdout, reason = "verify prints its results to stdout")]

use crate::{
    utils::cli,
    vcf::{DeNovoCpGCandidate, InCpG, Methylated},
};
use clio::ClioPath;
use color_eyre::eyre::{Result, WrapErr, ensure, eyre};
use rastair_types::{RegionString, SmolStr};
use rastair_vcf::VcfField as _;
use rust_htslib::bcf::{self, Read as _, header::HeaderView};
use rustc_hash::{FxHashMap, FxHashSet};
use std::{
    num::NonZeroU64,
    path::{Path, PathBuf},
    thread::available_parallelism,
};
use tracing::{info, instrument, warn};

// ─── CLI params ────────────────────────────────────────────────────────────

#[derive(Debug, clap::Args)]
pub struct VerifyParams {
    /// Predictions VCF file (output from rastair call)
    #[arg(help_heading = cli::sections::INPUT, value_hint = clap::ValueHint::FilePath)]
    predictions: ClioPath,

    /// Ground truth VCF file (e.g., GIAB)
    #[arg(long, help_heading = cli::sections::INPUT, value_hint = clap::ValueHint::FilePath)]
    truth: Option<ClioPath>,

    /// Competitor VCF file (e.g., DRAGEN or another Rastair version)
    #[arg(long, help_heading = cli::sections::INPUT, value_hint = clap::ValueHint::FilePath)]
    competitor: Option<ClioPath>,

    /// Regions to analyze (repeatable, e.g. -l chr1 -l chr2:100-200)
    #[arg(short = 'l', long = "region", help_heading = cli::sections::INPUT)]
    regions: Vec<RegionString>,

    /// Write JSON report to file
    #[arg(long = "output-json", help_heading = cli::sections::OUTPUT, value_hint = clap::ValueHint::FilePath)]
    output_json: Option<ClioPath>,

    /// Write interactive HTML report to file
    #[arg(long = "output-html", help_heading = cli::sections::OUTPUT, value_hint = clap::ValueHint::FilePath)]
    output_html: Option<ClioPath>,

    /// Enable experimental indel loading from VCF records.
    #[arg(long, default_value_t = false, help_heading = cli::sections::PROCESSING)]
    experimental_indels: bool,

    /// Number of threads
    #[arg(
        short = '@',
        long = "threads",
        env = "RASTAIR_THREADS",
        default_value_t = available_parallelism().map(|n| n.get()).unwrap_or(2).max(1),
        help_heading = cli::sections::PROCESSING
    )]
    threads: usize,
}

// ─── Core types ────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize)]
enum VariantCategory {
    CpG,
    DeNovo,
    Other,
    Insertion,
    Deletion,
}

impl std::fmt::Display for VariantCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VariantCategory::CpG => write!(f, "CpG"),
            VariantCategory::DeNovo => write!(f, "DeNovo"),
            VariantCategory::Other => write!(f, "Other"),
            VariantCategory::Insertion => write!(f, "Insertion"),
            VariantCategory::Deletion => write!(f, "Deletion"),
        }
    }
}

/// Position key for variant matching. Uses chromosome name (not rid) for cross-VCF compatibility.
/// Represents both SNVs (single-base alleles) and indels (multi-base alleles).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct FullPositionKey {
    chrom: SmolStr,
    pos: u64,
    ref_allele: SmolStr,
    alt_allele: SmolStr,
}

/// Position key for methylation beta matching.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct MethylationKey {
    chrom: SmolStr,
    pos: u64,
}

/// A beta value record read from a VCF `M5mC` FORMAT field.
struct BetaRecord {
    key: MethylationKey,
    beta: f64,
    is_cpg: bool,
    is_denovo: bool,
    has_variant: bool,
}

// ─── Output types (serde for JSON) ─────────────────────────────────────────

#[derive(Debug, serde::Serialize)]
struct VariantOverlap {
    truth_only: usize,
    predictions_only: usize,
    competitor_only: usize,
    truth_and_predictions: usize,
    truth_and_competitor: usize,
    predictions_and_competitor: usize,
    all_three: usize,
}

#[derive(Debug, serde::Serialize)]
struct VariantMetrics {
    precision: f64,
    recall: f64,
    f1: f64,
    tp: usize,
    fp: usize,
    fn_count: usize,
    fn_rate: f64,
}

/// Per-category precision: how many predictions in this category are in the truth set.
/// Recall is omitted because truth VCFs (e.g. GIAB) rarely carry CPG/CPGnovo flags,
/// so we cannot know how many truth variants belong to each Rastair category.
#[derive(Debug, serde::Serialize)]
struct CategoryPrecision {
    /// Predictions in this category found in truth
    tp: usize,
    /// Predictions in this category not found in truth
    fp: usize,
    /// Total predictions in this category
    n: usize,
    precision: f64,
}

/// Pred vs competitor set overlap for one variant category.
#[derive(Debug, serde::Serialize)]
struct CategoryOverlap {
    pred_only: usize,
    comp_only: usize,
    pred_and_comp: usize,
}

#[derive(Debug, serde::Serialize)]
struct CategoryMetrics {
    #[serde(skip_serializing_if = "Option::is_none")]
    predictions_vs_truth: Option<CategoryPrecision>,
    #[serde(skip_serializing_if = "Option::is_none")]
    competitor_vs_truth: Option<CategoryPrecision>,
    /// Pred vs competitor overlap (present when both callers have variants in this category).
    #[serde(skip_serializing_if = "Option::is_none")]
    overlap: Option<CategoryOverlap>,
}

#[derive(Debug, serde::Serialize)]
struct VariantReport {
    #[serde(skip_serializing_if = "Option::is_none")]
    overlap: Option<VariantOverlap>,
    #[serde(skip_serializing_if = "Option::is_none")]
    predictions_vs_truth: Option<VariantMetrics>,
    #[serde(skip_serializing_if = "Option::is_none")]
    competitor_vs_truth: Option<VariantMetrics>,
    by_category: FxHashMap<String, CategoryMetrics>,
}

/// Pre-binned 2D density grid for the beta correlation heatmap.
/// Row-major flat array indexed `[by * bins + bx]`, where `bx` is the predictions axis
/// and `by` is the competitor axis (both starting from beta=0).
#[derive(Debug, serde::Serialize)]
struct DensityGrid {
    bins: usize,
    max_count: u32,
    counts: Vec<u32>,
}

/// Pre-binned histogram of |Δβ| values.
#[derive(Debug, serde::Serialize)]
struct DiffHistogram {
    nbins: usize,
    max_diff: f64,
    counts: Vec<u32>,
}

#[derive(Debug, serde::Serialize)]
struct MethylationComparison {
    n_compared: usize,
    n_predictions_only: usize,
    n_competitor_only: usize,
    pearson_r: f64,
    r_squared: f64,
    mean_abs_diff: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    density: Option<DensityGrid>,
    #[serde(skip_serializing_if = "Option::is_none")]
    density_cpg: Option<DensityGrid>,
    #[serde(skip_serializing_if = "Option::is_none")]
    density_cpg_variant: Option<DensityGrid>,
    #[serde(skip_serializing_if = "Option::is_none")]
    density_cpg_no_variant: Option<DensityGrid>,
    #[serde(skip_serializing_if = "Option::is_none")]
    density_denovo: Option<DensityGrid>,
    #[serde(skip_serializing_if = "Option::is_none")]
    density_denovo_variant: Option<DensityGrid>,
    #[serde(skip_serializing_if = "Option::is_none")]
    density_denovo_no_variant: Option<DensityGrid>,
    #[serde(skip_serializing_if = "Option::is_none")]
    diff_histogram: Option<DiffHistogram>,
}

#[derive(Debug, serde::Serialize)]
struct Report {
    #[serde(skip_serializing_if = "Option::is_none")]
    variants: Option<VariantReport>,
    #[serde(skip_serializing_if = "Option::is_none")]
    methylation: Option<MethylationComparison>,
}

// ─── Entry point ───────────────────────────────────────────────────────────

#[instrument(level = "info", skip_all)]
pub fn verify(params: &VerifyParams) -> Result<()> {
    ensure!(
        params.truth.is_some() || params.competitor.is_some(),
        "At least one of --truth or --competitor must be provided"
    );

    let pred_path = params.predictions.path().to_path_buf();
    let truth_path = params.truth.as_ref().map(|t| t.path().to_path_buf());
    let comp_path = params.competitor.as_ref().map(|c| c.path().to_path_buf());
    let regions = &params.regions;
    let threads = params.threads;
    let load_betas_for_methyl = comp_path.is_some();

    let experimental_indels = params.experimental_indels;

    // Slots filled by each rayon task.
    let mut pred_variants: Result<FxHashMap<FullPositionKey, VariantCategory>> =
        Err(eyre!("pred variants not loaded"));
    let mut truth_variants: Result<Option<FxHashMap<FullPositionKey, VariantCategory>>> = Ok(None);
    let mut comp_variants: Result<Option<FxHashMap<FullPositionKey, VariantCategory>>> = Ok(None);
    let mut pred_betas: Result<Option<Vec<BetaRecord>>> = Ok(None);
    let mut comp_betas: Result<Option<Vec<BetaRecord>>> = Ok(None);

    rayon::scope(|s| {
        s.spawn(|_| {
            pred_variants = load_variants(&pred_path, regions, threads, experimental_indels)
                .wrap_err("Failed to load predictions variants");
        });
        if let Some(p) = truth_path.as_deref() {
            s.spawn(|_| {
                truth_variants = load_variants(p, regions, threads, experimental_indels)
                    .wrap_err("Failed to load truth variants")
                    .map(Some);
            });
        }
        if let Some(p) = comp_path.as_deref() {
            s.spawn(|_| {
                comp_variants = load_variants(p, regions, threads, experimental_indels)
                    .wrap_err("Failed to load competitor variants")
                    .map(Some);
            });
        }
        if load_betas_for_methyl {
            s.spawn(|_| {
                pred_betas = load_betas(&pred_path, regions, threads)
                    .wrap_err("Failed to load predictions betas")
                    .map(Some);
            });
            if let Some(p) = comp_path.as_deref() {
                s.spawn(|_| {
                    comp_betas = load_betas(p, regions, threads)
                        .wrap_err("Failed to load competitor betas")
                        .map(Some);
                });
            }
        }
    });

    let pred_variants = pred_variants?;
    let truth_variants = truth_variants?;
    let comp_variants = comp_variants?;
    let pred_betas = pred_betas?;
    let comp_betas = comp_betas?;

    info!(count = pred_variants.len(), "Loaded predictions variants");

    let variant_report =
        compute_variant_report(&pred_variants, truth_variants.as_ref(), comp_variants.as_ref());

    let methyl_comparison = match (pred_betas, comp_betas) {
        (Some(pb), Some(cb)) => Some(compare_methylation(pb, cb)),
        _ => None,
    };

    print_report(&variant_report, methyl_comparison.as_ref());

    let report = Report { variants: Some(variant_report), methylation: methyl_comparison };

    if let Some(json_path) = &params.output_json {
        let mut file = json_path.clone().create().wrap_err("Failed to create JSON output file")?;
        serde_json::to_writer_pretty(&mut file, &report).wrap_err("Failed to write JSON report")?;
        info!(path = %json_path.display(), "Wrote JSON report");
    }

    if let Some(html_path) = &params.output_html {
        let mut file = html_path.clone().create().wrap_err("Failed to create HTML output file")?;
        write_html_report(&report, &mut file)?;
        info!(path = %html_path.display(), "Wrote HTML report");
    }

    Ok(())
}

// ─── VCF loading ───────────────────────────────────────────────────────────

/// Load PASS variants from a VCF, returning position key → category.
/// With regions: uses `IndexedReader` (requires `.csi` index).
/// Without regions: streams sequentially (no index required).
#[instrument(level = "info", skip_all, fields(path = %path.display()))]
fn load_variants(
    path: &Path,
    regions: &[RegionString],
    threads: usize,
    experimental_indels: bool,
) -> Result<FxHashMap<FullPositionKey, VariantCategory>> {
    ensure!(path.exists(), "VCF file `{}` not found", path.display());

    let mut result = FxHashMap::default();

    if regions.is_empty() {
        let mut reader = bcf::Reader::from_path(path)
            .wrap_err_with(|| format!("Failed to open VCF: {}", path.display()))?;
        reader.set_threads(threads.max(2)).wrap_err("Failed to set reader threads")?;
        let header = reader.header().clone();
        for rec in reader.records() {
            match rec {
                Ok(r) => extract_variants(&r, &header, &mut result, experimental_indels),
                Err(e) => warn!(error = %e, "Failed to read VCF record"),
            }
        }
    } else {
        ensure_index_exists(path)?;
        let mut reader = bcf::IndexedReader::from_path(path)
            .wrap_err_with(|| format!("Failed to open indexed VCF: {}", path.display()))?;
        reader.set_threads(threads.max(2)).wrap_err("Failed to set reader threads")?;
        for region in regions {
            let header = reader.header().clone();
            let rid = header
                .name2rid(region.chromosome.as_bytes())
                .wrap_err_with(|| format!("Chromosome `{}` not found in VCF", region.chromosome))?;
            reader
                .fetch(
                    rid,
                    region.start.map(NonZeroU64::from).map(|x| x.get()).unwrap_or(0),
                    region.end.map(NonZeroU64::from).map(|x| x.get()),
                )
                .wrap_err_with(|| format!("Failed to fetch region {region}"))?;
            for rec in reader.records() {
                match rec {
                    Ok(r) => extract_variants(&r, &header, &mut result, experimental_indels),
                    Err(e) => warn!(error = %e, "Failed to read VCF record"),
                }
            }
        }
    }

    Ok(result)
}

fn extract_variants(
    record: &bcf::Record,
    header: &HeaderView,
    result: &mut FxHashMap<FullPositionKey, VariantCategory>,
    experimental_indels: bool,
) {
    if !record.has_filter("PASS".as_bytes()) {
        return;
    }

    let alleles = record.alleles();
    let ref_allele = match alleles.first() {
        Some(a) if !a.is_empty() => a,
        _ => return,
    };

    // Skip records where REF contains non-ACGT bases (e.g. N-masked regions).
    if ref_allele
        .iter()
        .any(|&b| !b.is_ascii_alphabetic() || !matches!(b, b'A' | b'C' | b'G' | b'T'))
    {
        return;
    }

    let ref_str = match std::str::from_utf8(ref_allele) {
        Ok(s) => SmolStr::from(s),
        Err(_) => return,
    };

    let chrom = match record.rid().and_then(|rid| header.rid2name(rid).ok()) {
        Some(name) => SmolStr::from(std::str::from_utf8(name).unwrap_or("unknown")),
        None => return,
    };

    let pos = record.pos() as u64;
    let is_cpg = record.info(InCpG::ID.as_bytes()).flag().unwrap_or(false);
    let is_denovo = record.info(DeNovoCpGCandidate::ID.as_bytes()).flag().unwrap_or(false);

    let snv_category = if is_denovo {
        VariantCategory::DeNovo
    } else if is_cpg {
        VariantCategory::CpG
    } else {
        VariantCategory::Other
    };

    for alt_allele in alleles.iter().skip(1) {
        if alt_allele.is_empty() {
            continue;
        }
        // Skip alt alleles with non-ACGT bases (e.g. symbolic like <*>).
        if alt_allele
            .iter()
            .any(|&b| !b.is_ascii_alphabetic() || !matches!(b, b'A' | b'C' | b'G' | b'T'))
        {
            continue;
        }
        let alt_str = match std::str::from_utf8(alt_allele) {
            Ok(s) => SmolStr::from(s),
            Err(_) => continue,
        };

        // When indels are not enabled, skip multi-base alleles entirely
        // (preserves the pre-indel behavior of only matching SNVs).
        if !experimental_indels && (ref_allele.len() != 1 || alt_allele.len() != 1) {
            continue;
        }

        // Classify by allele lengths (VCF anchor-base convention: first base is shared).
        let category = if ref_allele.len() == 1 && alt_allele.len() == 1 {
            // Single-base substitution → use INFO flag-based category.
            snv_category
        } else if ref_allele.len() == 1 && alt_allele.len() > 1 {
            // REF shorter than ALT → insertion after anchor base.
            VariantCategory::Insertion
        } else if ref_allele.len() > 1 && alt_allele.len() == 1 {
            // REF longer than ALT → deletion after anchor base.
            VariantCategory::Deletion
        } else {
            // Both multi-base (MNP or complex): treat as Other.
            VariantCategory::Other
        };

        result.insert(
            FullPositionKey {
                chrom: chrom.clone(),
                pos,
                ref_allele: ref_str.clone(),
                alt_allele: alt_str,
            },
            category,
        );
    }
}

/// Load `M5mC` beta values from a VCF.
#[instrument(level = "info", skip_all, fields(path = %path.display()))]
fn load_betas(path: &Path, regions: &[RegionString], threads: usize) -> Result<Vec<BetaRecord>> {
    ensure!(path.exists(), "VCF file `{}` not found", path.display());

    let mut result = Vec::new();

    if regions.is_empty() {
        let mut reader = bcf::Reader::from_path(path)
            .wrap_err_with(|| format!("Failed to open VCF: {}", path.display()))?;
        reader.set_threads(threads.max(2)).wrap_err("Failed to set reader threads")?;
        let header = reader.header().clone();
        for rec in reader.records() {
            match rec {
                Ok(r) => extract_beta(&r, &header, &mut result),
                Err(e) => warn!(error = %e, "Failed to read VCF record"),
            }
        }
    } else {
        ensure_index_exists(path)?;
        let mut reader = bcf::IndexedReader::from_path(path)
            .wrap_err_with(|| format!("Failed to open indexed VCF: {}", path.display()))?;
        reader.set_threads(threads.max(2)).wrap_err("Failed to set reader threads")?;
        for region in regions {
            let header = reader.header().clone();
            let rid = header
                .name2rid(region.chromosome.as_bytes())
                .wrap_err_with(|| format!("Chromosome `{}` not found in VCF", region.chromosome))?;
            reader
                .fetch(
                    rid,
                    region.start.map(NonZeroU64::from).map(|x| x.get()).unwrap_or(0),
                    region.end.map(NonZeroU64::from).map(|x| x.get()),
                )
                .wrap_err_with(|| format!("Failed to fetch region {region}"))?;
            for rec in reader.records() {
                match rec {
                    Ok(r) => extract_beta(&r, &header, &mut result),
                    Err(e) => warn!(error = %e, "Failed to read VCF record"),
                }
            }
        }
    }

    Ok(result)
}

fn extract_beta(record: &bcf::Record, header: &HeaderView, result: &mut Vec<BetaRecord>) {
    let beta = match record.format(Methylated::ID.as_bytes()).float() {
        Ok(v) => match v.first().and_then(|s| s.first().copied()) {
            Some(b) if !b.is_nan() => f64::from(b),
            _ => return,
        },
        Err(_) => return,
    };

    let chrom = match record.rid().and_then(|rid| header.rid2name(rid).ok()) {
        Some(name) => SmolStr::from(std::str::from_utf8(name).unwrap_or("unknown")),
        None => return,
    };

    let pos = record.pos() as u64;
    let is_cpg = record.info(InCpG::ID.as_bytes()).flag().unwrap_or(false);
    let is_denovo = record.info(DeNovoCpGCandidate::ID.as_bytes()).flag().unwrap_or(false);
    let has_variant =
        record.has_filter("PASS".as_bytes()) && record.alleles().iter().skip(1).any(|a| *a != b".");
    result.push(BetaRecord {
        key: MethylationKey { chrom, pos },
        beta,
        is_cpg,
        is_denovo,
        has_variant,
    });
}

fn ensure_index_exists(path: &Path) -> Result<()> {
    let csi = PathBuf::from(format!("{}.csi", path.display()));
    ensure!(
        csi.exists(),
        "VCF index not found: `{}`. Create with `bcftools index {}`",
        csi.display(),
        path.display()
    );
    Ok(())
}

// ─── Comparison logic ──────────────────────────────────────────────────────

fn compute_variant_report(
    pred: &FxHashMap<FullPositionKey, VariantCategory>,
    truth: Option<&FxHashMap<FullPositionKey, VariantCategory>>,
    competitor: Option<&FxHashMap<FullPositionKey, VariantCategory>>,
) -> VariantReport {
    let pred_keys: FxHashSet<FullPositionKey> = pred.keys().cloned().collect();
    let truth_keys: Option<FxHashSet<FullPositionKey>> = truth.map(|m| m.keys().cloned().collect());
    let comp_keys: Option<FxHashSet<FullPositionKey>> =
        competitor.map(|m| m.keys().cloned().collect());

    let overlap = compute_overlap(&pred_keys, truth_keys.as_ref(), comp_keys.as_ref());

    let pred_vs_truth = truth_keys.as_ref().map(|t| compute_metrics(&pred_keys, t));
    let comp_vs_truth = match (truth_keys.as_ref(), comp_keys.as_ref()) {
        (Some(t), Some(c)) => Some(compute_metrics(c, t)),
        _ => None,
    };

    let pred_by_cat = split_by_category(pred);
    let comp_by_cat = competitor.map(split_by_category);

    // Per-category precision: compare each category's predictions against the *full* truth set.
    // We cannot compute recall because GIAB-style truth VCFs carry no CPG/CPGnovo flags.
    let mut by_category: FxHashMap<String, CategoryMetrics> = FxHashMap::default();
    for cat in [
        VariantCategory::CpG,
        VariantCategory::DeNovo,
        VariantCategory::Other,
        VariantCategory::Insertion,
        VariantCategory::Deletion,
    ] {
        let pred_cat = pred_by_cat.get(&cat).cloned().unwrap_or_default();
        let comp_cat = comp_by_cat.as_ref().and_then(|m| m.get(&cat)).cloned().unwrap_or_default();

        let pred_vs_truth_cat = truth_keys.as_ref().and_then(|t| {
            if pred_cat.is_empty() { None } else { Some(compute_category_precision(&pred_cat, t)) }
        });
        let comp_vs_truth_cat = truth_keys.as_ref().and_then(|t| {
            if comp_cat.is_empty() { None } else { Some(compute_category_precision(&comp_cat, t)) }
        });
        let cat_overlap = if !pred_cat.is_empty() && !comp_cat.is_empty() {
            Some(CategoryOverlap {
                pred_only: pred_cat.difference(&comp_cat).count(),
                comp_only: comp_cat.difference(&pred_cat).count(),
                pred_and_comp: pred_cat.intersection(&comp_cat).count(),
            })
        } else {
            None
        };

        if pred_cat.is_empty() && comp_cat.is_empty() {
            continue;
        }
        if pred_vs_truth_cat.is_some() || comp_vs_truth_cat.is_some() || cat_overlap.is_some() {
            by_category.insert(
                cat.to_string(),
                CategoryMetrics {
                    predictions_vs_truth: pred_vs_truth_cat,
                    competitor_vs_truth: comp_vs_truth_cat,
                    overlap: cat_overlap,
                },
            );
        }
    }

    VariantReport {
        overlap,
        predictions_vs_truth: pred_vs_truth,
        competitor_vs_truth: comp_vs_truth,
        by_category,
    }
}

fn compute_overlap(
    pred: &FxHashSet<FullPositionKey>,
    truth: Option<&FxHashSet<FullPositionKey>>,
    competitor: Option<&FxHashSet<FullPositionKey>>,
) -> Option<VariantOverlap> {
    match (truth, competitor) {
        (None, None) => None,
        (Some(t), None) => Some(VariantOverlap {
            truth_only: t.difference(pred).count(),
            predictions_only: pred.difference(t).count(),
            competitor_only: 0,
            truth_and_predictions: t.intersection(pred).count(),
            truth_and_competitor: 0,
            predictions_and_competitor: 0,
            all_three: 0,
        }),
        (None, Some(c)) => Some(VariantOverlap {
            truth_only: 0,
            predictions_only: pred.difference(c).count(),
            competitor_only: c.difference(pred).count(),
            truth_and_predictions: 0,
            truth_and_competitor: 0,
            predictions_and_competitor: pred.intersection(c).count(),
            all_three: 0,
        }),
        (Some(t), Some(c)) => {
            let all: FxHashSet<&FullPositionKey> =
                t.iter().chain(pred.iter()).chain(c.iter()).collect();

            let mut truth_only = 0;
            let mut predictions_only = 0;
            let mut competitor_only = 0;
            let mut truth_and_predictions = 0;
            let mut truth_and_competitor = 0;
            let mut predictions_and_competitor = 0;
            let mut all_three = 0;

            for k in &all {
                match (t.contains(*k), pred.contains(*k), c.contains(*k)) {
                    (true, false, false) => truth_only += 1,
                    (false, true, false) => predictions_only += 1,
                    (false, false, true) => competitor_only += 1,
                    (true, true, false) => truth_and_predictions += 1,
                    (true, false, true) => truth_and_competitor += 1,
                    (false, true, true) => predictions_and_competitor += 1,
                    (true, true, true) => all_three += 1,
                    (false, false, false) => {}
                }
            }

            Some(VariantOverlap {
                truth_only,
                predictions_only,
                competitor_only,
                truth_and_predictions,
                truth_and_competitor,
                predictions_and_competitor,
                all_three,
            })
        }
    }
}

fn compute_metrics(
    caller: &FxHashSet<FullPositionKey>,
    truth: &FxHashSet<FullPositionKey>,
) -> VariantMetrics {
    let tp = caller.intersection(truth).count();
    let fp = caller.difference(truth).count();
    let fn_count = truth.difference(caller).count();
    let precision = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
    let recall = if tp + fn_count > 0 { tp as f64 / (tp + fn_count) as f64 } else { 0.0 };
    let f1 = if precision + recall > 0.0 {
        2.0 * precision * recall / (precision + recall)
    } else {
        0.0
    };
    let fn_rate = if tp + fn_count > 0 { fn_count as f64 / (tp + fn_count) as f64 } else { 0.0 };
    VariantMetrics { precision, recall, f1, tp, fp, fn_count, fn_rate }
}

/// Compute precision for a single category against the full truth set.
/// We check how many of the caller's variants in this category appear anywhere in truth.
fn compute_category_precision(
    caller_cat: &FxHashSet<FullPositionKey>,
    truth_all: &FxHashSet<FullPositionKey>,
) -> CategoryPrecision {
    let n = caller_cat.len();
    let tp = caller_cat.intersection(truth_all).count();
    let fp = n - tp;
    let precision = if n > 0 { tp as f64 / n as f64 } else { 0.0 };
    CategoryPrecision { tp, fp, n, precision }
}

fn split_by_category(
    map: &FxHashMap<FullPositionKey, VariantCategory>,
) -> FxHashMap<VariantCategory, FxHashSet<FullPositionKey>> {
    let mut result: FxHashMap<VariantCategory, FxHashSet<FullPositionKey>> = FxHashMap::default();
    for (key, &cat) in map {
        result.entry(cat).or_default().insert(key.clone());
    }
    result
}

const DENSITY_BINS: usize = 200;
const HIST_BINS: usize = 80;

fn compare_methylation(
    pred_betas: Vec<BetaRecord>,
    comp_betas: Vec<BetaRecord>,
) -> MethylationComparison {
    let comp_map: FxHashMap<&MethylationKey, &BetaRecord> =
        comp_betas.iter().map(|r| (&r.key, r)).collect();

    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let mut cpg_xs = Vec::new();
    let mut cpg_ys = Vec::new();
    let mut cpg_var_xs = Vec::new();
    let mut cpg_var_ys = Vec::new();
    let mut cpg_novar_xs = Vec::new();
    let mut cpg_novar_ys = Vec::new();
    let mut denovo_xs = Vec::new();
    let mut denovo_ys = Vec::new();
    let mut denovo_var_xs = Vec::new();
    let mut denovo_var_ys = Vec::new();
    let mut denovo_novar_xs = Vec::new();
    let mut denovo_novar_ys = Vec::new();
    let mut n_predictions_only = 0usize;

    for pred in &pred_betas {
        let Some(comp) = comp_map.get(&pred.key) else {
            n_predictions_only += 1;
            continue;
        };

        xs.push(pred.beta);
        ys.push(comp.beta);

        // Per-category split.
        if pred.is_denovo {
            denovo_xs.push(pred.beta);
            denovo_ys.push(comp.beta);
            if pred.has_variant {
                denovo_var_xs.push(pred.beta);
                denovo_var_ys.push(comp.beta);
            } else {
                denovo_novar_xs.push(pred.beta);
                denovo_novar_ys.push(comp.beta);
            }
        } else if pred.is_cpg {
            cpg_xs.push(pred.beta);
            cpg_ys.push(comp.beta);
            if pred.has_variant {
                cpg_var_xs.push(pred.beta);
                cpg_var_ys.push(comp.beta);
            } else {
                cpg_novar_xs.push(pred.beta);
                cpg_novar_ys.push(comp.beta);
            }
        }
    }

    let pred_keys: FxHashSet<&MethylationKey> = pred_betas.iter().map(|r| &r.key).collect();
    let n_competitor_only = comp_betas.iter().filter(|r| !pred_keys.contains(&r.key)).count();

    let n_compared = xs.len();
    let r = pearson_r(&xs, &ys);
    let mean_abs_diff = if n_compared > 0 {
        xs.iter().zip(&ys).map(|(x, y)| (x - y).abs()).sum::<f64>() / n_compared as f64
    } else {
        0.0
    };

    let (density, diff_histogram) = if n_compared > 0 {
        let density = compute_density_grid(&xs, &ys);
        let diff_histogram = compute_diff_histogram(&xs, &ys);
        (Some(density), Some(diff_histogram))
    } else {
        (None, None)
    };

    let density_cpg = (!cpg_xs.is_empty()).then(|| compute_density_grid(&cpg_xs, &cpg_ys));
    let density_cpg_variant =
        (!cpg_var_xs.is_empty()).then(|| compute_density_grid(&cpg_var_xs, &cpg_var_ys));
    let density_cpg_no_variant =
        (!cpg_novar_xs.is_empty()).then(|| compute_density_grid(&cpg_novar_xs, &cpg_novar_ys));
    let density_denovo =
        (!denovo_xs.is_empty()).then(|| compute_density_grid(&denovo_xs, &denovo_ys));
    let density_denovo_variant =
        (!denovo_var_xs.is_empty()).then(|| compute_density_grid(&denovo_var_xs, &denovo_var_ys));
    let density_denovo_no_variant = (!denovo_novar_xs.is_empty())
        .then(|| compute_density_grid(&denovo_novar_xs, &denovo_novar_ys));

    MethylationComparison {
        n_compared,
        n_predictions_only,
        n_competitor_only,
        pearson_r: r,
        r_squared: r * r,
        mean_abs_diff,
        density,
        density_cpg,
        density_cpg_variant,
        density_cpg_no_variant,
        density_denovo,
        density_denovo_variant,
        density_denovo_no_variant,
        diff_histogram,
    }
}

fn compute_density_grid(xs: &[f64], ys: &[f64]) -> DensityGrid {
    let mut counts = vec![0u32; DENSITY_BINS * DENSITY_BINS];
    for (&x, &y) in xs.iter().zip(ys) {
        let bx = ((x * DENSITY_BINS as f64) as usize).min(DENSITY_BINS - 1);
        let by = ((y * DENSITY_BINS as f64) as usize).min(DENSITY_BINS - 1);
        counts[by * DENSITY_BINS + bx] += 1;
    }
    let max_count = counts.iter().copied().max().unwrap_or(0);
    DensityGrid { bins: DENSITY_BINS, max_count, counts }
}

fn compute_diff_histogram(xs: &[f64], ys: &[f64]) -> DiffHistogram {
    let max_diff = xs.iter().zip(ys).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);

    let mut counts = vec![0u32; HIST_BINS];
    if max_diff > 0.0 {
        for (&x, &y) in xs.iter().zip(ys) {
            let d = (x - y).abs();
            let b = ((d / max_diff) * HIST_BINS as f64) as usize;
            counts[b.min(HIST_BINS - 1)] += 1;
        }
    } else {
        // All diffs are zero — put everything in the first bin
        counts[0] = xs.len() as u32;
    }

    DiffHistogram { nbins: HIST_BINS, max_diff, counts }
}

fn pearson_r(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len();
    if n == 0 {
        return 0.0;
    }
    let n_f = n as f64;
    let mean_x = xs.iter().sum::<f64>() / n_f;
    let mean_y = ys.iter().sum::<f64>() / n_f;
    let mut cov = 0.0f64;
    let mut var_x = 0.0f64;
    let mut var_y = 0.0f64;
    for (x, y) in xs.iter().zip(ys) {
        let dx = x - mean_x;
        let dy = y - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }
    let denom = (var_x * var_y).sqrt();
    if denom == 0.0 { 0.0 } else { cov / denom }
}

// ─── Output formatting ─────────────────────────────────────────────────────

fn fmt_n(n: usize) -> String {
    readable::num::Unsigned::from(n).to_string()
}

fn print_report(variant_report: &VariantReport, methyl: Option<&MethylationComparison>) {
    println!("\n=== Variant Verification ===\n");

    if let Some(overlap) = &variant_report.overlap {
        println!("Overlap:");
        if overlap.truth_and_predictions > 0
            || overlap.truth_only > 0
            || overlap.truth_and_competitor > 0
            || overlap.all_three > 0
        {
            println!("  Truth only:                 {:>15}", fmt_n(overlap.truth_only));
            println!("  Predictions only:           {:>15}", fmt_n(overlap.predictions_only));
            if overlap.competitor_only > 0
                || overlap.truth_and_competitor > 0
                || overlap.predictions_and_competitor > 0
                || overlap.all_three > 0
            {
                println!("  Competitor only:            {:>15}", fmt_n(overlap.competitor_only));
            }
            println!("  Truth ∩ Predictions:        {:>15}", fmt_n(overlap.truth_and_predictions));
            if overlap.competitor_only > 0
                || overlap.truth_and_competitor > 0
                || overlap.predictions_and_competitor > 0
                || overlap.all_three > 0
            {
                println!(
                    "  Truth ∩ Competitor:         {:>15}",
                    fmt_n(overlap.truth_and_competitor)
                );
                println!(
                    "  Predictions ∩ Competitor:   {:>15}",
                    fmt_n(overlap.predictions_and_competitor)
                );
                println!("  All three:                  {:>15}", fmt_n(overlap.all_three));
            }
        } else {
            println!("  Predictions only:           {:>15}", fmt_n(overlap.predictions_only));
            println!(
                "  Predictions ∩ Competitor:   {:>15}",
                fmt_n(overlap.predictions_and_competitor)
            );
            println!("  Competitor only:            {:>15}", fmt_n(overlap.competitor_only));
        }
        println!();
    }

    let has_competitor = variant_report.competitor_vs_truth.is_some();
    if variant_report.predictions_vs_truth.is_some() {
        if has_competitor {
            println!("{:<30} {:>12}  {:>12}", "Metrics vs truth:", "Predictions", "Competitor");
        } else {
            println!("{:<30} {:>12}", "Metrics vs truth:", "Predictions");
        }
        print_metric_row(
            "Precision",
            variant_report.predictions_vs_truth.as_ref(),
            variant_report.competitor_vs_truth.as_ref(),
            has_competitor,
            |m| m.precision,
        );
        print_metric_row(
            "Recall",
            variant_report.predictions_vs_truth.as_ref(),
            variant_report.competitor_vs_truth.as_ref(),
            has_competitor,
            |m| m.recall,
        );
        print_metric_row(
            "F1",
            variant_report.predictions_vs_truth.as_ref(),
            variant_report.competitor_vs_truth.as_ref(),
            has_competitor,
            |m| m.f1,
        );
        print_metric_row(
            "FN rate",
            variant_report.predictions_vs_truth.as_ref(),
            variant_report.competitor_vs_truth.as_ref(),
            has_competitor,
            |m| m.fn_rate,
        );
        println!();
    }

    if !variant_report.by_category.is_empty() {
        let has_comp_cat =
            variant_report.by_category.values().any(|m| m.competitor_vs_truth.is_some());
        println!("{:<10} {:>8}  {:>10}  {:>10}  {:>10}", "Category", "N", "TP", "FP", "Precision");
        for cat in [
            VariantCategory::CpG,
            VariantCategory::DeNovo,
            VariantCategory::Other,
            VariantCategory::Insertion,
            VariantCategory::Deletion,
        ] {
            let Some(cat_metrics) = variant_report.by_category.get(&cat.to_string()) else {
                continue;
            };
            match (&cat_metrics.predictions_vs_truth, &cat_metrics.competitor_vs_truth) {
                (Some(m), comp) => {
                    print!(
                        "  {cat:<8} {:>8}  {:>10}  {:>10}  {:>10.4}",
                        fmt_n(m.n),
                        fmt_n(m.tp),
                        fmt_n(m.fp),
                        m.precision
                    );
                    if has_comp_cat && let Some(c) = comp {
                        print!(
                            "   comp: n={} tp={} fp={} prec={:.4}",
                            fmt_n(c.n),
                            fmt_n(c.tp),
                            fmt_n(c.fp),
                            c.precision
                        );
                    }
                    println!();
                }
                (None, Some(c)) if has_comp_cat => {
                    println!(
                        "  {cat:<8} {:>8}  {:>10}  {:>10}  {:>10}   comp: n={} tp={} fp={} prec={:.4}",
                        "—",
                        "—",
                        "—",
                        "—",
                        fmt_n(c.n),
                        fmt_n(c.tp),
                        fmt_n(c.fp),
                        c.precision
                    );
                }
                (None, _) => {}
            }
        }
        println!();
    }

    if let Some(m) = methyl {
        println!("=== Methylation Comparison ===\n");
        println!("Positions compared:    {:>15}", fmt_n(m.n_compared));
        println!("Predictions only:      {:>15}", fmt_n(m.n_predictions_only));
        println!("Competitor only:       {:>15}", fmt_n(m.n_competitor_only));
        println!();
        println!("Pearson r:             {:>15.4}", m.pearson_r);
        println!("R²:                    {:>15.4}", m.r_squared);
        println!("Mean |Δβ|:             {:>15.4}", m.mean_abs_diff);
    }
}

fn print_metric_row(
    label: &str,
    pred: Option<&VariantMetrics>,
    comp: Option<&VariantMetrics>,
    has_competitor: bool,
    get: impl Fn(&VariantMetrics) -> f64,
) {
    let pred_val = pred.map(|m| format!("{:.4}", get(m))).unwrap_or_else(|| "N/A".to_string());
    if has_competitor {
        let comp_val = comp.map(|m| format!("{:.4}", get(m))).unwrap_or_else(|| "N/A".to_string());
        println!("  {label:<28} {pred_val:>12}  {comp_val:>12}");
    } else {
        println!("  {label:<28} {pred_val:>12}");
    }
}

// ─── HTML output ───────────────────────────────────────────────────────────

fn write_html_report(report: &Report, writer: &mut impl std::io::Write) -> Result<()> {
    let template = include_str!("verify_report.html");
    let data = serde_json::to_string(report).wrap_err("Failed to serialize report for HTML")?;
    let html = template.replace("{{DATA}}", &data);
    writer.write_all(html.as_bytes()).wrap_err("Failed to write HTML report")?;
    Ok(())
}

// ─── Unit tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn key(chrom: &str, pos: u64, ref_allele: &str, alt_allele: &str) -> FullPositionKey {
        FullPositionKey {
            chrom: SmolStr::from(chrom),
            pos,
            ref_allele: SmolStr::from(ref_allele),
            alt_allele: SmolStr::from(alt_allele),
        }
    }

    fn set_of(keys: impl IntoIterator<Item = FullPositionKey>) -> FxHashSet<FullPositionKey> {
        keys.into_iter().collect()
    }

    #[test]
    fn full_position_key_equality() {
        let k1 = key("chr1", 100, "C", "T");
        let k2 = key("chr1", 100, "C", "T");
        let k3 = key("chr2", 100, "C", "T");
        assert_eq!(k1, k2);
        assert_ne!(k1, k3);
    }

    #[test]
    fn three_way_venn_counts_are_correct() {
        let truth = set_of([
            key("chr1", 100, "C", "T"),
            key("chr1", 200, "G", "A"),
            key("chr1", 300, "A", "G"),
        ]);
        let predictions = set_of([
            key("chr1", 100, "C", "T"),
            key("chr1", 200, "G", "A"),
            key("chr1", 400, "C", "T"),
        ]);
        let competitor = set_of([
            key("chr1", 100, "C", "T"),
            key("chr1", 300, "A", "G"),
            key("chr1", 500, "G", "T"),
        ]);

        let overlap = compute_overlap(&predictions, Some(&truth), Some(&competitor)).unwrap();
        assert_eq!(overlap.all_three, 1, "pos 100");
        assert_eq!(overlap.truth_and_predictions, 1, "pos 200");
        assert_eq!(overlap.truth_and_competitor, 1, "pos 300");
        assert_eq!(overlap.predictions_only, 1, "pos 400");
        assert_eq!(overlap.competitor_only, 1, "pos 500");
        assert_eq!(overlap.truth_only, 0);
        assert_eq!(overlap.predictions_and_competitor, 0);
    }

    #[test]
    fn two_way_venn_truth_only() {
        let truth = set_of([key("chr1", 100, "C", "T"), key("chr1", 200, "G", "A")]);
        let predictions = set_of([key("chr1", 100, "C", "T"), key("chr1", 300, "C", "T")]);

        let overlap = compute_overlap(&predictions, Some(&truth), None).unwrap();
        assert_eq!(overlap.truth_and_predictions, 1);
        assert_eq!(overlap.truth_only, 1);
        assert_eq!(overlap.predictions_only, 1);
        assert_eq!(overlap.all_three, 0);
    }

    #[test]
    fn metrics_perfect_recall() {
        let truth = set_of([key("chr1", 100, "C", "T"), key("chr1", 200, "G", "A")]);
        let caller = set_of([
            key("chr1", 100, "C", "T"),
            key("chr1", 200, "G", "A"),
            key("chr1", 300, "A", "G"),
        ]);
        let m = compute_metrics(&caller, &truth);
        assert_eq!(m.recall, 1.0);
        assert!((m.precision - 2.0 / 3.0).abs() < 1e-10);
        assert_eq!(m.fn_count, 0);
        assert_eq!(m.fp, 1);
        assert_eq!(m.tp, 2);
    }

    #[test]
    fn metrics_empty_sets_no_panic() {
        let empty = FxHashSet::default();
        let m = compute_metrics(&empty, &empty);
        assert_eq!(m.precision, 0.0);
        assert_eq!(m.recall, 0.0);
        assert_eq!(m.f1, 0.0);
    }

    #[test]
    fn pearson_r_perfect_correlation() {
        let xs = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let ys = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        assert!((pearson_r(&xs, &ys) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn pearson_r_negative_correlation() {
        let xs = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let ys = vec![1.0, 0.75, 0.5, 0.25, 0.0];
        assert!((pearson_r(&xs, &ys) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn pearson_r_zero_variance() {
        let xs = vec![0.0, 0.5, 1.0];
        let ys = vec![0.5, 0.5, 0.5];
        assert!((pearson_r(&xs, &ys)).abs() < 1e-10);
    }

    #[test]
    fn pearson_r_empty() {
        assert_eq!(pearson_r(&[], &[]), 0.0);
    }

    #[test]
    fn methylation_comparison_includes_all_categories() {
        let pred_betas = vec![
            BetaRecord {
                key: MethylationKey { chrom: SmolStr::from("chr1"), pos: 100 },
                beta: 0.8,
                is_cpg: true,
                is_denovo: false,
                has_variant: false,
            },
            BetaRecord {
                key: MethylationKey { chrom: SmolStr::from("chr1"), pos: 200 },
                beta: 0.5,
                is_cpg: false,
                is_denovo: true,
                has_variant: false,
            },
        ];
        let comp_betas = vec![
            BetaRecord {
                key: MethylationKey { chrom: SmolStr::from("chr1"), pos: 100 },
                beta: 0.75,
                is_cpg: false,
                is_denovo: false,
                has_variant: false,
            },
            BetaRecord {
                key: MethylationKey { chrom: SmolStr::from("chr1"), pos: 200 },
                beta: 0.6,
                is_cpg: false,
                is_denovo: false,
                has_variant: false,
            },
        ];
        let result = compare_methylation(pred_betas, comp_betas);
        assert_eq!(result.n_compared, 2);
        assert!(result.density_cpg.is_some());
        assert!(result.density_denovo.is_some());
    }

    #[test]
    fn methylation_overlap_counts() {
        let pred_betas = vec![
            BetaRecord {
                key: MethylationKey { chrom: SmolStr::from("chr1"), pos: 100 },
                beta: 0.8,
                is_cpg: false,
                is_denovo: false,
                has_variant: false,
            },
            BetaRecord {
                key: MethylationKey { chrom: SmolStr::from("chr1"), pos: 300 },
                beta: 0.4,
                is_cpg: false,
                is_denovo: false,
                has_variant: false,
            },
        ];
        let comp_betas = vec![
            BetaRecord {
                key: MethylationKey { chrom: SmolStr::from("chr1"), pos: 100 },
                beta: 0.75,
                is_cpg: false,
                is_denovo: false,
                has_variant: false,
            },
            BetaRecord {
                key: MethylationKey { chrom: SmolStr::from("chr1"), pos: 200 },
                beta: 0.5,
                is_cpg: false,
                is_denovo: false,
                has_variant: false,
            },
        ];
        let result = compare_methylation(pred_betas, comp_betas);
        assert_eq!(result.n_compared, 1); // only pos 100
        assert_eq!(result.n_predictions_only, 1); // pos 300
        assert_eq!(result.n_competitor_only, 1); // pos 200
    }
}
