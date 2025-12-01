#![allow(clippy::print_stdout, reason = "we're not writing a file to stdout here")]

use crate::{
    call::ml::DEFAULT_ML_THRESHOLD,
    train::{PositionKey, load_truth_vcf},
    utils::{IntoF64, cli},
    vcf::{DeNovoCpGCandidate, InCpG},
};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, Result, ensure};
use rastair_types::{Base, Probability, RegionString};
use rastair_vcf::VcfField as _;
use rayon::prelude::*;
use rust_htslib::bcf::{self, Read as _};
use rustc_hash::{FxHashMap, FxHashSet};
use std::{
    num::NonZeroU64,
    ops::{Add, AddAssign},
    path::PathBuf,
    thread::available_parallelism,
};
use tracing::{info, instrument, warn};

#[derive(Debug, clap::Args)]
pub struct VerifyParams {
    /// Path to predictions VCF file (output from rastair2 call)
    #[arg(help_heading = cli::sections::INPUT, value_hint=clap::ValueHint::FilePath)]
    predictions: ClioPath,

    /// Path to ground truth VCF file
    #[arg(help_heading = cli::sections::INPUT, value_hint=clap::ValueHint::FilePath)]
    truth: ClioPath,

    /// Region to analyze (e.g., chr16)
    #[arg(short = 'l', long = "region")]
    #[arg(help_heading = cli::sections::INPUT)]
    region: Option<RegionString>,

    /// ML score threshold for calling variants
    #[arg(long = "ml-threshold", default_value_t = DEFAULT_ML_THRESHOLD)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    ml_threshold: Probability,

    /// Number of threads to use
    #[arg(short='@', long = "threads", env = "RASTAIR_THREADS", default_value_t = available_parallelism().map(|n|n.get()).unwrap_or(2).max(1))]
    #[arg(help_heading = cli::sections::PROCESSING)]
    threads: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum VariantCategory {
    CpG,
    DeNovo,
    Other,
}

impl std::fmt::Display for VariantCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VariantCategory::CpG => write!(f, "CpG"),
            VariantCategory::DeNovo => write!(f, "De-novo CpG"),
            VariantCategory::Other => write!(f, "Other"),
        }
    }
}

#[derive(Debug)]
struct PredictionRecord {
    key: PositionKey,
    score: f64,
    category: VariantCategory,
}

#[derive(Debug, Clone, Default)]
struct ConfusionMatrix {
    true_positives: usize,
    true_negatives: usize,
    false_positives: usize,
    false_negatives: usize,
}

impl Add for ConfusionMatrix {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            true_positives: self.true_positives + other.true_positives,
            true_negatives: self.true_negatives + other.true_negatives,
            false_positives: self.false_positives + other.false_positives,
            false_negatives: self.false_negatives + other.false_negatives,
        }
    }
}

impl AddAssign for ConfusionMatrix {
    fn add_assign(&mut self, other: Self) {
        self.true_positives += other.true_positives;
        self.true_negatives += other.true_negatives;
        self.false_positives += other.false_positives;
        self.false_negatives += other.false_negatives;
    }
}

impl ConfusionMatrix {
    fn add_result(&mut self, is_variant: bool, predicted_variant: bool) {
        match (is_variant, predicted_variant) {
            (true, true) => self.true_positives += 1,
            (true, false) => self.false_negatives += 1,
            (false, true) => self.false_positives += 1,
            (false, false) => self.true_negatives += 1,
        }
    }

    fn total(&self) -> usize {
        self.true_positives + self.true_negatives + self.false_positives + self.false_negatives
    }

    fn accuracy(&self) -> f64 {
        let total = self.total();
        if total == 0 {
            return 0.0;
        }
        (self.true_positives + self.true_negatives) as f64 / total as f64
    }

    fn sensitivity(&self) -> f64 {
        let positives = self.true_positives + self.false_negatives;
        if positives == 0 {
            return 0.0;
        }
        self.true_positives as f64 / positives as f64
    }

    fn specificity(&self) -> f64 {
        let negatives = self.true_negatives + self.false_positives;
        if negatives == 0 {
            return 0.0;
        }
        self.true_negatives as f64 / negatives as f64
    }

    fn precision(&self) -> f64 {
        let predicted_positive = self.true_positives + self.false_positives;
        if predicted_positive == 0 {
            return 0.0;
        }
        self.true_positives as f64 / predicted_positive as f64
    }

    fn negative_predictive_value(&self) -> f64 {
        let predicted_negative = self.true_negatives + self.false_negatives;
        if predicted_negative == 0 {
            return 0.0;
        }
        self.true_negatives as f64 / predicted_negative as f64
    }

    fn f1_score(&self) -> f64 {
        let precision = self.precision();
        let sensitivity = self.sensitivity();
        if precision + sensitivity == 0.0 {
            return 0.0;
        }
        2.0 * (precision * sensitivity) / (precision + sensitivity)
    }
}

#[instrument(level = "info", skip_all)]
pub fn verify(params: &VerifyParams) -> Result<()> {
    // Determine region
    let region = params.region.clone().unwrap_or_else(|| RegionString {
        chromosome: "chr1".into(),
        start: None,
        end: None,
    });

    info!("Loading truth and predictions VCFs in parallel...");
    let (truth_result, predictions_result) = rayon::join(
        || load_truth_vcf(&params.truth, &region, params.threads),
        || load_predictions_vcf(&params.predictions, &region, params.threads),
    );

    let truth_variants = truth_result.wrap_err("Failed to load truth VCF")?;
    let predictions = predictions_result.wrap_err("Failed to load predictions VCF")?;

    info!("Loaded {} true variants from truth VCF", truth_variants.len());
    info!("Loaded {} prediction records", predictions.len());

    // Calculate confusion matrices by category in parallel
    let (matrices, overall) = predictions
        .par_iter()
        .fold(
            || {
                (
                    FxHashMap::<VariantCategory, ConfusionMatrix>::default(),
                    ConfusionMatrix::default(),
                )
            },
            |(mut matrices, mut overall), pred| {
                let is_true_variant = truth_variants.contains(&pred.key);
                let is_predicted_variant = pred.score >= *params.ml_threshold;

                // Update category-specific matrix
                matrices
                    .entry(pred.category)
                    .or_default()
                    .add_result(is_true_variant, is_predicted_variant);

                // Update overall matrix
                overall.add_result(is_true_variant, is_predicted_variant);

                (matrices, overall)
            },
        )
        .reduce(
            || {
                (
                    FxHashMap::<VariantCategory, ConfusionMatrix>::default(),
                    ConfusionMatrix::default(),
                )
            },
            |(mut m1, o1), (m2, o2)| {
                // Merge hashmaps
                for (cat, matrix) in m2 {
                    *m1.entry(cat).or_default() += matrix;
                }
                (m1, o1 + o2)
            },
        );

    // Build FxHashSet of prediction keys for fast membership testing
    let prediction_keys: FxHashSet<_> = predictions.iter().map(|p| &p.key).collect();

    // Calculate how many truth variants were not in predictions at all (parallel)
    let uncallable = truth_variants.par_iter().filter(|key| !prediction_keys.contains(key)).count();

    // Output results
    println!("\n=== Model Verification Results ===\n");
    println!("Region: {}", region);
    println!("ML Threshold: {:.2}", params.ml_threshold);
    println!("True variants in truth set: {}", truth_variants.len());
    println!("Positions with evidence in predictions: {}", predictions.len());
    println!("Truth variants without evidence (uncallable): {}\n", uncallable);

    println!("=== Overall Performance ===\n");
    print_metrics(&overall);

    println!("\n=== Performance by Variant Category ===\n");
    for category in [VariantCategory::CpG, VariantCategory::DeNovo, VariantCategory::Other] {
        if let Some(matrix) = matrices.get(&category)
            && matrix.total() > 0
        {
            println!("--- {} ---", category);
            print_metrics(matrix);
            println!();
        }
    }

    // Additional breakdown: true variants by category
    println!("=== True Variants Called by Category ===\n");
    for category in [VariantCategory::CpG, VariantCategory::DeNovo, VariantCategory::Other] {
        if let Some(matrix) = matrices.get(&category) {
            let total_true = matrix.true_positives + matrix.false_negatives;
            if total_true > 0 {
                println!(
                    "{}: {}/{} called ({:.2}%)",
                    category,
                    matrix.true_positives,
                    total_true,
                    100.0 * matrix.true_positives as f64 / total_true as f64
                );
            }
        }
    }

    println!("\n=== False Positives by Category ===\n");
    for category in [VariantCategory::CpG, VariantCategory::DeNovo, VariantCategory::Other] {
        if let Some(matrix) = matrices.get(&category) {
            let total_ref = matrix.true_negatives + matrix.false_positives;
            if total_ref > 0 {
                println!(
                    "{}: {}/{} false positives ({:.4}%)",
                    category,
                    matrix.false_positives,
                    total_ref,
                    100.0 * matrix.false_positives as f64 / total_ref as f64
                );
            }
        }
    }

    Ok(())
}

fn print_metrics(matrix: &ConfusionMatrix) {
    println!("Confusion Matrix:");
    println!("  True Positives:  {}", matrix.true_positives);
    println!("  True Negatives:  {}", matrix.true_negatives);
    println!("  False Positives: {}", matrix.false_positives);
    println!("  False Negatives: {}", matrix.false_negatives);
    println!();
    println!("Metrics:");
    println!("  Accuracy:    {:.4}", matrix.accuracy());
    println!("  Sensitivity: {:.4} (recall for variants)", matrix.sensitivity());
    println!("  Specificity: {:.4} (recall for reference)", matrix.specificity());
    println!("  Precision:   {:.4} (PPV)", matrix.precision());
    println!("  NPV:         {:.4}", matrix.negative_predictive_value());
    println!("  F1 Score:    {:.4}", matrix.f1_score());
}

#[instrument(level = "info", skip_all)]
fn load_predictions_vcf(
    vcf_path: &ClioPath,
    region: &RegionString,
    threads: usize,
) -> Result<Vec<PredictionRecord>> {
    // Ensure index file exists
    let index_path = PathBuf::from(format!("{}.csi", vcf_path.path().display()));
    ensure!(
        index_path.exists(),
        "Predictions VCF file `{}` not found. Please create an index with `bcftools index {}`",
        index_path.display(),
        vcf_path.display(),
    );

    let mut reader = bcf::IndexedReader::from_path(vcf_path.path())
        .wrap_err_with(|| format!("Failed to open predictions VCF file: {}", vcf_path.display()))?;
    reader
        .set_threads(threads.max(2))
        .wrap_err("Failed to set threads for predictions VCF reader")?;

    let header = reader.header().clone();

    reader
        .fetch(
            header.name2rid(region.chromosome.as_bytes()).wrap_err_with(|| {
                format!("Failed to get rid for chromosome {} in predictions VCF", region.chromosome)
            })?,
            region.start.map(|x| NonZeroU64::from(x).get()).unwrap_or_default(),
            region.end.map(|x| NonZeroU64::from(x).get()),
        )
        .wrap_err("Failed to fetch region from predictions VCF")?;

    let mut predictions = Vec::new();

    for result in reader.records() {
        let record = match result {
            Ok(record) => record,
            Err(error) => {
                warn!(error=%error, "Failed to read record from predictions VCF");
                continue;
            }
        };

        // Extract basic info
        let pos = record.pos() as u64;
        let alleles = record.alleles();

        if alleles.is_empty() {
            continue;
        }

        let ref_allele = alleles[0];
        if ref_allele.len() != 1 {
            continue; // Skip non-SNP reference
        }

        let ref_base = Base::from(ref_allele[0]);
        if ref_base == Base::Unknown {
            continue;
        }

        // Check CpG status
        let is_cpg = record.info(InCpG::ID.as_bytes()).flag().unwrap_or(false);
        let is_denovo = record.info(DeNovoCpGCandidate::ID.as_bytes()).flag().unwrap_or(false);

        // Process each alternate allele
        for (alt_idx, alt_allele) in alleles.iter().skip(1).enumerate() {
            if alt_allele.len() != 1 {
                continue; // Skip non-SNP
            }

            let alt_base = Base::from(alt_allele[0]);
            if alt_base == Base::Unknown {
                continue;
            }

            // Get ML score for this allele from FORMAT field (for first sample)
            let score = record
                .format(b"ML")
                .float()
                .ok()
                .and_then(|v| v.first().and_then(|scores| scores.get(alt_idx).copied()))
                .unwrap_or(0.0)
                .f();

            let category = if is_denovo {
                VariantCategory::DeNovo
            } else if is_cpg {
                VariantCategory::CpG
            } else {
                VariantCategory::Other
            };

            predictions.push(PredictionRecord {
                key: PositionKey { pos, ref_base, alt_base },
                score,
                category,
            });
        }
    }

    Ok(predictions)
}
