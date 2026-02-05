use crate::{
    call::{
        ml::DEFAULT_ML_THRESHOLD,
        process::{PileupMappingParams, calculate_pileup_metrics, get_pileups},
        variant_calling::VariantCallingParams,
    },
    metrics::{
        PileupMetrics,
        ml::{
            features::FeatureCalculator,
            types::{MlFeatureSet, RastairModel},
        },
    },
    sequence::{ChunkRegion, ReaderParams, Readers, SegmentationParams},
    utils::{PileupMetricsIteratorExt, cli},
};
use biosphere::{MaxFeatures, RandomForest, RandomForestParameters};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, ContextCompat, Result, bail, ensure};
use lz4::EncoderBuilder;
use ndarray::{Array1, Array2, Axis};
use rand::prelude::*;
use rastair_types::{Base, Probability, RegionString, SmallVec};
use rayon::prelude::*;
use rust_htslib::bcf::{self, Read as _};
use std::{collections::HashSet, num::NonZeroU64, path::PathBuf, thread::available_parallelism};
use tracing::{info, instrument, trace, warn};

#[derive(Debug, clap::Args)]
pub struct TrainModelParams {
    #[command(flatten)]
    reader: ReaderParams,

    /// Path to the ground truth file (VCF) to train with
    #[arg(help_heading = cli::sections::INPUT, value_hint=clap::ValueHint::FilePath)]
    truth: ClioPath,

    /// Output directory for trained models
    #[arg(short = 'o', long = "output", default_value = "./models")]
    #[arg(help_heading = cli::sections::OUTPUT, value_hint=clap::ValueHint::FilePath)]
    output: ClioPath,

    #[command(flatten)]
    model_params: ModelParameters,

    /// ML threshold for model evaluation (used for reporting metrics)
    #[arg(long = "ml", default_value_t = DEFAULT_ML_THRESHOLD, default_missing_value = "0.8", num_args = 0..=1)]
    #[arg(help_heading = cli::sections::TRAINING)]
    ml: Probability,

    #[arg(long, default_value_t = MlFeatureSet::Standard)]
    ml_features: MlFeatureSet,

    /// Number of threads to use
    #[arg(short='@', long = "threads", env = "RASTAIR_THREADS", default_value_t = available_parallelism().map(|n|n.get()).unwrap_or(2).max(1))]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub threads: usize,
}

#[derive(Debug, clap::Args)]
struct ModelParameters {
    /// Number of trees in the random forest
    #[arg(long = "n-trees", default_value = "800")]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub n_trees: usize,

    /// Number of features to consider at each split (mtry parameter)
    #[arg(long = "max-features", default_value = "2")]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub max_features: usize,

    /// Number of positive examples (SNPs) to sample for training
    #[arg(long = "n-positive", default_value = "4000")]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub n_positive: usize,

    /// Number of negative examples (REF positions) to sample for training
    #[arg(long = "n-negative", default_value = "16000")]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub n_negative: usize,
}

/// Key for indexing positions in truth set
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct PositionKey {
    pub pos: u64,
    pub ref_base: Base,
    pub alt_base: Base,
}

/// Training data for a specific model type (CpG, denovo, or other)
struct TrainingData {
    features: Vec<Array2<f64>>,
    labels: Vec<f64>,
}

impl TrainingData {
    fn new() -> Self {
        Self { features: Vec::new(), labels: Vec::new() }
    }

    fn add_example(&mut self, features: Array2<f64>, label: f64) {
        self.features.push(features);
        self.labels.push(label);
    }

    fn merge(&mut self, other: TrainingData) {
        self.features.extend(other.features);
        self.labels.extend(other.labels);
    }

    fn len(&self) -> usize {
        self.labels.len()
    }

    fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }
}

#[instrument(level = "info", skip_all)]
pub fn train_model(params: &TrainModelParams) -> Result<()> {
    // Create output directory if it doesn't exist
    params
        .output
        .parent()
        .wrap_err("output path invalid")
        .and_then(|p| {
            std::fs::create_dir_all(p).wrap_err("Failed to create output parent directory")
        })
        .wrap_err_with(|| {
            format!("Failed to create output directory: {}", params.output.display())
        })?;

    // Load truth VCF and index variants
    info!("Loading truth VCF from {}", params.truth.display());
    let region = params.reader.region.clone().unwrap_or_else(|| RegionString {
        chromosome: "chr12".into(),
        start: None,
        end: None,
    });
    let truth_variants = load_truth_vcf(&params.truth, &region, params.threads)
        .wrap_err("Failed to load truth VCF")?;
    info!("Loaded {} true variants from truth VCF", truth_variants.len());

    // Get segments to process
    let segmentation = SegmentationParams::default();
    let regions: Vec<ChunkRegion> = params
        .reader
        .readers()
        .wrap_err("Failed to read BAM/FASTA files")?
        .segments(segmentation.segment_max_length, segmentation.segment_overlap)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();

    if regions.is_empty() {
        bail!("No segments found in BAM file");
    }

    info!("Processing {} segments to collect training data", regions.len());

    // Process segments in parallel to collect training data
    let results: Vec<(TrainingData, TrainingData, TrainingData)> =
        rayon::ThreadPoolBuilder::new()
            .thread_name(|idx| format!("training-worker-{idx}"))
            .num_threads(params.threads)
            .start_handler(|idx| trace!(idx, "Starting training worker thread"))
            .exit_handler(|idx| trace!(idx, "Closing training worker thread"))
            .build()
            .wrap_err("Failed to create thread pool for rayon")?
            .install(move || {
                thread_local! {
                    /// Readers for the BAM and FASTA files, initialized per thread to avoid
                    /// re-opening files or having a lock
                    static READERS: std::cell::RefCell<Option<Readers>> = const { std::cell::RefCell::new(None) };
                }

                regions
                    .par_iter()
                    .map(|chunk_region| {
                        // Use thread-local readers to avoid re-opening files in each thread
                        READERS.with(|local_readers| -> (TrainingData, TrainingData, TrainingData) {
                            let mut local_readers = local_readers.borrow_mut();
                            let readers = {
                                // Initialize thread-local readers first time the thread accesses them
                                if local_readers.is_none() {
                                    match params.reader.readers() {
                                        Ok(readers) => {
                                            *local_readers = Some(readers);
                                        }
                                        Err(e) => {
                                            warn!(
                                                error = format!("{e:#}"),
                                                "Failed to open readers in worker thread"
                                            );
                                            return (
                                                TrainingData::new(),
                                                TrainingData::new(),
                                                TrainingData::new(),
                                            );
                                        }
                                    }
                                }
                                match local_readers.as_mut() {
                                    Some(readers) => readers,
                                    None => {
                                        warn!("Failed to access thread-local readers");
                                        return (
                                            TrainingData::new(),
                                            TrainingData::new(),
                                            TrainingData::new(),
                                        );
                                    }
                                }
                            };

                            // Collect training data from this segment
                            match collect_training_data_from_segment(
                                chunk_region,
                                readers,
                                &truth_variants,
                                params.ml_features.get_calculator(),
                            ) {
                                Ok(data) => data,
                                Err(e) => {
                                    warn!(
                                        error = format!("{e:#}"),
                                        "Failed to collect training data from segment"
                                    );
                                    (TrainingData::new(), TrainingData::new(), TrainingData::new())
                                }
                            }
                        })
                    })
                    .collect()
            });

    // Merge all results from parallel processing
    let mut cpg_data = TrainingData::new();
    let mut denovo_data = TrainingData::new();
    let mut other_data = TrainingData::new();

    for (cpg, denovo, other) in results {
        cpg_data.merge(cpg);
        denovo_data.merge(denovo);
        other_data.merge(other);
    }

    info!(
        cpg = cpg_data.len(),
        denovo = denovo_data.len(),
        other = other_data.len(),
        "Collected training examples",
    );

    let model = RastairModel {
        feature_set: params.ml_features,
        cpg: train_and_save_model("cpg", cpg_data, params).wrap_err("Failed to train CpG model")?,
        denovo: train_and_save_model("denovo", denovo_data, params)
            .wrap_err("Failed to train de-novo CpG model")?,
        others: train_and_save_model("other", other_data, params)
            .wrap_err("Failed to train other model")?,
        cpg_platt: Default::default(),
        denovo_platt: Default::default(),
        others_platt: Default::default(),
    };

    serialize_model(&model, params.output.clone())
        .wrap_err_with(|| format!("Failed to serialize model to {}", params.output.display()))?;

    info!(path=%params.output, "Saved model");

    Ok(())
}

/// Load truth VCF and create an index of variant positions
#[instrument(level = "info", skip_all)]
pub fn load_truth_vcf(
    vcf_path: &ClioPath,
    region: &RegionString,
    threads: usize,
) -> Result<HashSet<PositionKey>> {
    ensure!(vcf_path.exists(), "Predictions VCF file `{vcf_path:?}` not found.");
    let index_path = PathBuf::from(format!("{}.csi", vcf_path.path().display()));
    ensure!(
        index_path.exists(),
        "Predictions VCF index `{index_path:?}` not found. Please create an index with `bcftools index {vcf_path}`",
    );

    let mut reader = bcf::IndexedReader::from_path(vcf_path.path())
        .wrap_err_with(|| format!("Failed to open truth VCF file: {}", vcf_path.display()))?;
    reader.set_threads(threads.max(2)).wrap_err("Failed to set threads for truth VCF reader")?;

    let mut variants = HashSet::new();
    let header = reader.header();

    reader
        .fetch(
            header.name2rid(region.chromosome.as_bytes()).wrap_err_with(|| {
                format!("Failed to get rid for chromosome {} in truth VCF", region.chromosome)
            })?,
            region.start.map(|x| NonZeroU64::from(x).get()).unwrap_or_default(),
            region.end.map(|x| NonZeroU64::from(x).get()),
        )
        .wrap_err("Failed to fetch region from truth VCF")?;

    // Read records
    for result in reader.records() {
        let record = match result {
            Ok(record) => record,
            Err(error) => {
                warn!(error=%error, "Failed to read record from truth VCF");
                continue;
            }
        };

        variants.extend(process_truth_record(&record));
    }

    Ok(variants)
}

/// Process a truth VCF record and extract variant information.
/// Returns a list of PositionKeys for each valid SNP alt allele.
/// Multi-allelic sites produce multiple keys (one per alt).
fn process_truth_record(record: &bcf::Record) -> SmallVec<PositionKey, 2> {
    // Filter: only PASS variants
    // Note: has_filter returns true if the filter is NOT PASS
    if !record.has_filter("PASS".as_bytes()) {
        return SmallVec::new();
    }

    let alleles = record.alleles();
    if alleles.is_empty() {
        return SmallVec::new();
    }

    let ref_allele = alleles.first().expect("alleles is not empty");

    // Reference must be single base (SNP)
    if ref_allele.len() != 1 {
        return SmallVec::new();
    }

    let ref_base = Base::from(ref_allele[0]);
    if ref_base == Base::Unknown {
        return SmallVec::new();
    }

    let pos = record.pos() as u64;

    // Process each alt allele (indices 1 onwards)
    alleles
        .iter()
        .skip(1)
        .filter_map(|alt_allele| {
            // Alt must be single base (SNP)
            if alt_allele.len() != 1 {
                return None;
            }
            let alt_base = Base::from(alt_allele[0]);
            if alt_base == Base::Unknown {
                return None;
            }
            Some(PositionKey { pos, ref_base, alt_base })
        })
        .collect()
}

/// Collect training data from a single segment
fn collect_training_data_from_segment(
    chunk_region: &ChunkRegion,
    readers: &mut Readers,
    truth_variants: &HashSet<PositionKey>,
    calculator: Box<dyn FeatureCalculator>,
) -> Result<(TrainingData, TrainingData, TrainingData)> {
    // Create local training data for this segment
    let mut cpg_data = TrainingData::new();
    let mut denovo_data = TrainingData::new();
    let mut other_data = TrainingData::new();
    // Build pileups
    let mapping_params = PileupMappingParams { variant_calling: VariantCallingParams::default() };
    let (segment, pileup_iter) =
        get_pileups(readers, chunk_region, &mapping_params).wrap_err("Failed to build pileups")?;

    let metrics = calculate_pileup_metrics(pileup_iter, &segment);

    // Process each position with metrics
    metrics
        .filter_map(|x: Result<PileupMetrics>| match x {
            Err(e) => {
                warn!(error = format!("{e:#}"), "Failed to calculate pileup metrics");
                None
            }
            Ok(x) => Some(x),
        })
        .map_surrounding(|before, current, after| {
            let pos = u64::from(current.pileup.pos);

            // Process each alt allele
            for alt in &current.alts {
                let ref_base = current.pileup.reference_base;
                let alt_base = alt.base;

                // Skip Unknown bases
                if ref_base == Base::Unknown || alt_base == Base::Unknown {
                    continue;
                }

                // Determine label: is this position in truth set?
                let key = PositionKey { pos, ref_base, alt_base };
                let label = if truth_variants.contains(&key) {
                    1.0 // True variant
                } else {
                    0.0 // Reference position
                };

                // Create MetricsForAlt for this alternative allele
                let alt_metrics_for_ml = current.alt_metrics(alt_base);

                if let Some(alt_m) = alt_metrics_for_ml {
                    // Generate features based on position type
                    if alt_m.is_evidence_for_methylation() {
                        if let Ok(features) = calculator.calculate_cpg(&alt_m, before, after)
                            && !features.is_any_nan()
                        {
                            cpg_data.add_example(features, label);
                        }
                    } else if *alt.metrics.denovo {
                        if let Ok(features) = calculator.calculate_denovo_cpg(&alt_m, before, after)
                            && !features.is_any_nan()
                        {
                            denovo_data.add_example(features, label);
                        }
                    } else if let Ok(features) = calculator.calculate_others(&alt_m, before, after)
                        && !features.is_any_nan()
                    {
                        other_data.add_example(features, label);
                    }
                }
            }

            Ok(())
        })
        .for_each(|_x| {});

    Ok((cpg_data, denovo_data, other_data))
}

/// Train a model and save it to disk
fn train_and_save_model(
    model_name: &str,
    data: TrainingData,
    params: &TrainModelParams,
) -> Result<RandomForest> {
    ensure!(!data.is_empty(), "No training data for {model_name} model, skipping");

    info!("Training {} model with {} examples", model_name, data.len());

    // Subsample data
    let (features, labels) = subsample_training_data(
        data,
        params.model_params.n_positive,
        params.model_params.n_negative,
    )?;

    info!(
        "Subsampled to {} examples ({} positive, {} negative)",
        labels.len(),
        labels.iter().filter(|&&l| l == 1.0).count(),
        labels.iter().filter(|&&l| l == 0.0).count()
    );

    // Train model
    let rf_params = RandomForestParameters::default()
        .with_max_features(MaxFeatures::Value(params.model_params.max_features))
        .with_n_estimators(params.model_params.n_trees)
        .with_max_depth(None)
        .with_n_jobs(i32::try_from(params.threads).ok());

    let mut model = RandomForest::new(rf_params);
    model.fit(&features.view(), &labels.view());

    info!("Finished training {} model", model_name);

    Ok(model)
}

/// Subsample training data to balance positive and negative examples
fn subsample_training_data(
    data: TrainingData,
    n_positive: usize,
    n_negative: usize,
) -> Result<(Array2<f64>, Array1<f64>)> {
    let mut rng = rand::rng();

    // Separate positive and negative indices
    let mut positive_indices = Vec::new();
    let mut negative_indices = Vec::new();

    for (i, &label) in data.labels.iter().enumerate() {
        if label == 1.0 {
            positive_indices.push(i);
        } else {
            negative_indices.push(i);
        }
    }

    // Sample indices
    let n_pos_actual = positive_indices.len().min(n_positive);
    let n_neg_actual = negative_indices.len().min(n_negative);

    ensure!(n_pos_actual > 0, "No positive examples available for training");
    ensure!(n_neg_actual > 0, "No negative examples available for training");

    positive_indices.shuffle(&mut rng);
    negative_indices.shuffle(&mut rng);

    let selected_pos = &positive_indices[..n_pos_actual];
    let selected_neg = &negative_indices[..n_neg_actual];

    // Combine and sort indices
    let mut all_indices: Vec<usize> =
        selected_pos.iter().chain(selected_neg.iter()).copied().collect();
    all_indices.sort_unstable();

    // Extract features and labels for selected indices
    let mut feature_rows = Vec::with_capacity(all_indices.len());
    let mut label_vec = Vec::with_capacity(all_indices.len());

    for &idx in &all_indices {
        feature_rows.push(data.features[idx].row(0).to_owned());
        label_vec.push(data.labels[idx]);
    }

    // Stack feature rows into a single matrix
    let feature_views: Vec<_> = feature_rows.iter().map(|r| r.view()).collect();
    let features =
        ndarray::stack(Axis(0), &feature_views).wrap_err("Failed to stack feature arrays")?;

    let labels = Array1::from_vec(label_vec);

    Ok((features, labels))
}

/// Serialize a model to disk with LZ4 compression
fn serialize_model(model: &RastairModel, path: ClioPath) -> Result<()> {
    let file = path.create().wrap_err("Failed to create output file for model serialization")?;
    let mut encoder =
        EncoderBuilder::new().level(16).build(file).wrap_err("Failed to create LZ4 encoder")?;

    rmp_serde::encode::write(&mut encoder, &model).wrap_err("Failed to serialize model")?;

    let (_output, result) = encoder.finish();
    result.wrap_err("Failed to finalize LZ4 compression")?;

    Ok(())
}
