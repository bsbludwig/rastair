//! Train random forest classifiers for variant filtering.
//!
//! Rastair uses five separate RF models: CpG, de-novo CpG, "other" (SNVs),
//! insertion, and deletion.
//! Training works like this:
//!
//! 1. Collect: iterate pileups, compute ML for alts & indels, put into model buckets
//! 2. Sample: pick `n_positive` and `n_negative` examples
//! 3. Train a `RandomForest` (params from CLI) on the samples
//! 4. Scaling: do Platt scaling on everything but the sampled data
//! 5. Export: Build `RastairFlatModel` and write to file

use crate::{
    call::{
        ml::DEFAULT_ML_THRESHOLD,
        pileup::indels::IndelAllele,
        process::{PileupMappingParams, calculate_pileup_metrics, get_pileups},
        variant_calling::indel_calling::{self, IndelParams},
    },
    metrics::{
        MetricsForIndel, PileupMetrics,
        ml::{
            features::FeatureCalculator,
            types::{MlFeatureSet, PlattScaling, RastairFlatModel},
        },
    },
    rayon_all,
    sequence::{ChunkRegion, ReaderParams, Readers, SegmentationParams},
    utils::{PileupMetricsIteratorExt, cli},
};
use biosphere::{FlatForest, MaxFeatures, RandomForest, RandomForestParameters};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, ContextCompat, Result, bail, ensure};
use lz4::EncoderBuilder;
use ndarray::{Array1, Array2, Axis};
use rand::prelude::*;
use rayon::prelude::*;
use rust_htslib::bcf::{self, Read as _};
use seqair_types::{Base, Probability, RegionString, SmallVec, SmolStr};
use std::{
    collections::HashSet,
    fmt,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    thread::available_parallelism,
};
use tracing::{error, info, instrument, trace, warn};

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

    /// Export collected features and labels as TSV files to this directory.
    /// One file per model type: `cpg_features.tsv`, `denovo_features.tsv`,
    /// `other_features.tsv`, `insertion_features.tsv`, `deletion_features.tsv`.
    #[arg(long, help_heading = cli::sections::OUTPUT, value_hint=clap::ValueHint::DirPath)]
    export_features: Option<ClioPath>,

    /// Export features importances as TSV files to this directory.
    #[arg(long, help_heading = cli::sections::OUTPUT, value_hint=clap::ValueHint::DirPath)]
    feature_analytics: Option<ClioPath>,

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
    #[arg(long = "n-trees", default_value_t = 800)]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub n_trees: usize,

    /// Number of features to consider at each split (mtry parameter)
    #[arg(long = "max-features", default_value_t = 2)]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub max_features: usize,

    /// Maximum tree depth. `0` grows trees until pure (unbounded), which can
    /// produce very large models on noisy/poorly-separable data such as indels.
    /// A cap of ~20 typically removes the noise-memorising depth at negligible
    /// accuracy cost. See <https://scikit-learn.org/stable/modules/tree.html>.
    #[arg(long = "max-depth", default_value_t = 20)]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub max_depth: usize,

    /// Minimum number of samples required at each leaf. Larger values prevent
    /// the forest from memorising individual samples, regularising noisy data
    /// and shrinking the model. scikit-learn suggests 5 as a starting value;
    /// use 1 for fully-grown leaves (best for clean, well-separated classes).
    #[arg(long = "min-samples-leaf", default_value_t = 5)]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub min_samples_leaf: usize,

    /// Number of positive examples (SNPs) to sample for training
    #[arg(long = "n-positive", default_value_t = 8_000)]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub n_positive: usize,

    /// Number of negative examples (REF positions) to sample for training
    #[arg(long = "n-negative", default_value_t = 20_000)]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub n_negative: usize,

    /// Random seed for reproducibility (subsampling and forest training).
    /// Omit for a random seed.
    #[arg(long)]
    #[arg(help_heading = cli::sections::TRAINING)]
    pub seed: Option<u64>,
}

/// Key for indexing positions in truth set
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct PositionKey {
    pub pos: u64,
    pub ref_base: Base,
    pub alt_base: Base,
}

/// Key for indexing indel positions in truth set
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct IndelKey {
    pub pos: u64,
    pub allele: IndelAllele,
}

/// Training data for a specific model type
struct TrainingData {
    features: Vec<Array2<f64>>,
    labels: Vec<f64>,
    positions: Vec<(SmolStr, u64)>,
}

impl TrainingData {
    fn new() -> Self {
        Self { features: Vec::new(), labels: Vec::new(), positions: Vec::new() }
    }

    fn add_example(&mut self, features: Array2<f32>, label: f64, chrom: SmolStr, pos: u64) {
        // Features are computed in f32 (matching the f32 inference forests);
        // biosphere's RandomForest fits on f64, so widen at this boundary only.
        self.features.push(features.mapv(f64::from));
        self.labels.push(label);
        self.positions.push((chrom, pos));
    }

    fn merge(&mut self, other: TrainingData) {
        self.features.extend(other.features);
        self.labels.extend(other.labels);
        self.positions.extend(other.positions);
    }

    fn len(&self) -> usize {
        self.labels.len()
    }

    fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    fn positives(&self) -> usize {
        self.labels.iter().filter(|&&l| l == 1.0).count()
    }
}

type SegmentResult = (TrainingData, TrainingData, TrainingData, TrainingData, TrainingData);

fn empty_segment_result() -> SegmentResult {
    (
        TrainingData::new(),
        TrainingData::new(),
        TrainingData::new(),
        TrainingData::new(),
        TrainingData::new(),
    )
}

#[instrument(level = "debug", skip_all)]
pub fn train_model(params: &TrainModelParams) -> Result<()> {
    let seed = params.model_params.seed.unwrap_or_else(rand::random);
    info!(
        seed,
        n_trees = params.model_params.n_trees,
        max_features = params.model_params.max_features,
        n_positive = params.model_params.n_positive,
        n_negative = params.model_params.n_negative,
        "Training parameters",
    );

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

    // Load truth VCF and index variants across all requested regions.
    // If no regions are specified, default to the entire chr12 (common for training).
    let regions =
        params.reader.regions.as_ref().map(|input| input.regions().to_vec()).unwrap_or_else(|| {
            vec![RegionString { chromosome: "chr12".into(), start: None, end: None }]
        });

    let mut snp_variants = HashSet::new();
    let mut indel_variants = HashSet::new();
    for region in &regions {
        let (snps, indels) = load_truth_vcf(&params.truth, region, params.threads)
            .wrap_err_with(|| format!("Failed to load truth VCF for region {region}"))?;
        snp_variants.extend(snps);
        indel_variants.extend(indels);
    }

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

    // Use default IndelParams with experimental_indels enabled so we can
    // collect training examples for insertion/deletion models.
    let indel_params = IndelParams { experimental_indels: true, ..IndelParams::default() };
    let calculator = params.ml_features.get_calculator();

    // Process segments in parallel to collect training data
    let results: Vec<SegmentResult> = rayon::ThreadPoolBuilder::new()
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
                    let _span = tracing::info_span!("collect_segment", region = %chunk_region.region).entered();

                    // Use thread-local readers to avoid re-opening files in each thread
                    READERS.with(|local_readers| -> SegmentResult {
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
                                        return empty_segment_result();
                                    }
                                }
                            }
                            match local_readers.as_mut() {
                                Some(readers) => readers,
                                None => {
                                    warn!("Failed to access thread-local readers");
                                    return empty_segment_result();
                                }
                            }
                        };

                        // Collect training data from this segment
                        match collect_training_data_from_segment(
                            chunk_region,
                            readers,
                            &snp_variants,
                            &indel_variants,
                            &indel_params,
                            &*calculator,
                        ) {
                            Ok(data) => data,
                            Err(e) => {
                                warn!(
                                    error = format!("{e:#}"),
                                    "Failed to collect training data from segment"
                                );
                                empty_segment_result()
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
    let mut insertion_data = TrainingData::new();
    let mut deletion_data = TrainingData::new();

    for (cpg, denovo, other, insertion, deletion) in results {
        cpg_data.merge(cpg);
        denovo_data.merge(denovo);
        other_data.merge(other);
        insertion_data.merge(insertion);
        deletion_data.merge(deletion);
    }

    info!(
        cpg = cpg_data.len(),
        cpg_pos = cpg_data.positives(),
        denovo = denovo_data.len(),
        denovo_pos = denovo_data.positives(),
        other = other_data.len(),
        other_pos = other_data.positives(),
        insertion = insertion_data.len(),
        insertion_pos = insertion_data.positives(),
        deletion = deletion_data.len(),
        deletion_pos = deletion_data.positives(),
        "Collected training examples",
    );

    let features = params.ml_features.get_calculator().feature_num();
    let feature_names = params.ml_features.get_calculator().feature_names();

    if let Some(ref export_dir) = params.export_features {
        std::fs::create_dir_all(export_dir.path()).wrap_err_with(|| {
            format!("Failed to create export directory: {}", export_dir.display())
        })?;
        info!(dir = %export_dir.display(), "Exporting features as TSV");
        export_features_tsv(&cpg_data, "cpg", &feature_names.cpg, export_dir.path())?;
        export_features_tsv(&denovo_data, "denovo", &feature_names.denovo_cpg, export_dir.path())?;
        export_features_tsv(&other_data, "other", &feature_names.others, export_dir.path())?;
        export_features_tsv(
            &insertion_data,
            "insertion",
            &feature_names.insertion,
            export_dir.path(),
        )?;
        export_features_tsv(
            &deletion_data,
            "deletion",
            &feature_names.deletion,
            export_dir.path(),
        )?;
    }

    // Derive independent seeds for each model from the base seed
    let mut seed_rng = rand::rngs::StdRng::seed_from_u64(seed);
    let cpg_seed: u64 = seed_rng.random();
    let denovo_seed: u64 = seed_rng.random();
    let others_seed: u64 = seed_rng.random();
    let insertion_seed: u64 = seed_rng.random();
    let deletion_seed: u64 = seed_rng.random();

    info!("Training all 5 models in parallel");
    let (cpg_result, denovo_result, others_result, insertion_result, deletion_result) = rayon_all!(
        train_and_save("cpg", cpg_data, params, cpg_seed),
        train_and_save("denovo", denovo_data, params, denovo_seed),
        train_and_save("other", other_data, params, others_seed),
        train_and_save("insertion", insertion_data, params, insertion_seed),
        train_and_save("deletion", deletion_data, params, deletion_seed),
    );

    let (cpg, cpg_platt) = cpg_result.wrap_err("Failed to train CpG model")?;
    if let Some(path) = params.feature_analytics.as_ref() {
        export_feature_importances(
            &cpg,
            &feature_names.cpg,
            &path.path().join("cpg_feature_importances.csv"),
        )
        .wrap_err("Failed to export CpG feature importances")?;
    }
    let cpg = FlatForest::from_forest(&cpg, features.cpg);

    let (denovo, denovo_platt) = denovo_result.wrap_err("Failed to train de-novo CpG model")?;
    if let Some(path) = params.feature_analytics.as_ref() {
        export_feature_importances(
            &denovo,
            &feature_names.denovo_cpg,
            &path.path().join("denovo_feature_importances.csv"),
        )
        .wrap_err("Failed to export de-novo CpG feature importances")?;
    }
    let denovo = FlatForest::from_forest(&denovo, features.denovo_cpg);

    let (others, others_platt) = others_result.wrap_err("Failed to train other model")?;
    if let Some(path) = params.feature_analytics.as_ref() {
        export_feature_importances(
            &others,
            &feature_names.others,
            &path.path().join("other_feature_importances.csv"),
        )
        .wrap_err("Failed to export other feature importances")?;
    }
    let others = FlatForest::from_forest(&others, features.others);

    let (insertion, insertion_platt) =
        insertion_result.wrap_err("Failed to train insertion model")?;
    if let Some(path) = params.feature_analytics.as_ref() {
        export_feature_importances(
            &insertion,
            &feature_names.insertion,
            &path.path().join("insertion_feature_importances.csv"),
        )
        .wrap_err("Failed to export insertion feature importances")?;
    }
    let insertion = FlatForest::from_forest(&insertion, features.insertion);

    let (deletion, deletion_platt) = deletion_result.wrap_err("Failed to train deletion model")?;
    if let Some(path) = params.feature_analytics.as_ref() {
        export_feature_importances(
            &deletion,
            &feature_names.deletion,
            &path.path().join("deletion_feature_importances.csv"),
        )
        .wrap_err("Failed to export deletion feature importances")?;
    }
    let deletion = FlatForest::from_forest(&deletion, features.deletion);

    #[derive(Debug, serde::Serialize)]
    struct ModelReport<'a> {
        forest: &'a biosphere::ForestMeta,
        scaling: PlattScaling,
    }

    let report = std::collections::BTreeMap::from([
        ("cpg", ModelReport { forest: &cpg.meta, scaling: cpg_platt }),
        ("denovo", ModelReport { forest: &denovo.meta, scaling: denovo_platt }),
        ("others", ModelReport { forest: &others.meta, scaling: others_platt }),
        ("insertion", ModelReport { forest: &insertion.meta, scaling: insertion_platt }),
        ("deletion", ModelReport { forest: &deletion.meta, scaling: deletion_platt }),
    ]);
    info!(report = ?report, "Trained models");

    let model = RastairFlatModel {
        feature_set: params.ml_features,
        cpg,
        cpg_platt,
        denovo,
        denovo_platt,
        others,
        others_platt,
        insertion,
        insertion_platt,
        deletion,
        deletion_platt,
    };

    serialize_model(&model, params.output.clone())
        .wrap_err_with(|| format!("Failed to serialize model to {}", params.output.display()))?;

    info!(path=%params.output, "Saved model");

    Ok(())
}

/// Load truth VCF and create an index of variant positions (SNPs and indels).
#[instrument(level = "info", skip_all)]
pub fn load_truth_vcf(
    vcf_path: &ClioPath,
    region: &RegionString,
    threads: usize,
) -> Result<(HashSet<PositionKey>, HashSet<IndelKey>)> {
    info!(path=%vcf_path, %region, "Loading truth vcf");

    ensure!(vcf_path.exists(), "Predictions VCF file `{vcf_path:?}` not found.");
    let index_path = PathBuf::from(format!("{}.csi", vcf_path.path().display()));
    ensure!(
        index_path.exists(),
        "Predictions VCF index `{index_path:?}` not found. Please create an index with `bcftools index {vcf_path}`",
    );

    let mut reader = bcf::IndexedReader::from_path(vcf_path.path())
        .wrap_err_with(|| format!("Failed to open truth VCF file: {}", vcf_path.display()))?;
    reader.set_threads(threads.max(2)).wrap_err("Failed to set threads for truth VCF reader")?;

    let mut snp_variants = HashSet::new();
    let mut indel_variants = HashSet::new();
    let header = reader.header();

    reader
        .fetch(
            header.name2rid(region.chromosome.as_bytes()).wrap_err_with(|| {
                format!("Failed to get rid for chromosome {} in truth VCF", region.chromosome)
            })?,
            region.start.map(|x: seqair_types::Pos1| x.as_u64()).unwrap_or_default(),
            region.end.map(|x: seqair_types::Pos1| x.as_u64()),
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

        let (snps, indels) = process_truth_record(&record);
        snp_variants.extend(snps);
        indel_variants.extend(indels);
    }

    info!(snps = snp_variants.len(), indels = indel_variants.len(), "Loaded true variants");

    Ok((snp_variants, indel_variants))
}

/// Process a truth VCF record and extract variant information.
///
/// Returns SNP [`PositionKey`]s for single-base substitutions and
/// [`IndelKey`]s for insertions/deletions.
/// Multi-allelic sites produce multiple keys (one per alt).
fn process_truth_record(record: &bcf::Record) -> (SmallVec<PositionKey, 2>, SmallVec<IndelKey, 2>) {
    // Filter: only PASS variants
    if !record.has_filter("PASS".as_bytes()) {
        return (SmallVec::new(), SmallVec::new());
    }

    let alleles = record.alleles();
    if alleles.is_empty() {
        return (SmallVec::new(), SmallVec::new());
    }

    let ref_allele = alleles.first().expect("alleles is not empty");
    let pos = record.pos() as u64;

    let mut snps = SmallVec::new();
    let mut indels = SmallVec::new();

    let ref_base = if ref_allele.is_empty() { Base::Unknown } else { Base::from(ref_allele[0]) };

    for alt_allele in alleles.iter().skip(1) {
        if alt_allele.is_empty() {
            continue;
        }

        // Both ref and alt are single base → SNP
        if ref_allele.len() == 1 && alt_allele.len() == 1 {
            let alt_base = Base::from(alt_allele[0]);
            if ref_base != Base::Unknown && alt_base != Base::Unknown {
                snps.push(PositionKey { pos, ref_base, alt_base });
            }
            continue;
        }

        // Multi-base allele → indel
        // The first base of REF and ALT must match (VCF anchor base).
        // The remainder determines whether it's an insertion or deletion.
        if ref_allele.len() < 2 && alt_allele.len() < 2 {
            // One of them is empty or just the anchor — can't determine
            continue;
        }

        // Parse the first base from each
        let ref_first = ref_allele.first().copied().unwrap_or(b'N');
        let alt_first = alt_allele.first().copied().unwrap_or(b'N');
        if ref_first != alt_first {
            // Anchor base mismatch — skip (complex variant, not a simple indel)
            continue;
        }

        let ref_rest = &ref_allele[1..];
        let alt_rest = &alt_allele[1..];

        let allele = if ref_rest.len() > alt_rest.len() {
            // Deletion: REF has extra bases
            let del_bases: SmallVec<Base, 4> = ref_rest.iter().map(|&b| Base::from(b)).collect();
            if del_bases.contains(&Base::Unknown) {
                continue;
            }
            IndelAllele::Deletion(del_bases)
        } else if alt_rest.len() > ref_rest.len() {
            // Insertion: ALT has extra bases
            let ins_bases: SmallVec<Base, 4> = alt_rest.iter().map(|&b| Base::from(b)).collect();
            if ins_bases.contains(&Base::Unknown) {
                continue;
            }
            IndelAllele::Insertion(ins_bases)
        } else {
            // Same length but multi-base (e.g. MNP) — skip for now
            continue;
        };

        indels.push(IndelKey { pos, allele });
    }

    (snps, indels)
}

/// Collect training data from a single segment
fn collect_training_data_from_segment(
    chunk_region: &ChunkRegion,
    readers: &mut Readers,
    snp_truth: &HashSet<PositionKey>,
    indel_truth: &HashSet<IndelKey>,
    indel_params: &IndelParams,
    calculator: &dyn FeatureCalculator,
) -> Result<SegmentResult> {
    // Create local training data for this segment
    let mut cpg_data = TrainingData::new();
    let mut denovo_data = TrainingData::new();
    let mut other_data = TrainingData::new();
    let mut insertion_data = TrainingData::new();
    let mut deletion_data = TrainingData::new();

    // Build pileups
    let mapping_params = PileupMappingParams::default();
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

            // -- SNP alt alleles --
            for alt in &current.alts {
                let ref_base = current.pileup.reference_base;
                let alt_base = alt.base;

                // Skip Unknown bases
                if ref_base == Base::Unknown || alt_base == Base::Unknown {
                    continue;
                }

                // Determine label: is this position in truth set?
                let key = PositionKey { pos, ref_base, alt_base };
                let label = if snp_truth.contains(&key) { 1.0 } else { 0.0 };

                // Create MetricsForAlt for this alternative allele
                let alt_metrics_for_ml = current.alt_metrics(alt_base);

                if let Some(alt_m) = alt_metrics_for_ml {
                    let chrom = chunk_region.contig.clone();
                    // Generate features based on position type
                    if alt_m.is_evidence_for_methylation() {
                        if let Ok(features) = calculator.calculate_cpg(&alt_m, before, after)
                            && !features.is_any_nan()
                        {
                            cpg_data.add_example(features, label, chrom.clone(), pos);
                        }
                    } else if *alt.metrics.denovo {
                        if let Ok(features) = calculator.calculate_denovo_cpg(&alt_m, before, after)
                            && !features.is_any_nan()
                        {
                            denovo_data.add_example(features, label, chrom.clone(), pos);
                        }
                    } else if let Ok(features) = calculator.calculate_others(&alt_m, before, after)
                        && !features.is_any_nan()
                    {
                        other_data.add_example(features, label, chrom.clone(), pos);
                    }
                }
            }

            // -- Indel alleles --
            if !current.indels.is_empty() {
                let indel_calls = indel_calling::call_indels(&current.indels, indel_params, true);
                for call in &indel_calls {
                    let indel_key = IndelKey { pos, allele: call.allele.clone() };
                    let label = if indel_truth.contains(&indel_key) { 1.0 } else { 0.0 };

                    let indel_m = MetricsForIndel { metrics: current, indel: call };

                    match &call.allele {
                        IndelAllele::Insertion(_) => {
                            if let Ok(features) = calculator.calculate_insertion(&indel_m)
                                && !features.is_any_nan()
                            {
                                insertion_data.add_example(
                                    features,
                                    label,
                                    chunk_region.contig.clone(),
                                    pos,
                                );
                            }
                        }
                        IndelAllele::Deletion(_) => {
                            if let Ok(features) = calculator.calculate_deletion(&indel_m)
                                && !features.is_any_nan()
                            {
                                deletion_data.add_example(
                                    features,
                                    label,
                                    chunk_region.contig.clone(),
                                    pos,
                                );
                            }
                        }
                    }
                }
            }

            Ok(())
        })
        .for_each(|_x| {});

    Ok((cpg_data, denovo_data, other_data, insertion_data, deletion_data))
}

#[instrument(level = "info", skip_all, fields(model=%model_name))]
fn train_and_save(
    model_name: &str,
    data: TrainingData,
    params: &TrainModelParams,
    seed: u64,
) -> Result<(RandomForest, PlattScaling)> {
    ensure!(!data.is_empty(), "No training data for {model_name} model, skipping");

    info!(seed, examples = data.len(), "Training model");

    // Subsample for training, keep held-out data for Platt calibration
    let (train_features, train_labels, holdout_features, holdout_labels) = subsample_training_data(
        &data,
        params.model_params.n_positive,
        params.model_params.n_negative,
        seed,
    )?;

    info!(
        training = train_labels.len(),
        positive = train_labels.iter().filter(|&&l| l == 1.0).count(),
        negative = train_labels.iter().filter(|&&l| l == 0.0).count(),
        holdout = holdout_labels.len(),
        "Subsampled"
    );

    // Train model. `max_depth == 0` means unbounded (grow until pure).
    let max_depth = (params.model_params.max_depth != 0).then_some(params.model_params.max_depth);
    let rf_params = RandomForestParameters::default()
        .with_max_features(MaxFeatures::Value(params.model_params.max_features))
        .with_n_estimators(params.model_params.n_trees)
        .with_max_depth(max_depth)
        .with_min_samples_leaf(params.model_params.min_samples_leaf)
        .with_n_jobs(i32::try_from(params.threads).ok())
        .with_seed(seed);

    let mut model = RandomForest::new(rf_params);
    model.fit(&train_features.view(), &train_labels.view());

    // Fit Platt scaling on held-out predictions
    let raw_scores = model.predict(&holdout_features.view());
    let platt = fit_platt_scaling(
        raw_scores.as_slice().unwrap_or(&[]),
        holdout_labels.as_slice().unwrap_or(&[]),
    );

    let raw = raw_scores.as_slice().unwrap_or(&[]);
    let labels = holdout_labels.as_slice().unwrap_or(&[]);
    if let (Some(before), Some(after)) = (
        calibration(raw, labels),
        calibration(&raw.iter().map(|&s| *platt.calibrate_score(s)).collect::<Vec<_>>(), labels),
    ) {
        // The holdout is "everything not sampled for training", so its base rate
        // is a side effect of the subsampling: training consumes `n_negative`
        // negatives against `n_positive` positives, which leaves the holdout
        // enriched in positives wherever those counts are a large share of the
        // data. That is negligible for the CpG models and material for the indel
        // ones, and it biases the fit — worth seeing next to the numbers.
        info!(%before, %after, "Calibration on holdout (raw forest votes vs Platt-scaled)");
        if after.ece > before.ece {
            warn!(
                model = model_name,
                before_ece = before.ece,
                after_ece = after.ece,
                "Platt scaling made calibration worse. A random forest's vote fraction is \
                 already a bounded probability, so a two-parameter sigmoid has little to fix \
                 and can distort it; isotonic regression is the usual alternative at this \
                 holdout size."
            );
        }
    }

    if platt.a == 1.0 && platt.b == 0.0 {
        error!(
            "Platt scaling is identity (a=1.0, b=0.0) — model likely failed to learn. \
             Check class balance, feature quality, and consider reducing n-negative."
        );
    }

    Ok((model, platt))
}

/// Split into a calibration set that preserves the natural class prior, and a
/// class-balanced training sample drawn from what is left.
///
/// Returns `(train_features, train_labels, holdout_features, holdout_labels)`.
///
/// The order matters, and it used to be the other way round: training took its
/// `n_positive`/`n_negative` from the whole set and Platt was fitted on the
/// leftovers. That makes the calibration set's class balance an artifact of the
/// sampling counts rather than a property of the data — and it moves when those
/// counts move. Measured on chr1/chr6/chr11, the insertion holdout was 45.9%
/// positive at 8000/20000 and 53.0% at 40000/60000, against a natural 41.7%;
/// deletion swung 39.0% -> 25.8% against a natural 35.4%. Platt calibrates the
/// mean predicted probability onto whatever prior it is shown, so a distorted
/// holdout yields a model whose probabilities are systematically off, by an
/// amount that differs per model. That is what stops one `--ml` threshold from
/// meaning the same thing for a CpG call and an indel call.
///
/// Carving the calibration set out first, uniformly at random, gives it the
/// natural prior by construction. Training still balances the classes from the
/// remainder, which is the right thing at CpG's 1% prevalence.
fn subsample_training_data(
    data: &TrainingData,
    n_positive: usize,
    n_negative: usize,
    seed: u64,
) -> Result<(Array2<f64>, Array1<f64>, Array2<f64>, Array1<f64>)> {
    const MAX_HOLDOUT: usize = 100_000;
    /// Share of the data reserved for calibration. Capped by `MAX_HOLDOUT`, and
    /// small enough that the training pool is never starved for the models whose
    /// candidate sets are small (the indel ones, ~100-150k examples).
    const HOLDOUT_FRACTION: f64 = 0.2;

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let mut all_indices: Vec<usize> = (0..data.len()).collect();
    all_indices.shuffle(&mut rng);

    let n_holdout = ((data.len() as f64 * HOLDOUT_FRACTION) as usize).min(MAX_HOLDOUT);
    let (holdout_slice, train_pool) = all_indices.split_at(n_holdout);
    let holdout: HashSet<usize> = holdout_slice.iter().copied().collect();

    // Balance the classes for training, from the pool the calibration set did
    // not claim.
    let (mut positive_indices, mut negative_indices) = (Vec::new(), Vec::new());
    for &i in train_pool {
        match data.labels.get(i) {
            Some(&label) if label == 1.0 => positive_indices.push(i),
            Some(_) => negative_indices.push(i),
            None => {}
        }
    }

    let n_pos_actual = positive_indices.len().min(n_positive);
    let n_neg_actual = negative_indices.len().min(n_negative);

    ensure!(n_pos_actual > 0, "No positive examples available for training");
    ensure!(n_neg_actual > 0, "No negative examples available for training");

    positive_indices.shuffle(&mut rng);
    negative_indices.shuffle(&mut rng);

    let mut train_indices: Vec<usize> = positive_indices
        .get(..n_pos_actual)
        .unwrap_or_default()
        .iter()
        .chain(negative_indices.get(..n_neg_actual).unwrap_or_default())
        .copied()
        .collect();

    let train_matrix = build_matrix(data, &mut train_indices)?;
    let holdout_matrix = build_matrix(data, &mut holdout.into_iter().collect::<Vec<_>>())?;

    Ok((train_matrix.0, train_matrix.1, holdout_matrix.0, holdout_matrix.1))
}

fn build_matrix(data: &TrainingData, indices: &mut [usize]) -> Result<(Array2<f64>, Array1<f64>)> {
    ensure!(
        !indices.is_empty(),
        "Cannot build matrix from empty indices — no holdout examples available. \
         This happens when all training examples are consumed for the training set, \
         leaving none for Platt calibration. Consider reducing --n-positive / --n-negative \
         or providing more training data."
    );

    indices.sort_unstable();

    let mut feature_rows = Vec::with_capacity(indices.len());
    let mut label_vec = Vec::with_capacity(indices.len());

    for &idx in indices.iter() {
        feature_rows.push(data.features[idx].row(0).to_owned());
        label_vec.push(data.labels[idx]);
    }

    let feature_views: Vec<_> = feature_rows.iter().map(|r| r.view()).collect();
    let features = ndarray::stack(Axis(0), &feature_views)
        .wrap_err_with(|| format!("Failed to stack feature arrays: {}", feature_rows.len()))?;

    let labels = Array1::from_vec(label_vec);

    Ok((features, labels))
}

/// How well a set of scores matches the outcomes they claim to predict.
///
/// Reported for the raw forest votes and again after Platt scaling, because
/// calibration is a step that can make things *worse* and nothing was checking.
/// Platt fits a two-parameter sigmoid, which is the right correction when the
/// distortion is sigmoid-shaped (SVMs, boosting) but not obviously right for a
/// random forest, whose vote fraction is already a bounded, roughly calibrated
/// probability.
#[derive(Debug, Clone, Copy)]
struct Calibration {
    /// Mean squared error against the outcome. Lower is better; sensitive to
    /// both calibration and discrimination.
    brier: f64,
    /// Expected calibration error: mean over bins of |mean predicted − observed
    /// frequency|, weighted by bin population. Purely a calibration measure.
    ece: f64,
    /// Mean prediction against the actual positive rate. A gap between these two
    /// is systematic over- or under-confidence.
    mean_score: f64,
    base_rate: f64,
}

impl fmt::Display for Calibration {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "brier={:.4} ece={:.4} mean={:.3} vs base={:.3}",
            self.brier, self.ece, self.mean_score, self.base_rate
        )
    }
}

fn calibration(scores: &[f64], labels: &[f64]) -> Option<Calibration> {
    const BINS: usize = 20;
    let n = scores.len();
    if n == 0 || n != labels.len() {
        return None;
    }

    let mut bin_sum = [0.0_f64; BINS];
    let mut bin_pos = [0.0_f64; BINS];
    let mut bin_n = [0.0_f64; BINS];
    let (mut brier, mut score_sum, mut pos) = (0.0, 0.0, 0.0);

    for (&s, &y) in scores.iter().zip(labels) {
        let s = s.clamp(0.0, 1.0);
        brier += (s - y) * (s - y);
        score_sum += s;
        pos += y;
        // The top edge belongs to the last bin rather than to a BINS+1'th one.
        let idx = ((s * BINS as f64) as usize).min(BINS - 1);
        bin_sum[idx] += s;
        bin_pos[idx] += y;
        bin_n[idx] += 1.0;
    }

    let ece = (0..BINS)
        .filter(|&i| bin_n[i] > 0.0)
        .map(|i| (bin_n[i] / n as f64) * (bin_sum[i] / bin_n[i] - bin_pos[i] / bin_n[i]).abs())
        .sum();

    Some(Calibration {
        brier: brier / n as f64,
        ece,
        mean_score: score_sum / n as f64,
        base_rate: pos / n as f64,
    })
}

/// Fit Platt scaling parameters A and B so that
/// `P(y=1|f) = 1 / (1 + exp(A*f + B))` is a well-calibrated probability.
///
/// Uses Newton's method with backtracking line search and Bayesian-smoothed
/// targets, following Lin, Lin, and Weng (2007).
fn fit_platt_scaling(scores: &[f64], labels: &[f64]) -> PlattScaling {
    let n = scores.len();
    if n == 0 {
        return PlattScaling::default();
    }

    let n_pos = labels.iter().filter(|&&y| y > 0.5).count() as f64;
    let n_neg = n as f64 - n_pos;

    if n_pos == 0.0 || n_neg == 0.0 {
        return PlattScaling::default();
    }

    // Bayesian-smoothed targets avoid log(0)
    let hi_target = (n_pos + 1.0) / (n_pos + 2.0);
    let lo_target = 1.0 / (n_neg + 2.0);

    let mut a = 0.0_f64;
    let mut b = ((n_neg + 1.0) / (n_pos + 1.0)).ln();

    // Initial objective value
    let mut fval = 0.0;
    for i in 0..n {
        let t = if labels[i] > 0.5 { hi_target } else { lo_target };
        let z = scores[i] * a + b;
        fval += if z >= 0.0 {
            t * z + (1.0 + (-z).exp()).ln()
        } else {
            (t - 1.0) * z + (1.0 + z.exp()).ln()
        };
    }

    const MAX_ITER: usize = 100;
    const MIN_STEP: f64 = 1e-10;
    const SIGMA: f64 = 1e-12;

    for _ in 0..MAX_ITER {
        let mut h11 = SIGMA;
        let mut h22 = SIGMA;
        let mut h12 = 0.0_f64;
        let mut g1 = 0.0_f64;
        let mut g2 = 0.0_f64;

        for i in 0..n {
            let t = if labels[i] > 0.5 { hi_target } else { lo_target };
            let z = scores[i] * a + b;
            let (p, q) = if z >= 0.0 {
                let ez = (-z).exp();
                (ez / (1.0 + ez), 1.0 / (1.0 + ez))
            } else {
                let ez = z.exp();
                (1.0 / (1.0 + ez), ez / (1.0 + ez))
            };
            let d2 = p * q;
            h11 += scores[i] * scores[i] * d2;
            h22 += d2;
            h12 += scores[i] * d2;
            let d1 = t - p;
            g1 += scores[i] * d1;
            g2 += d1;
        }

        // Newton step: H * [dA, dB] = -[g1, g2]
        let det = h11 * h22 - h12 * h12;
        let da = -(h22 * g1 - h12 * g2) / det;
        let db = -(-h12 * g1 + h11 * g2) / det;
        let gd = g1 * da + g2 * db;

        // Backtracking line search with Armijo condition
        let mut stepsize = 1.0_f64;
        while stepsize >= MIN_STEP {
            let new_a = a + stepsize * da;
            let new_b = b + stepsize * db;

            let mut newf = 0.0;
            for i in 0..n {
                let t = if labels[i] > 0.5 { hi_target } else { lo_target };
                let z = scores[i] * new_a + new_b;
                newf += if z >= 0.0 {
                    t * z + (1.0 + (-z).exp()).ln()
                } else {
                    (t - 1.0) * z + (1.0 + z.exp()).ln()
                };
            }

            if newf < fval + 0.0001 * stepsize * gd {
                a = new_a;
                b = new_b;
                fval = newf;
                break;
            }
            stepsize /= 2.0;
        }

        if stepsize < MIN_STEP || (g1.abs() < 1e-5 && g2.abs() < 1e-5) {
            break;
        }
    }

    PlattScaling { a, b }
}

/// Serialize a model to disk with LZ4 compression
fn serialize_model(model: &RastairFlatModel, path: ClioPath) -> Result<()> {
    let file = path.create().wrap_err("Failed to create output file for model serialization")?;
    let mut encoder =
        EncoderBuilder::new().level(16).build(file).wrap_err("Failed to create LZ4 encoder")?;

    rmp_serde::encode::write(&mut encoder, &model).wrap_err("Failed to serialize model")?;

    let (_output, result) = encoder.finish();
    result.wrap_err("Failed to finalize LZ4 compression")?;

    Ok(())
}

#[instrument(level = "info", skip_all, fields(model=%model_name))]
fn export_features_tsv(
    data: &TrainingData,
    model_name: &str,
    feature_names: &[&str],
    dir: &std::path::Path,
) -> Result<()> {
    if data.is_empty() {
        warn!("No examples to export — skipping TSV");
        return Ok(());
    }

    let path = dir.join(format!("{model_name}_features.tsv"));
    let file =
        File::create(&path).wrap_err_with(|| format!("Failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    write!(writer, "chrom\tpos\tlabel")?;
    for name in feature_names {
        write!(writer, "\t{name}")?;
    }
    writeln!(writer)?;

    for ((chrom, pos), (features, &label)) in
        data.positions.iter().zip(data.features.iter().zip(data.labels.iter()))
    {
        write!(writer, "{chrom}\t{pos}\t{label}")?;
        let row = features.row(0);
        for &v in row.iter() {
            write!(writer, "\t{v}")?;
        }
        writeln!(writer)?;
    }

    writer.flush().wrap_err("Failed to flush TSV writer")?;
    info!(examples = data.len(), path = %path.display(), "Exported features");
    Ok(())
}

fn export_feature_importances(model: &RandomForest, names: &[&str], path: &Path) -> Result<()> {
    let file = File::create(path).wrap_err_with(|| {
        format!("Failed to create feature importance file: {}", path.display())
    })?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "index\tfeature\timportance")
        .wrap_err("Failed to write feature importance header")?;
    for (idx, importance) in model.feature_importances().iter().enumerate() {
        // Fall back to the index if a name is missing, so a names/model length
        // mismatch is visible rather than silently truncating the output.
        let name = names.get(idx).copied().unwrap_or("<unknown>");
        writeln!(writer, "{idx}\t{name}\t{importance}")
            .wrap_err("Failed to write feature importance row")?;
    }

    info!(path = %path.display(), "Exported feature importances");

    Ok(())
}
