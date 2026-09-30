//! Train the five random forests a model file carries: CpG, de-novo CpG,
//! other SNVs, insertion and deletion.
//!
//! 1. Collect: run the column pipeline over every segment, label each SNV and
//!    indel candidate against the truth set, and file it under its model.
//! 2. Sample: draw `--n-positive`/`--n-negative` examples per model.
//! 3. Fit: a random forest per model, Platt-scaled on examples it did not see.
//! 4. Export: write the `RastairFlatModel`.

#[cfg(not(feature = "experimental-seqair"))]
use crate::call::process::calculate_pileup_metrics;
use crate::{
    call::{
        ml::DEFAULT_ML_THRESHOLD,
        pileup::indels::IndelAllele,
        process::{PileupMappingParams, get_pileups},
        variant_calling::indel_calling::{self, IndelParams, IndelPathway},
    },
    metrics::{
        MetricsForIndel, PileupMetrics,
        ml::{
            features::FeatureCalculator,
            types::{ByModel, MlFeatureSet, MlModel, PlattScaling, RastairFlatModel},
        },
    },
    regions::ConfidentRegions,
    sequence::{ChunkRegion, PileupReaders, ReaderParams, ReaderSource, SegmentationParams},
    utils::{cli, map_surrounding},
};
use biosphere::{FlatForest, MaxFeatures, RandomForest, RandomForestParameters};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, ContextCompat as _, Result, ensure, eyre};
use lz4::EncoderBuilder;
use ndarray::{Array1, Array2, Axis};
use rand::prelude::*;
use rayon::prelude::*;
use rust_htslib::bcf::{self, Read as _};
use seqair_types::{Base, Probability, RegionString, SmallVec, SmolStr};
use std::{
    collections::{BTreeMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    thread::available_parallelism,
};
use tracing::{debug, error, info, instrument, trace, warn};

#[derive(Debug, clap::Args)]
pub struct TrainModelParams {
    #[command(flatten)]
    reader: ReaderParams,

    /// Path to the ground truth file (VCF) to train with
    #[arg(help_heading = cli::sections::INPUT, value_hint=clap::ValueHint::FilePath)]
    truth: ClioPath,

    /// Restrict training to the intervals of a BED file (e.g. a GIAB
    /// high-confidence region file).
    ///
    /// Strongly recommended whenever the truth set has one. Outside its
    /// high-confidence regions a truth VCF asserts nothing, so a real variant
    /// there is simply missing from it and would be labelled *negative*. Those
    /// mislabelled negatives concentrate in the repetitive regions the truth set
    /// excludes — which is where indels live — so training without this teaches
    /// the model that indel-like signal in a repeat is false, using examples
    /// where the truth is merely unknown.
    ///
    /// Candidates outside the intervals are dropped, not relabelled.
    #[arg(long = "regions-file", short = 'R', help_heading = cli::sections::INPUT, value_hint = clap::ValueHint::FilePath)]
    regions_file: Option<ClioPath>,

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

/// The order the models' forest seeds are drawn from the run seed in.
const SEED_ORDER: [MlModel; MlModel::COUNT] =
    [MlModel::Cpg, MlModel::DenovoCpg, MlModel::Others, MlModel::Insertion, MlModel::Deletion];

/// Key for indexing positions in truth set
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct PositionKey {
    pos: u64,
    ref_base: Base,
    alt_base: Base,
}

/// Key for indexing indel positions in truth set
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct IndelKey {
    pos: u64,
    allele: IndelAllele,
}

/// The labels a candidate is checked against.
struct Truth {
    snps: HashSet<PositionKey>,
    indels: HashSet<IndelKey>,
    /// Outside these intervals the truth set makes no claim, so candidates
    /// there are dropped rather than labelled negative.
    confident: Option<ConfidentRegions>,
}

impl Truth {
    fn load(params: &TrainModelParams, regions: &[RegionString]) -> Result<Self> {
        let mut snps = HashSet::new();
        let mut indels = HashSet::new();
        for region in regions {
            let (region_snps, region_indels) =
                load_truth_vcf(&params.truth, region, params.threads)
                    .wrap_err_with(|| format!("Failed to load truth VCF for region {region}"))?;
            snps.extend(region_snps);
            indels.extend(region_indels);
        }

        let confident = params
            .regions_file
            .as_ref()
            .map(|p| ConfidentRegions::load(p.path()))
            .transpose()
            .wrap_err("Failed to load training regions")?;
        if confident.is_none() {
            warn!(
                "Training without --regions-file: candidates outside the truth set's \
                 high-confidence regions will be labelled negative even though the truth set \
                 makes no claim there, which is most severe for indels."
            );
        }

        Ok(Self { snps, indels, confident })
    }

    fn claims(&self, contig: &str, pos: u64) -> bool {
        self.confident.as_ref().is_none_or(|r| r.contains(contig, pos))
    }
}

/// Training data for a specific model type
#[derive(Default)]
struct TrainingData {
    features: Vec<Array2<f64>>,
    labels: Vec<f64>,
    positions: Vec<(SmolStr, u64)>,
}

impl TrainingData {
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

type Collected = ByModel<TrainingData>;

fn merge_collected(mut acc: Collected, other: Collected) -> Collected {
    for (model, data) in other {
        acc[model].merge(data);
    }
    acc
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

    create_output_dirs(params)?;

    let collected = collect(params)?;
    for (model, data) in collected.iter() {
        info!(
            model = model.name(),
            examples = data.len(),
            positives = data.positives(),
            "Collected training examples"
        );
    }

    if let Some(dir) = params.export_features.as_ref() {
        info!(dir = %dir.display(), "Exporting features as TSV");
        let names = params.ml_features.get_calculator().feature_names();
        for (model, data) in collected.iter() {
            export_features_tsv(data, model.name(), names.get(model), dir.path())?;
        }
    }

    let model = fit(params, seed, collected)?;
    serialize_model(&model, params.output.clone())
        .wrap_err_with(|| format!("Failed to serialize model to {}", params.output.display()))?;
    info!(path=%params.output, "Saved model");

    Ok(())
}

/// Every directory a run writes into, created before any work: failing here
/// costs seconds, failing after fitting costs the whole run.
fn create_output_dirs(params: &TrainModelParams) -> Result<()> {
    let model_dir = params.output.parent().wrap_err("output path invalid")?;
    std::fs::create_dir_all(model_dir).wrap_err_with(|| {
        format!("Failed to create output directory: {}", params.output.display())
    })?;
    for dir in [&params.feature_analytics, &params.export_features].into_iter().flatten() {
        std::fs::create_dir_all(dir.path())
            .wrap_err_with(|| format!("Failed to create directory: {}", dir.display()))?;
    }
    Ok(())
}

/// Label every candidate in the requested regions and file it under its model.
fn collect(params: &TrainModelParams) -> Result<Collected> {
    // Without explicit regions the truth set is read for chr12 only.
    let regions =
        params.reader.regions.as_ref().map(|input| input.regions().to_vec()).unwrap_or_else(|| {
            vec![RegionString { chromosome: "chr12".into(), start: None, end: None }]
        });
    let truth = Truth::load(params, &regions)?;

    let segmentation = SegmentationParams::default();
    let readers = params.reader.pileup_readers().wrap_err("Failed to read BAM/FASTA files")?;
    let segments: Vec<ChunkRegion> = readers
        .segments(segmentation.segment_max_length, segmentation.segment_overlap)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();
    ensure!(!segments.is_empty(), "No segments found in BAM file");
    info!("Processing {} segments to collect training data", segments.len());

    let collector = SegmentCollector {
        truth: &truth,
        calculator: &*params.ml_features.get_calculator(),
        // Training examples come from the ML pathway, not the hard-filter chain.
        indel_params: IndelParams {
            experimental_indels: Some(IndelPathway::Ml),
            ..IndelParams::default()
        },
    };
    let readers = ReaderSource::from(readers);
    let collected = rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("training-worker-{idx}"))
        .num_threads(params.threads)
        .start_handler(|idx| trace!(idx, "Starting training worker thread"))
        .exit_handler(|idx| trace!(idx, "Closing training worker thread"))
        .build()
        .wrap_err("Failed to create thread pool for rayon")?
        .install(|| {
            segments
                .par_iter()
                .map_init(
                    || readers.fork(),
                    |readers, segment| {
                        let _span =
                            tracing::info_span!("collect_segment", region = %segment.region)
                                .entered();
                        let readers = readers
                            .as_mut()
                            .map_err(|e| eyre!("Failed to open readers in worker thread: {e:#}"))?;
                        collector.segment(readers, segment)
                    },
                )
                .filter_map(|result| {
                    result
                        .inspect_err(|e| {
                            warn!(
                                error = format!("{e:#}"),
                                "Failed to collect training data from segment"
                            );
                        })
                        .ok()
                })
                .collect::<Vec<_>>()
        });

    Ok(collected.into_iter().fold(ByModel::from_fn(|_| TrainingData::default()), merge_collected))
}

/// The columns of one segment, as the feature extractors read them.
fn segment_columns(
    readers: &mut PileupReaders,
    segment: &ChunkRegion,
) -> Result<Vec<PileupMetrics>> {
    let mapping = PileupMappingParams { call_indels: true, ..Default::default() };

    #[cfg(feature = "experimental-seqair")]
    let columns = {
        let (_segment, columns) =
            get_pileups(readers, segment, &mapping).wrap_err("Failed to build pileups")?;
        columns.collect()
    };
    #[cfg(not(feature = "experimental-seqair"))]
    let columns = {
        let (segment, pileups) =
            get_pileups(readers, segment, &mapping).wrap_err("Failed to build pileups")?;
        calculate_pileup_metrics(pileups, &segment)
            .filter_map(|metrics| {
                metrics
                    .inspect_err(|e| {
                        warn!(error = format!("{e:#}"), "Failed to calculate pileup metrics");
                    })
                    .ok()
            })
            .collect()
    };

    Ok(columns)
}

/// Turns one segment's candidates into labelled examples.
struct SegmentCollector<'a> {
    truth: &'a Truth,
    calculator: &'a dyn FeatureCalculator,
    indel_params: IndelParams,
}

impl SegmentCollector<'_> {
    fn segment(&self, readers: &mut PileupReaders, segment: &ChunkRegion) -> Result<Collected> {
        let mut columns = segment_columns(readers, segment)?;
        let mut collected = ByModel::from_fn(|_| TrainingData::default());
        map_surrounding(
            &mut columns,
            |before, column, after| {
                let pos = u64::from(column.pos);
                if !self.truth.claims(&segment.contig, pos) {
                    return Ok(());
                }
                self.snvs(&mut collected, &segment.contig, before, column, after);
                self.indels(&mut collected, &segment.contig, column);
                Ok(())
            },
            "failed to extract training features, skipping",
        );
        Ok(collected)
    }

    fn snvs(
        &self,
        collected: &mut Collected,
        contig: &SmolStr,
        before: Option<&PileupMetrics>,
        column: &PileupMetrics,
        after: Option<&PileupMetrics>,
    ) {
        let pos = u64::from(column.pos);
        let ref_base = column.reference_base;
        for alt in &column.alts {
            let alt_base = alt.base;
            if ref_base == Base::Unknown || alt_base == Base::Unknown {
                continue;
            }
            let Some(candidate) = column.alt_metrics(alt_base) else { continue };
            let label = self.truth.snps.contains(&PositionKey { pos, ref_base, alt_base });

            let (model, features) = if candidate.is_evidence_for_methylation() {
                (MlModel::Cpg, self.calculator.calculate_cpg(&candidate, before, after))
            } else if *alt.metrics.denovo {
                (
                    MlModel::DenovoCpg,
                    self.calculator.calculate_denovo_cpg(&candidate, before, after),
                )
            } else {
                (MlModel::Others, self.calculator.calculate_others(&candidate, before, after))
            };
            add_candidate(&mut collected[model], model, features, label, contig, pos);
        }
    }

    fn indels(&self, collected: &mut Collected, contig: &SmolStr, column: &PileupMetrics) {
        let Some(indel_data) = column.indel_data.as_ref() else { return };
        let pos = u64::from(column.pos);
        let tract = u32::from(indel_data.homopolymer_run.max(indel_data.dinucleotide_run));
        let calls =
            indel_calling::call_indels(&indel_data.counts, &self.indel_params, true, tract, false);
        for call in &calls {
            let label = self.truth.indels.contains(&IndelKey { pos, allele: call.allele.clone() });
            let candidate = MetricsForIndel { metrics: column, indel: call };
            let (model, features) = match call.allele {
                IndelAllele::Insertion(_) => {
                    (MlModel::Insertion, self.calculator.calculate_insertion(&candidate))
                }
                IndelAllele::Deletion(_) => {
                    (MlModel::Deletion, self.calculator.calculate_deletion(&candidate))
                }
            };
            add_candidate(&mut collected[model], model, features, label, contig, pos);
        }
    }
}

/// Keep a candidate whose features could be computed and are all finite.
fn add_candidate(
    data: &mut TrainingData,
    model: MlModel,
    features: Result<Array2<f32>>,
    in_truth: bool,
    contig: &SmolStr,
    pos: u64,
) {
    match features {
        Ok(features) if !features.is_any_nan() => {
            let label = if in_truth { 1.0 } else { 0.0 };
            data.add_example(features, label, contig.clone(), pos);
        }
        Ok(_) => {}
        Err(error) => {
            debug!(
                model = model.name(),
                error = format!("{error:#}"),
                "No features for training candidate"
            );
        }
    }
}

/// Load truth VCF and create an index of variant positions (SNPs and indels).
#[instrument(level = "info", skip_all)]
fn load_truth_vcf(
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

/// The SNP and indel keys of one PASS truth record.
fn process_truth_record(record: &bcf::Record) -> (SmallVec<PositionKey, 2>, SmallVec<IndelKey, 2>) {
    let mut keys = (SmallVec::new(), SmallVec::new());
    if !record.has_filter("PASS".as_bytes()) {
        return keys;
    }
    let Ok(pos) = u64::try_from(record.pos()) else { return keys };
    let alleles = record.alleles();
    let Some((&ref_allele, alts)) = alleles.split_first() else { return keys };
    for &alt_allele in alts {
        match TruthKey::of(pos, ref_allele, alt_allele) {
            Some(TruthKey::Snv(key)) => keys.0.push(key),
            Some(TruthKey::Indel(key)) => keys.1.push(key),
            None => {}
        }
    }
    keys
}

/// What one ALT allele of a truth record labels.
#[derive(Debug, PartialEq, Eq)]
enum TruthKey {
    Snv(PositionKey),
    Indel(IndelKey),
}

impl TruthKey {
    /// The allele's key once the suffix it shares with REF past the anchor
    /// base is stripped, which is how a multi-allelic record pads its shorter
    /// alleles.
    ///
    /// A single-base substitution or an insertion or deletion after the shared
    /// anchor base has one; MNPs, complex and symbolic alleles do not.
    fn of(pos: u64, ref_allele: &[u8], alt_allele: &[u8]) -> Option<Self> {
        let (&ref_anchor, ref_rest) = ref_allele.split_first()?;
        let (&alt_anchor, alt_rest) = alt_allele.split_first()?;
        let (ref_rest, alt_rest) = strip_common_suffix(ref_rest, alt_rest);
        let known = |bases: &[u8]| -> Option<SmallVec<Base, 4>> {
            bases.iter().map(|&b| Some(Base::from(b)).filter(|&b| b != Base::Unknown)).collect()
        };
        match (ref_rest.is_empty(), alt_rest.is_empty()) {
            (true, true) if ref_anchor != alt_anchor => {
                let (ref_base, alt_base) = (Base::from(ref_anchor), Base::from(alt_anchor));
                (ref_base != Base::Unknown && alt_base != Base::Unknown)
                    .then_some(Self::Snv(PositionKey { pos, ref_base, alt_base }))
            }
            (false, true) if ref_anchor == alt_anchor => known(ref_rest)
                .map(|bases| Self::Indel(IndelKey { pos, allele: IndelAllele::Deletion(bases) })),
            (true, false) if ref_anchor == alt_anchor => known(alt_rest)
                .map(|bases| Self::Indel(IndelKey { pos, allele: IndelAllele::Insertion(bases) })),
            _ => None,
        }
    }
}

/// Both alleles without the longest suffix they share.
fn strip_common_suffix<'a>(mut left: &'a [u8], mut right: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    while let (Some((l, left_rest)), Some((r, right_rest))) =
        (left.split_last(), right.split_last())
        && l == r
    {
        (left, right) = (left_rest, right_rest);
    }
    (left, right)
}

/// Fit every model's forest, one after the other: biosphere parallelises
/// each fit over `--threads` in a pool of its own.
fn fit(params: &TrainModelParams, seed: u64, collected: Collected) -> Result<RastairFlatModel> {
    let seeds = forest_seeds(seed);
    let calculator = params.ml_features.get_calculator();
    let (counts, names) = (calculator.feature_num(), calculator.feature_names());

    let trained = collected.try_map(|model, data| -> Result<_> {
        let (forest, platt) = fit_model(model, &data, params, seeds[model])
            .wrap_err_with(|| format!("Failed to train the {} model", model.name()))?;
        if let Some(dir) = params.feature_analytics.as_ref() {
            let path = dir.path().join(format!("{}_feature_importances.csv", model.name()));
            export_feature_importances(&forest, names.get(model), &path).wrap_err_with(|| {
                format!("Failed to export {} feature importances", model.name())
            })?;
        }
        Ok((FlatForest::from_forest(&forest, counts.get(model)), platt))
    })?;

    #[derive(Debug, serde::Serialize)]
    struct ModelReport<'a> {
        forest: &'a biosphere::ForestMeta,
        scaling: PlattScaling,
    }
    let report: BTreeMap<_, _> = trained
        .iter()
        .map(|(model, (forest, scaling))| {
            (model.name(), ModelReport { forest: &forest.meta, scaling: *scaling })
        })
        .collect();
    info!(report = ?report, "Trained models");

    Ok(RastairFlatModel::from_trained(params.ml_features, trained))
}

/// One forest seed per model, drawn from the run seed in [`SEED_ORDER`].
fn forest_seeds(seed: u64) -> ByModel<u64> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut seeds = ByModel::from_fn(|_| 0);
    for model in SEED_ORDER {
        seeds[model] = rng.random();
    }
    seeds
}

#[instrument(level = "info", skip_all, fields(model = model.name()))]
fn fit_model(
    model: MlModel,
    data: &TrainingData,
    params: &TrainModelParams,
    seed: u64,
) -> Result<(RandomForest, PlattScaling)> {
    ensure!(
        !data.is_empty(),
        "No training data collected for the {} model, so no model file can be written",
        model.name()
    );

    info!(seed, examples = data.len(), "Training model");

    // Subsample for training, keep held-out data for Platt calibration
    let (train_features, train_labels, holdout_features, holdout_labels) = subsample_training_data(
        data,
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

    // `max_depth == 0` means unbounded (grow until pure).
    let max_depth = (params.model_params.max_depth != 0).then_some(params.model_params.max_depth);
    let rf_params = RandomForestParameters::default()
        .with_max_features(MaxFeatures::Value(params.model_params.max_features))
        .with_n_estimators(params.model_params.n_trees)
        .with_max_depth(max_depth)
        .with_min_samples_leaf(params.model_params.min_samples_leaf)
        .with_n_jobs(i32::try_from(params.threads).ok())
        .with_seed(seed);

    let mut forest = RandomForest::new(rf_params);
    forest.fit(&train_features.view(), &train_labels.view());

    let raw_scores = forest.predict(&holdout_features.view());
    let platt = fit_platt_scaling(
        raw_scores.as_slice().wrap_err("Holdout scores are not contiguous")?,
        holdout_labels.as_slice().wrap_err("Holdout labels are not contiguous")?,
    );

    if platt.a == 1.0 && platt.b == 0.0 {
        error!(
            "Platt scaling is identity (a=1.0, b=0.0) — model likely failed to learn. \
             Check class balance, feature quality, and consider reducing n-negative."
        );
    }

    Ok((forest, platt))
}

/// Subsample training data to balance positive and negative examples.
///
/// Returns `(train_features, train_labels, holdout_features, holdout_labels)`.
/// The held-out set is capped at `MAX_HOLDOUT` to keep memory bounded while
/// still providing enough data for a stable Platt scaling fit.
fn subsample_training_data(
    data: &TrainingData,
    n_positive: usize,
    n_negative: usize,
    seed: u64,
) -> Result<(Array2<f64>, Array1<f64>, Array2<f64>, Array1<f64>)> {
    const MAX_HOLDOUT: usize = 100_000;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

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

    // Sample indices for training
    let n_pos_actual = positive_indices.len().min(n_positive);
    let n_neg_actual = negative_indices.len().min(n_negative);

    ensure!(n_pos_actual > 0, "No positive examples available for training");
    ensure!(n_neg_actual > 0, "No negative examples available for training");

    positive_indices.shuffle(&mut rng);
    negative_indices.shuffle(&mut rng);

    let selected_pos = &positive_indices[..n_pos_actual];
    let selected_neg = &negative_indices[..n_neg_actual];

    let train_indices: HashSet<usize> =
        selected_pos.iter().chain(selected_neg.iter()).copied().collect();

    // Build training matrix
    let train_matrix = build_matrix(data, &mut train_indices.iter().copied().collect::<Vec<_>>())?;

    // Build held-out matrix from remaining indices, capped for memory
    let mut holdout_indices: Vec<usize> =
        (0..data.len()).filter(|i| !train_indices.contains(i)).collect();
    holdout_indices.shuffle(&mut rng);
    holdout_indices.truncate(MAX_HOLDOUT);

    let holdout_matrix = build_matrix(data, &mut holdout_indices)?;

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

/// Fit Platt scaling parameters A and B so that
/// `P(y=1|f) = 1 / (1 + exp(A*f + B))` is a well-calibrated probability.
///
/// Uses Newton's method with backtracking line search and Bayesian-smoothed
/// targets, following Lin, Lin, and Weng (2007).
fn fit_platt_scaling(scores: &[f64], labels: &[f64]) -> PlattScaling {
    const MAX_ITER: usize = 100;
    const MIN_STEP: f64 = 1e-10;
    const SIGMA: f64 = 1e-12;

    let samples = || scores.iter().copied().zip(labels.iter().map(|&y| y > 0.5));
    let (n_pos, n) = samples().fold((0.0_f64, 0.0_f64), |(pos, all), (_, y)| {
        (if y { pos + 1.0 } else { pos }, all + 1.0)
    });
    let n_neg = n - n_pos;
    if n_pos == 0.0 || n_neg == 0.0 {
        return PlattScaling::default();
    }

    // Bayesian-smoothed targets avoid log(0)
    let hi_target = (n_pos + 1.0) / (n_pos + 2.0);
    let lo_target = 1.0 / (n_neg + 2.0);
    let target = |positive: bool| if positive { hi_target } else { lo_target };
    let objective = |a: f64, b: f64| {
        samples().fold(0.0, |sum, (score, y)| {
            let t = target(y);
            let z = score * a + b;
            sum + if z >= 0.0 {
                t * z + (1.0 + (-z).exp()).ln()
            } else {
                (t - 1.0) * z + (1.0 + z.exp()).ln()
            }
        })
    };

    let mut a = 0.0_f64;
    let mut b = ((n_neg + 1.0) / (n_pos + 1.0)).ln();
    let mut fval = objective(a, b);

    for _ in 0..MAX_ITER {
        let mut h11 = SIGMA;
        let mut h22 = SIGMA;
        let mut h12 = 0.0_f64;
        let mut g1 = 0.0_f64;
        let mut g2 = 0.0_f64;

        for (score, y) in samples() {
            let z = score * a + b;
            let (p, q) = if z >= 0.0 {
                let ez = (-z).exp();
                (ez / (1.0 + ez), 1.0 / (1.0 + ez))
            } else {
                let ez = z.exp();
                (1.0 / (1.0 + ez), ez / (1.0 + ez))
            };
            let d2 = p * q;
            h11 += score * score * d2;
            h22 += d2;
            h12 += score * d2;
            let d1 = target(y) - p;
            g1 += score * d1;
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
            let (new_a, new_b) = (a + stepsize * da, b + stepsize * db);
            let newf = objective(new_a, new_b);
            if newf < fval + 0.0001 * stepsize * gd {
                (a, b, fval) = (new_a, new_b, newf);
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
    dir: &Path,
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A retrain reproduces a model only if every forest gets the seed it got
    /// before, so the draw order is pinned.
    #[test]
    fn forest_seeds_are_drawn_in_training_order() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let expected: Vec<u64> = (0..MlModel::COUNT).map(|_| rng.random()).collect();
        let seeds = forest_seeds(7);
        let order = [
            MlModel::Cpg,
            MlModel::DenovoCpg,
            MlModel::Others,
            MlModel::Insertion,
            MlModel::Deletion,
        ];
        assert_eq!(order.map(|model| seeds[model]).to_vec(), expected);
    }

    fn truth_key(ref_allele: &str, alt_allele: &str) -> Option<TruthKey> {
        TruthKey::of(10, ref_allele.as_bytes(), alt_allele.as_bytes())
    }

    fn deletion(bases: &[Base]) -> Option<TruthKey> {
        Some(TruthKey::Indel(IndelKey { pos: 10, allele: IndelAllele::Deletion(bases.into()) }))
    }

    /// `REF=ATT ALT=A,AT` deletes two bases and one; the second allele is
    /// padded with the T it shares with REF.
    #[test]
    fn a_padded_allele_of_a_multi_allelic_record_keys_its_own_indel() {
        assert_eq!(truth_key("ATT", "A"), deletion(&[Base::T, Base::T]));
        assert_eq!(truth_key("ATT", "AT"), deletion(&[Base::T]));
        assert_eq!(truth_key("AA", "A"), deletion(&[Base::A]));
        assert_eq!(
            truth_key("A", "ATG"),
            Some(TruthKey::Indel(IndelKey {
                pos: 10,
                allele: IndelAllele::Insertion([Base::T, Base::G].as_slice().into()),
            }))
        );
    }

    /// `REF=AT ALT=GT` is an A>G substitution padded to the record's REF.
    #[test]
    fn a_padded_substitution_keys_an_snv() {
        let snv =
            |ref_base, alt_base| Some(TruthKey::Snv(PositionKey { pos: 10, ref_base, alt_base }));
        assert_eq!(truth_key("AT", "GT"), snv(Base::A, Base::G));
        assert_eq!(truth_key("C", "T"), snv(Base::C, Base::T));
    }

    #[test]
    fn complex_and_symbolic_alleles_key_nothing() {
        for (ref_allele, alt_allele) in [
            ("ATG", "AC"),
            ("AT", "GC"),
            ("AT", "AT"),
            ("A", "<DEL>"),
            ("A", "*"),
            ("TA", "A"),
            ("A", "N"),
        ] {
            assert_eq!(truth_key(ref_allele, alt_allele), None, "{ref_allele}>{alt_allele}");
        }
    }
}
