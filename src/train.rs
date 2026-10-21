//! Train the three random forests a model file carries: CpG, de-novo CpG and
//! other SNVs.
//!
//! 1. Collect: run the column pipeline over every segment, label each SNV
//!    candidate against the truth set, and file it under its model.
//!    Each model keeps a bounded uniform sample of what it is offered.
//! 2. Sample: draw `--n-positive`/`--n-negative` examples per model.
//! 3. Fit: a random forest per model, Platt-scaled on examples it did not see.
//! 4. Export: write the `RastairFlatModel`.
//!
//! ```bash
//! rastair ml train -r hg38.fa.gz taps.bam HG001_benchmark.vcf.gz \
//!   -R HG001_benchmark.bed -l "chr1 chr6 chr11" --seed <n> -o models/mymodel.rff.mpk.lz4 -@ 8
//! ```
//!
//! - Positionals are BAM first, then the truth VCF. The defaults are the
//!   recipe of the bundled model (`models/rastair_with_indels.rff.mpk.lz4`),
//!   so a bare command line retrains it.
//! - The truth VCF needs a `.csi` (a `.tbi` alone is refused), and `-R` is a
//!   plain-text BED. Without `-R`, candidates outside the truth set's
//!   confident regions are labelled negative.
//! - Train on GIAB HG001 and evaluate against Platinum Genomes on a held-out
//!   chromosome. All five forests must train, or no model is written.
//! - The sample is a keyed reservoir per model and class, seeded per segment
//!   and model from `--seed` by `KeySource`, so it depends neither on
//!   scheduling nor on other models' pools. The draw takes at most four fifths
//!   of a class; the holdout keeps the population's class ratio, or Platt
//!   scaling would shift what an ML threshold means.
//! - The `Collected training examples` log line reports pool, kept, requested
//!   and drawn per model. `--export-features` writes a 1-based `pos`.

mod reservoir;

use crate::{
    call::{
        SegmentColumns, ml::MachineLearningParams, process::PileupMappingParams, segment_columns,
    },
    metrics::{
        PileupMetrics,
        ml::{
            features::{FeatureCalculator, FeatureNum},
            types::{
                ByModel, MachineLearning, MlFeatureSet, MlModel, PlattScaling, RastairFlatModel,
            },
        },
    },
    regions::ConfidentRegions,
    sequence::{ChunkRegion, ReaderParams, Region, SegmentationParams, SeqairReaders},
    utils::{cli, map_surrounding},
    verify::for_each_record_in_regions,
};
use biosphere::{FlatForest, MaxFeatures, RandomForest, RandomForestParameters};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, ContextCompat as _, Result, ensure, eyre};
use lz4::EncoderBuilder;
use ndarray::Array2;
use rand::prelude::*;
use rayon::prelude::*;
use reservoir::{ByLabel, KeySource, Label, SamplingRequest, TrainingData};
use rust_htslib::bcf::{self, Read as _};
use seqair_types::{Base, Pos0, RegionString, SmolStr};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    thread::available_parallelism,
};
use tracing::{debug, info, instrument, trace, warn};

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

    /// Model file to write
    #[arg(short = 'o', long = "output", default_value = "models/rastair.rff.mpk.lz4")]
    #[arg(help_heading = cli::sections::OUTPUT, value_hint=clap::ValueHint::FilePath)]
    output: ClioPath,

    /// Export collected features and labels as TSV files to this directory.
    /// One file per model type: `cpg_features.tsv`, `denovo_features.tsv`,
    /// `other_features.tsv`.
    #[arg(long, help_heading = cli::sections::OUTPUT, value_hint=clap::ValueHint::DirPath)]
    export_features: Option<ClioPath>,

    /// Export features importances as TSV files to this directory.
    #[arg(long, help_heading = cli::sections::OUTPUT, value_hint=clap::ValueHint::DirPath)]
    feature_analytics: Option<ClioPath>,

    #[command(flatten)]
    model_params: ModelParameters,

    #[arg(long, default_value_t = MlFeatureSet::Standard)]
    ml_features: MlFeatureSet,

    /// Also fit the forests `call --snv-engine` re-scores SNV candidates
    /// with; the model file then carries both sets.
    ///
    /// The candidates are collected a second time the way the cascade hands
    /// them on: the forests fitted first score every candidate, the pair-HMM
    /// engine scores the ones they keep on the clusters
    /// `call --experimental-indels` builds, and only those are fitted on.
    #[arg(long, help_heading = cli::sections::TRAINING)]
    snv_engine: bool,

    /// Number of threads to use
    #[arg(short='@', long = "threads", env = "RASTAIR_THREADS", default_value_t = available_parallelism().map(|n|n.get()).unwrap_or(2).max(1))]
    #[arg(help_heading = cli::sections::PROCESSING)]
    threads: usize,
}

#[derive(Debug, clap::Args)]
struct ModelParameters {
    /// Number of trees in the random forest
    #[arg(long = "n-trees", default_value_t = 800)]
    #[arg(help_heading = cli::sections::TRAINING)]
    n_trees: usize,

    /// Number of features to consider at each split (mtry parameter)
    #[arg(long = "max-features", default_value_t = 4)]
    #[arg(help_heading = cli::sections::TRAINING)]
    max_features: usize,

    /// Maximum tree depth. `0` grows trees until pure (unbounded), which can
    /// produce very large models on noisy/poorly-separable data such as indels.
    /// See <https://scikit-learn.org/stable/modules/tree.html>.
    #[arg(long = "max-depth", default_value_t = 40)]
    #[arg(help_heading = cli::sections::TRAINING)]
    max_depth: usize,

    /// Minimum number of samples required at each leaf. Larger values prevent
    /// the forest from memorising individual samples, regularising noisy data
    /// and shrinking the model; use 1 for fully-grown leaves.
    #[arg(long = "min-samples-leaf", default_value_t = 10)]
    #[arg(help_heading = cli::sections::TRAINING)]
    min_samples_leaf: usize,

    /// Number of positive examples to draw for each SNV model (CpG, de-novo
    /// CpG, other).
    #[arg(long = "n-positive", default_value_t = 8_000)]
    #[arg(help_heading = cli::sections::TRAINING)]
    n_positive: usize,

    /// Number of negative examples to draw for each SNV model.
    ///
    /// Almost every SNV candidate is a sequencing or conversion artefact, and
    /// the draw should be about as negative-heavy as the candidates a call
    /// meets.
    #[arg(long = "n-negative", default_value_t = 200_000)]
    #[arg(help_heading = cli::sections::TRAINING)]
    n_negative: usize,

    /// Random seed for reproducibility (subsampling and forest training).
    /// Omit for a random seed.
    #[arg(long)]
    #[arg(help_heading = cli::sections::TRAINING)]
    seed: Option<u64>,
}

impl ModelParameters {
    /// The draw every model asks for.
    fn sampling_request(&self) -> SamplingRequest {
        ByLabel { positive: self.n_positive, negative: self.n_negative }
    }
}

/// The order the models' forest seeds are drawn from the run seed in.
const SEED_ORDER: [MlModel; MlModel::COUNT] = [MlModel::Cpg, MlModel::DenovoCpg, MlModel::Others];

/// Key for indexing positions in truth set
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
struct PositionKey {
    pos: u64,
    ref_base: Base,
    alt_base: Base,
}

/// The labels a candidate is checked against.
struct Truth {
    /// Only the contigs the truth VCF declares: it makes no claim about any
    /// other, so candidates there are not trained on.
    contigs: HashMap<SmolStr, ContigTruth>,
    /// Outside these intervals the truth set makes no claim, so candidates
    /// there are dropped rather than labelled negative.
    confident: Option<ConfidentRegions>,
}

/// The true variants of one contig.
#[derive(Debug, Default)]
struct ContigTruth {
    snps: HashSet<PositionKey>,
}

impl Truth {
    fn load(params: &TrainModelParams, regions: &[RegionString]) -> Result<Self> {
        let contigs = truth_keys(params.truth.path(), regions, params.threads)
            .wrap_err_with(|| format!("Failed to load truth VCF {}", params.truth.display()))?;

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
                 makes no claim there."
            );
        }

        Ok(Self { contigs, confident })
    }

    fn claims(&self, contig: &str, pos: u64) -> bool {
        self.confident.as_ref().is_none_or(|r| r.contains(contig, pos))
    }
}

type Collected = ByModel<TrainingData>;

fn merge_collected(mut acc: Collected, other: Collected) -> Result<Collected> {
    for (model, data) in other {
        acc[model].merge(data)?;
    }
    Ok(acc)
}

#[instrument(level = "debug", skip_all)]
pub fn train_model(params: &TrainModelParams) -> Result<()> {
    let seed = params.model_params.seed.unwrap_or_else(rand::random);
    let request = params.model_params.sampling_request();
    info!(
        seed,
        n_trees = params.model_params.n_trees,
        max_features = params.model_params.max_features,
        max_depth = params.model_params.max_depth,
        min_samples_leaf = params.model_params.min_samples_leaf,
        n_positive = request.positive,
        n_negative = request.negative,
        "Training parameters",
    );

    create_output_dirs(params)?;

    let collected = collect(params, request, seed, None)?;

    if let Some(dir) = params.export_features.as_ref() {
        info!(dir = %dir.display(), "Exporting features as TSV");
        let names = params.ml_features.get_calculator().feature_names();
        for (model, data) in collected.iter() {
            export_features_tsv(data, model.name(), names.get(model), dir.path())?;
        }
    }

    let mut model = fit(params, request, seed, &collected, "")?;
    if params.snv_engine {
        // The engine forests are fitted on what the cascade hands them: the
        // alts these forests keep, with the engine's features from their
        // clusters. Fitted on every gated alt instead, they would learn
        // `ev_scored` 0 as "the gate refused".
        let stage1 = MachineLearningParams::default().with_forests(model.clone(), 1)?;
        let collected = collect(params, request, seed, Some(&stage1))?;
        model.snv_engine = Some(Box::new(fit(params, request, seed, &collected, "snv_engine_")?));
    }
    serialize_model(&model, params.output.clone())
        .wrap_err_with(|| format!("Failed to serialize model to {}", params.output.display()))?;
    info!(path=%params.output, "Saved model");

    Ok(())
}

/// Every directory a run writes into, created before any work: failing here
/// costs seconds, failing after fitting costs the whole run.
fn create_output_dirs(params: &TrainModelParams) -> Result<()> {
    ensure!(
        !params.output.is_dir(),
        "--output {} is a directory; name the model file to write, e.g. {}",
        params.output.display(),
        params.output.path().join("rastair.rff.mpk.lz4").display(),
    );
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
/// With `stage1`, only the candidates the cascade re-scores after those
/// forests (the engine forests' rows).
fn collect(
    params: &TrainModelParams,
    request: SamplingRequest,
    seed: u64,
    stage1: Option<&MachineLearning>,
) -> Result<Collected> {
    let segmentation = SegmentationParams::default();
    let readers = params.reader.pileup_readers().wrap_err("Failed to read BAM/FASTA files")?;
    let segments: Vec<ChunkRegion> = readers
        .segments(segmentation.segment_max_length, segmentation.segment_overlap)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();
    ensure!(!segments.is_empty(), "No segments found in BAM file");
    info!("Processing {} segments to collect training data", segments.len());
    let truth = Truth::load(params, &covered_regions(&segments)?)?;

    let calculator = params.ml_features.get_calculator();
    let collector = SegmentCollector {
        truth: &truth,
        stage1,
        calculator: &*calculator,
        feature_num: calculator.feature_num(),
        request,
    };
    // Reduced rather than collected: holding every segment's sample until a
    // final merge would keep `segments x cap` examples alive.
    let (mut collected, failed_segments) = rayon::ThreadPoolBuilder::new()
        .thread_name(|idx| format!("training-worker-{idx}"))
        .num_threads(params.threads)
        .start_handler(|idx| trace!(idx, "Starting training worker thread"))
        .exit_handler(|idx| trace!(idx, "Closing training worker thread"))
        .build()
        .wrap_err("Failed to create thread pool for rayon")?
        .install(|| {
            segments
                .par_iter()
                .enumerate()
                .map_init(
                    || readers.fork(),
                    |readers, (index, segment)| {
                        let _span =
                            tracing::info_span!("collect_segment", region = %segment.region)
                                .entered();
                        let readers = readers
                            .as_mut()
                            .map_err(|e| eyre!("Failed to open readers in worker thread: {e:#}"))?;
                        collector.segment(
                            readers,
                            segment,
                            &mut ByModel::from_fn(|model| {
                                KeySource::for_segment(seed, index, model)
                            }),
                        )
                    },
                )
                .map(|result| -> Result<(Collected, usize)> {
                    Ok(result.map_or_else(
                        |e| {
                            warn!(
                                error = format!("{e:#}"),
                                "Failed to collect training data from segment"
                            );
                            (collector.empty(), 1)
                        },
                        |collected| (collected, 0),
                    ))
                })
                .try_reduce(
                    || (collector.empty(), 0),
                    |(left, left_failed), (right, right_failed)| {
                        Ok((merge_collected(left, right)?, left_failed + right_failed))
                    },
                )
        })?;

    // The reduce merged in completion order; finishing makes the sample a
    // function of the seed alone.
    for model in MlModel::ALL {
        collected[model].finish();
    }
    for (model, data) in collected.iter() {
        let (seen, kept) = (data.seen(), data.kept());
        let draw = request.draw_from(kept);
        info!(
            model = model.name(),
            forests = if stage1.is_some() { "snv_engine_" } else { "" },
            pool_positive = seen.positive,
            pool_negative = seen.negative,
            kept_positive = kept.positive,
            kept_negative = kept.negative,
            requested_positive = request.positive,
            requested_negative = request.negative,
            draw_positive = draw.positive,
            draw_negative = draw.negative,
            rejected = data.rejected(),
            "Collected training examples"
        );
        ensure!(
            kept.positive >= 2 && kept.negative >= 2,
            "The {} model collected {} positive and {} negative examples ({failed_segments} \
             segments failed); it needs at least two of each, one to fit on and one to \
             calibrate with",
            model.name(),
            kept.positive,
            kept.negative,
        );
    }
    Ok(collected)
}

/// Turns one segment's candidates into labelled examples.
struct SegmentCollector<'a> {
    truth: &'a Truth,
    /// With the forests of the cascade's first stage, the candidates are the
    /// ones the SNV engine scores after them, as in `call --snv-engine`.
    stage1: Option<&'a MachineLearning>,
    calculator: &'a dyn FeatureCalculator,
    feature_num: FeatureNum,
    request: SamplingRequest,
}

impl SegmentCollector<'_> {
    fn empty(&self) -> Collected {
        ByModel::from_fn(|model| {
            TrainingData::for_request(self.feature_num.get(model), self.request)
        })
    }

    /// Candidates are offered in column order, which together with `keys`
    /// fixes which of them the reservoirs keep.
    fn segment(
        &self,
        readers: &mut SeqairReaders,
        segment: &ChunkRegion,
        keys: &mut ByModel<KeySource>,
    ) -> Result<Collected> {
        let Some(truth) = self.truth.contigs.get(&segment.contig) else {
            return Ok(self.empty());
        };
        // With the first stage, the columns `call --snv-engine
        // --experimental-indels` builds, scored as it scores them; training
        // has no indel flag of its own.
        let cascade = self.stage1.is_some();
        let mapping = PileupMappingParams {
            snv_engine: cascade,
            call_indels: cascade,
            ..PileupMappingParams::default()
        };
        let SegmentColumns { segment, mut columns, held } =
            segment_columns(readers, segment, &mapping)?;
        if let Some(stage1) = self.stage1 {
            crate::call::process::score_on_cpu(&mut columns, stage1, |_, _| true)?;
            held.score_snv_clusters(&segment, &mut columns, &mapping);
        }
        let mut collected = self.empty();
        let mut examples = Examples { collected: &mut collected, contig: &segment.contig, keys };
        // The overlap columns stay in as neighbours, as in `call`, but only
        // the core offers candidates: the neighbouring segment offers the rest.
        map_surrounding(
            &mut columns,
            |before, column, after| {
                let pos = u64::from(column.pos);
                if !segment.is_core(pos) || !self.truth.claims(&segment.contig, pos) {
                    return Ok(());
                }
                self.snvs(truth, &mut examples, before, column, after)
            },
            "failed to extract training features, skipping",
        );
        Ok(collected)
    }

    fn snvs(
        &self,
        truth: &ContigTruth,
        examples: &mut Examples<'_>,
        before: Option<&PileupMetrics>,
        column: &PileupMetrics,
        after: Option<&PileupMetrics>,
    ) -> Result<()> {
        let pos = u64::from(column.pos);
        let ref_base = column.reference_base;
        for alt in &column.alts {
            let alt_base = alt.base;
            if ref_base == Base::Unknown || alt_base == Base::Unknown {
                continue;
            }
            if self.stage1.is_some() && !column.engine_scored(alt_base) {
                continue;
            }
            let Some(candidate) = column.alt_metrics(alt_base) else { continue };
            let label = truth.snps.contains(&PositionKey { pos, ref_base, alt_base });

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
            examples.add(model, features, label, pos)?;
        }
        Ok(())
    }
}

/// Where one segment's examples go.
struct Examples<'a> {
    collected: &'a mut Collected,
    contig: &'a SmolStr,
    keys: &'a mut ByModel<KeySource>,
}

impl Examples<'_> {
    /// Offer a candidate whose features could be computed and are all finite,
    /// and count any other as rejected.
    fn add(
        &mut self,
        model: MlModel,
        features: Result<Array2<f32>>,
        in_truth: bool,
        pos: u64,
    ) -> Result<()> {
        let features = match features {
            Ok(features) if features.iter().all(|value| value.is_finite()) => features,
            Ok(_) => {
                self.collected[model].reject();
                return Ok(());
            }
            Err(error) => {
                debug!(
                    model = model.name(),
                    error = format!("{error:#}"),
                    "No features for training candidate"
                );
                self.collected[model].reject();
                return Ok(());
            }
        };
        let row = features.as_slice().wrap_err("Feature row is not contiguous")?;
        self.collected[model].add_example(
            row,
            Label::of(in_truth),
            self.contig.clone(),
            pos,
            &mut self.keys[model],
        )
    }
}

/// The stretches of genome the segments cover, merged across their overlaps.
fn covered_regions(segments: &[ChunkRegion]) -> Result<Vec<RegionString>> {
    let mut spans: Vec<Region> = Vec::new();
    for segment in segments {
        match spans.last_mut() {
            Some(span)
                if span.contig == segment.contig && segment.start <= span.end.saturating_add(1) =>
            {
                span.end = span.end.max(segment.end);
            }
            _ => spans.push(segment.region.clone()),
        }
    }
    spans
        .iter()
        .map(|span| {
            let one_based = |pos: u64| {
                Pos0::try_from(pos)
                    .ok()
                    .and_then(|pos| pos.to_one_based().ok())
                    .wrap_err_with(|| format!("Segment {span} lies outside the addressable range"))
            };
            Ok(RegionString {
                chromosome: span.contig.clone(),
                start: Some(one_based(span.start)?),
                end: Some(one_based(span.end)?),
            })
        })
        .collect()
}

/// The SNP keys of every PASS record of the truth VCF in `regions`,
/// for each contig of `regions` the VCF declares.
#[instrument(level = "info", skip_all, fields(path = %path.display()))]
fn truth_keys(
    path: &Path,
    regions: &[RegionString],
    threads: usize,
) -> Result<HashMap<SmolStr, ContigTruth>> {
    let header =
        bcf::Reader::from_path(path).wrap_err("Failed to read the header")?.header().clone();
    let (declared, undeclared): (Vec<RegionString>, Vec<RegionString>) = regions
        .iter()
        .cloned()
        .partition(|region| header.name2rid(region.chromosome.as_bytes()).is_ok());
    let undeclared: BTreeSet<&str> = undeclared.iter().map(|r| r.chromosome.as_str()).collect();
    if !undeclared.is_empty() {
        warn!(
            contigs = ?undeclared,
            "The truth set does not declare these contigs; their candidates are not trained on"
        );
    }

    let mut contigs: HashMap<SmolStr, ContigTruth> =
        declared.iter().map(|r| (r.chromosome.clone(), ContigTruth::default())).collect();
    for_each_record_in_regions(path, &declared, threads, |record, header| {
        let Some(truth) = record
            .rid()
            .and_then(|rid| header.rid2name(rid).ok())
            .and_then(|name| std::str::from_utf8(name).ok())
            .and_then(|name| contigs.get_mut(name))
        else {
            return;
        };
        if !record.has_filter("PASS".as_bytes()) {
            return;
        }
        let Ok(pos) = u64::try_from(record.pos()) else { return };
        let alleles = record.alleles();
        let Some((&ref_allele, alts)) = alleles.split_first() else { return };
        for &alt_allele in alts {
            if let Some(key) = PositionKey::of_truth(pos, ref_allele, alt_allele) {
                truth.snps.insert(key);
            }
        }
    })?;
    info!(snps = contigs.values().map(|t| t.snps.len()).sum::<usize>(), "Loaded true variants");
    Ok(contigs)
}

impl PositionKey {
    /// The key of one ALT allele of a truth record once the suffix it shares
    /// with REF past the first base is stripped, which is how a multi-allelic
    /// record pads its shorter alleles.
    ///
    /// Only a single-base substitution has one; indels, MNPs, complex and
    /// symbolic alleles do not.
    fn of_truth(pos: u64, ref_allele: &[u8], alt_allele: &[u8]) -> Option<Self> {
        let (&ref_anchor, ref_rest) = ref_allele.split_first()?;
        let (&alt_anchor, alt_rest) = alt_allele.split_first()?;
        let (ref_rest, alt_rest) = strip_common_suffix(ref_rest, alt_rest);
        let (ref_base, alt_base) = (Base::from(ref_anchor), Base::from(alt_anchor));
        (ref_rest.is_empty()
            && alt_rest.is_empty()
            && ref_base != alt_base
            && ref_base != Base::Unknown
            && alt_base != Base::Unknown)
            .then_some(Self { pos, ref_base, alt_base })
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
/// each fit over `--threads` in a pool of its own. `prefix` names the forest
/// set in the feature importance files.
fn fit(
    params: &TrainModelParams,
    request: SamplingRequest,
    seed: u64,
    collected: &Collected,
    prefix: &str,
) -> Result<RastairFlatModel> {
    let seeds = forest_seeds(seed);
    let calculator = params.ml_features.get_calculator();
    let (counts, names) = (calculator.feature_num(), calculator.feature_names());

    let trained = ByModel::from_fn(|model| model).try_map(|model, _| -> Result<_> {
        let (forest, platt) = fit_model(model, &collected[model], params, request, seeds[model])
            .wrap_err_with(|| format!("Failed to train the {} model", model.name()))?;
        if let Some(dir) = params.feature_analytics.as_ref() {
            let path = dir.path().join(format!("{prefix}{}_feature_importances.csv", model.name()));
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
    info!(report = ?report, forests = prefix, "Trained models");

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
    request: SamplingRequest,
    seed: u64,
) -> Result<(RandomForest, PlattScaling)> {
    info!(seed, examples = data.len(), "Training model");

    let reservoir::Split { train, holdout } = reservoir::split(data, request, seed)?;
    let train_positives = train.labels.iter().filter(|&&l| l == Label::Positive.weight()).count();
    info!(
        training = train.labels.len(),
        positive = train_positives,
        negative = train.labels.len().saturating_sub(train_positives),
        holdout = holdout.labels.len(),
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
    forest.fit(&train.features.view(), &train.labels.view());

    let raw_scores = forest.predict(&holdout.features.view());
    let platt = fit_platt_scaling(
        raw_scores.as_slice().wrap_err("Holdout scores are not contiguous")?,
        holdout.labels.as_slice().wrap_err("Holdout labels are not contiguous")?,
    )?;

    Ok((forest, platt))
}

/// Fit Platt scaling parameters A and B so that
/// `P(y=1|f) = 1 / (1 + exp(A*f + B))` is a well-calibrated probability.
///
/// Uses Newton's method with backtracking line search and Bayesian-smoothed
/// targets, following Lin, Lin, and Weng (2007).
fn fit_platt_scaling(scores: &[f64], labels: &[f64]) -> Result<PlattScaling> {
    const MAX_ITER: usize = 100;
    const MIN_STEP: f64 = 1e-10;
    const SIGMA: f64 = 1e-12;

    let samples = || scores.iter().copied().zip(labels.iter().map(|&y| y > 0.5));
    let (n_pos, n) = samples().fold((0.0_f64, 0.0_f64), |(pos, all), (_, y)| {
        (if y { pos + 1.0 } else { pos }, all + 1.0)
    });
    let n_neg = n - n_pos;
    ensure!(
        n_pos > 0.0 && n_neg > 0.0,
        "Platt calibration needs both classes, but the holdout has {n_pos} positive and \
         {n_neg} negative examples"
    );

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

    Ok(PlattScaling { a, b })
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

    for example in data.examples() {
        let row = data.row_of(example).wrap_err("Training example has no feature row")?;
        let pos = example.pos.saturating_add(1);
        write!(writer, "{}\t{pos}\t{}", example.chrom, example.label.weight())?;
        for v in row {
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
    use clap::Parser as _;

    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        model: ModelParameters,
    }

    fn model_parameters(args: &[&str]) -> ModelParameters {
        Cli::try_parse_from(std::iter::once("train").chain(args.iter().copied())).unwrap().model
    }

    #[derive(clap::Parser)]
    struct TrainCli {
        #[command(flatten)]
        train: TrainModelParams,
    }

    /// The model is written after collecting and fitting everything, so a
    /// path that cannot take it must fail before any of that.
    #[test]
    fn an_output_directory_is_refused_before_any_work() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let output = dir.path().to_str().wrap_err("temporary path is not UTF-8")?;
        let params = TrainCli::try_parse_from([
            "train",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "truth.vcf.gz",
            "--output",
            output,
        ])?
        .train;
        let error = create_output_dirs(&params).expect_err("a directory is not a model file");
        assert!(error.to_string().contains("is a directory"), "{error}");
        Ok(())
    }

    /// `call` scores the columns `ml train` labels, so a candidate's shared
    /// features are the same on both paths.
    #[test]
    fn training_and_calling_see_the_same_common_features() -> Result<()> {
        use crate::metrics::ml::features::shared::{
            CommonFeatures, CommonSectionA, CommonSectionB,
        };

        #[derive(clap::Parser)]
        struct CallCli {
            #[command(flatten)]
            call: crate::CallParams,
        }
        let call = CallCli::try_parse_from([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=bacteriophage_lambda_CpG:1001-4000",
        ])?
        .call;
        let common = |column: &PileupMetrics, alt: Base| -> Option<Vec<f32>> {
            let features = CommonFeatures::extract(&column.alt_metrics(alt)?);
            let mut row = CommonSectionA::from_common(&features).as_row().to_vec();
            row.extend_from_slice(CommonSectionB::from_common(&features).as_row());
            Some(row)
        };

        let mut readers = call.segments.pileup_readers()?;
        let segment = readers.segments(10_000, 0)?.next().wrap_err("no segment")?;
        let training =
            segment_columns(&mut readers, &segment, &PileupMappingParams::default())?.columns;
        let mut calling = segment_columns(&mut readers, &segment, &call.pileup_mapping())?.columns;
        calling.retain(|column| call.record_filters.pre_filter(column));

        let training: HashMap<_, _> = training.iter().map(|column| (column.pos, column)).collect();
        let mut compared = 0;
        for column in &calling {
            let trained = training.get(&column.pos).wrap_err("call column missing in training")?;
            for alt in column.alts() {
                assert_eq!(common(trained, alt), common(column, alt), "{} {alt}", column.pos);
                compared += 1;
            }
        }
        assert!(compared > 100, "only {compared} candidates compared");
        Ok(())
    }

    /// A bare command line trains the bundled model's recipe.
    #[test]
    fn the_defaults_are_the_bundled_recipe() {
        let params = model_parameters(&[]);
        assert_eq!((params.max_features, params.max_depth, params.min_samples_leaf), (4, 40, 10));
        assert_eq!(params.sampling_request(), ByLabel { positive: 8_000, negative: 200_000 });
    }

    /// A retrain reproduces a model only if every forest gets the seed it got
    /// before, so the draw order is pinned.
    #[test]
    fn forest_seeds_are_drawn_in_training_order() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let expected: Vec<u64> = (0..MlModel::COUNT).map(|_| rng.random()).collect();
        let seeds = forest_seeds(7);
        let order = [MlModel::Cpg, MlModel::DenovoCpg, MlModel::Others];
        assert_eq!(order.map(|model| seeds[model]).to_vec(), expected);
    }

    fn truth_key(ref_allele: &str, alt_allele: &str) -> Option<PositionKey> {
        PositionKey::of_truth(10, ref_allele.as_bytes(), alt_allele.as_bytes())
    }

    /// `REF=AT ALT=GT` is an A>G substitution padded to the record's REF.
    #[test]
    fn a_padded_substitution_keys_an_snv() {
        let snv = |ref_base, alt_base| Some(PositionKey { pos: 10, ref_base, alt_base });
        assert_eq!(truth_key("AT", "GT"), snv(Base::A, Base::G));
        assert_eq!(truth_key("C", "T"), snv(Base::C, Base::T));
    }

    #[test]
    fn indel_complex_and_symbolic_alleles_key_nothing() {
        for (ref_allele, alt_allele) in [
            ("ATT", "AT"),
            ("A", "ATG"),
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

    /// A calibration fit on one class would be the identity, which inverts
    /// the forest's scores; no model may be written with it.
    #[test]
    fn a_holdout_of_one_class_cannot_be_calibrated() {
        assert!(fit_platt_scaling(&[0.9, 0.1], &[0.0, 0.0]).is_err());
        assert!(fit_platt_scaling(&[0.9, 0.1], &[1.0, 1.0]).is_err());
        let platt = fit_platt_scaling(&[0.9, 0.8, 0.2, 0.1], &[1.0, 1.0, 0.0, 0.0]).unwrap();
        assert!(platt.a < 0.0, "higher scores must mean higher probability: {platt:?}");
    }

    /// A row biosphere cannot fit on is counted, not collected.
    #[test]
    fn a_non_finite_or_missing_row_is_rejected_and_counted() -> Result<()> {
        let request = SamplingRequest { positive: 10, negative: 10 };
        let mut collected = ByModel::from_fn(|_| TrainingData::for_request(2, request));
        let mut keys = ByModel::from_fn(|model| KeySource::for_segment(0, 0, model));
        let contig = SmolStr::from("chr1");
        let mut examples = Examples { collected: &mut collected, contig: &contig, keys: &mut keys };
        let row = |values: [f32; 2]| Ok(Array2::from_shape_vec((1, 2), values.to_vec())?);

        examples.add(MlModel::Cpg, row([1.0, 2.0]), true, 1)?;
        examples.add(MlModel::Cpg, row([1.0, f32::INFINITY]), true, 2)?;
        examples.add(MlModel::Cpg, row([f32::NAN, 2.0]), false, 3)?;
        examples.add(MlModel::Cpg, Err(eyre!("no neighbour")), false, 4)?;

        assert_eq!(collected[MlModel::Cpg].len(), 1);
        assert_eq!(collected[MlModel::Cpg].rejected(), 3);
        assert_eq!(collected[MlModel::Others].rejected(), 0);
        Ok(())
    }

    /// `--export-features` is joined against VCFs, whose positions are 1-based.
    #[test]
    fn exported_features_carry_1_based_positions() -> Result<()> {
        let mut data = TrainingData::for_request(2, SamplingRequest { positive: 1, negative: 1 });
        let mut keys = KeySource::for_segment(0, 0, MlModel::Cpg);
        data.add_example(&[0.5, 1.5], Label::Positive, "chr1".into(), 99, &mut keys)?;
        let dir = tempfile::tempdir()?;
        export_features_tsv(&data, "cpg", &["a", "b"], dir.path())?;
        let tsv = std::fs::read_to_string(dir.path().join("cpg_features.tsv"))?;
        assert_eq!(tsv, "chrom\tpos\tlabel\ta\tb\nchr1\t100\t1\t0.5\t1.5\n");
        Ok(())
    }

    /// Neighbouring segments share their overlap, and each of its candidates
    /// must be offered once, by the segment whose core holds it.
    #[test]
    fn overlapping_segments_offer_every_candidate_once() -> Result<()> {
        let truth = Truth {
            contigs: HashMap::from([("bacteriophage_lambda_CpG".into(), ContigTruth::default())]),
            confident: None,
        };
        let calculator = MlFeatureSet::Standard.get_calculator();
        let collector = SegmentCollector {
            truth: &truth,
            stage1: None,
            calculator: &*calculator,
            feature_num: calculator.feature_num(),
            request: model_parameters(&[]).sampling_request(),
        };
        let mut params = ReaderParams::test_data();
        params.regions = Some(crate::utils::CliRegionInput::from_region(
            "bacteriophage_lambda_CpG:1001-4000".parse()?,
        ));
        let offered = |max_length: u64, overlap: u64| -> Result<u64> {
            let mut readers = params.pileup_readers()?;
            let segments: Vec<ChunkRegion> = readers.segments(max_length, overlap)?.collect();
            let mut offered = 0;
            for (index, segment) in segments.iter().enumerate() {
                let collected = collector.segment(
                    &mut readers,
                    segment,
                    &mut ByModel::from_fn(|model| KeySource::for_segment(0, index, model)),
                )?;
                offered += collected.iter().map(|(_, data)| data.seen().negative).sum::<u64>();
            }
            Ok(offered)
        };

        let whole = offered(10_000, 0)?;
        assert!(whole > 0);
        assert_eq!(offered(1_000, 200)?, whole);
        Ok(())
    }

    /// With the first stage's forests, a candidate is offered only when the
    /// SNV engine scored it after them: fewer than every alt, more than none.
    #[test]
    fn the_engine_forests_see_only_what_the_cascade_rescores() -> Result<()> {
        let truth = Truth {
            contigs: HashMap::from([("bacteriophage_lambda_CpG".into(), ContigTruth::default())]),
            confident: None,
        };
        let calculator = MlFeatureSet::Standard.get_calculator();
        let stage1 = MachineLearningParams::default().init(1)?;
        let params = ReaderParams::test_data();
        let offered = |stage1: Option<&MachineLearning>| -> Result<u64> {
            let collector = SegmentCollector {
                truth: &truth,
                stage1,
                calculator: &*calculator,
                feature_num: calculator.feature_num(),
                request: model_parameters(&[]).sampling_request(),
            };
            let mut readers = params.pileup_readers()?;
            let segments: Vec<ChunkRegion> = readers.segments(100_000, 0)?.collect();
            let mut offered = 0;
            for (index, segment) in segments.iter().enumerate() {
                let collected = collector.segment(
                    &mut readers,
                    segment,
                    &mut ByModel::from_fn(|model| KeySource::for_segment(0, index, model)),
                )?;
                offered += collected.iter().map(|(_, data)| data.seen().negative).sum::<u64>();
            }
            Ok(offered)
        };
        let (every_alt, rescored) = (offered(None)?, offered(Some(&stage1))?);
        assert!(rescored > 0, "the engine scored nothing on the test BAM");
        assert!(rescored < every_alt, "{rescored} of {every_alt} alts were rescored");
        Ok(())
    }

    fn chunk(contig: &str, start: u64, end: u64) -> ChunkRegion {
        ChunkRegion {
            region: Region { contig: contig.into(), start, end },
            last_position: 1_000_000,
            overlap_start: 0,
            overlap_end: 0,
        }
    }

    /// Without `-l` the segments cover the whole BAM, and every contig they
    /// reach must be labelled against its own truth.
    #[test]
    fn truth_is_read_for_every_stretch_the_segments_cover() -> Result<()> {
        let segments = [chunk("chr1", 99, 299), chunk("chr1", 200, 499), chunk("chr2", 0, 99)];
        assert_eq!(
            covered_regions(&segments)?,
            vec!["chr1:100-500".parse::<RegionString>()?, "chr2:1-100".parse()?]
        );
        Ok(())
    }

    /// Write a bgzipped, indexed truth VCF declaring `chr1` and `chr2`.
    fn truth_vcf(dir: &Path, records: &[u8]) -> Result<std::path::PathBuf> {
        use std::io::Write as _;
        let path = dir.join("truth.vcf.gz");
        let mut vcf = rust_htslib::bgzf::Writer::from_path(&path)?;
        vcf.write_all(
            b"##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
              ##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        )?;
        vcf.write_all(records)?;
        drop(vcf);
        bcf::index::build(&path, None, 1, bcf::index::Type::Csi(14))?;
        Ok(path)
    }

    /// A truth variant on a region's first base labels its candidate positive.
    #[test]
    fn a_truth_variant_on_the_first_base_of_a_region_is_loaded() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = truth_vcf(
            dir.path(),
            b"chr1\t100\t.\tA\tG\t.\tPASS\t.\nchr1\t301\t.\tC\tT\t.\tPASS\t.\n",
        )?;

        let contigs = truth_keys(&path, &covered_regions(&[chunk("chr1", 99, 299)])?, 1)?;
        let first = PositionKey { pos: 99, ref_base: Base::A, alt_base: Base::G };
        let chr1 = contigs.get("chr1").wrap_err("chr1 is declared")?;
        assert_eq!(chr1.snps, HashSet::from([first]));
        Ok(())
    }

    /// Positions repeat across contigs, so a truth variant must label only
    /// the contig it is on.
    #[test]
    fn a_truth_variant_labels_only_its_own_contig() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = truth_vcf(
            dir.path(),
            b"chr1\t100\t.\tA\tG\t.\tPASS\t.\nchr1\t200\t.\tA\tAT\t.\tPASS\t.\n",
        )?;

        let segments = [chunk("chr1", 0, 999), chunk("chr2", 0, 999)];
        let contigs = truth_keys(&path, &covered_regions(&segments)?, 1)?;
        let chr1 = contigs.get("chr1").wrap_err("chr1 is declared")?;
        let chr2 = contigs.get("chr2").wrap_err("chr2 is declared")?;
        assert!(chr1.snps.contains(&PositionKey { pos: 99, ref_base: Base::A, alt_base: Base::G }));
        assert_eq!(chr1.snps.len(), 1, "the insertion is not an SNV");
        assert!(chr2.snps.is_empty());
        Ok(())
    }

    /// A contig the truth VCF does not declare is one it makes no claim
    /// about, so its candidates must not all become negatives.
    #[test]
    fn a_contig_the_truth_does_not_declare_is_not_trained_on() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = truth_vcf(dir.path(), b"chr1\t100\t.\tA\tG\t.\tPASS\t.\n")?;

        let segments = [chunk("chr1", 0, 999), chunk("chrM", 0, 999)];
        let contigs = truth_keys(&path, &covered_regions(&segments)?, 1)?;
        assert!(contigs.contains_key("chr1"));
        assert!(!contigs.contains_key("chrM"));
        Ok(())
    }
}
