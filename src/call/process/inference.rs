//! A dedicated thread for GPU inference.
//!
//! Every worker used to fork its own copy of the five forests and dispatch its
//! own region's rows, so the GPU saw `regions × models` dispatches carrying a
//! few hundred rows each while nine threads shared one device per model. Here
//! one thread owns one set of forests, and whatever batches are already queued
//! when it wakes ride along in the same dispatch.
//!
//! Submitting and waiting are separate ([`InferenceStage::submit`] returns a
//! [`Ticket`]) so a caller can do other work while the GPU has its rows.

use super::calc_ml::{
    MlRows, MlScores, MlTargets, apply_ml_scores, extract_ml_rows, submit_and_collect,
};
use crate::metrics::{
    PileupMetrics,
    ml::types::{ByModel, GpuRastairModel, MachineLearning},
};
use color_eyre::eyre::{Report, Result, eyre};
use crossbeam_channel::{Receiver, Sender, bounded};
use ndarray::{Array2, s};
use seqair_types::Probability;
use std::{
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    thread,
    time::{Duration, Instant},
};
use tracing::{debug, trace, warn};

/// A dispatch failed for a whole round, so every job in it hears about it. A
/// [`Report`] is not `Clone`, and this is only ever read to log and fall back.
type Scored = Result<MlScores, Arc<Report>>;

/// The GPU-owning thread, plus the channel that feeds it.
pub struct InferenceStage {
    /// `None` only while [`Drop`] closes the channel to stop the thread.
    jobs: Option<Sender<Job>>,
    thread: Option<thread::JoinHandle<()>>,
    /// Set by the first region the GPU fails to score. From then on every
    /// region goes straight to the CPU: a GPU that has failed once tends to
    /// keep failing, and each retry costs a full collect timeout serialised
    /// through this one thread, plus a second feature extraction for the
    /// fallback.
    failed: AtomicBool,
}

struct Job {
    rows: MlRows,
    scored: Sender<Scored>,
}

/// Rows handed to the stage, not yet scored.
#[must_use = "a submitted batch leaves the region unscored until it is applied"]
pub struct Ticket {
    targets: MlTargets,
    scored: Receiver<Scored>,
}

impl InferenceStage {
    /// Take ownership of `gpu` and start scoring.
    ///
    /// `workers` sizes the queue. Two per worker is enough for every worker to
    /// have one batch in flight and one waiting, which is the most that can be
    /// outstanding even once callers pipeline.
    pub fn spawn(gpu: GpuRastairModel, workers: usize) -> Result<Self> {
        let (jobs, incoming) = bounded(workers.saturating_mul(2).max(1));
        let thread = thread::Builder::new()
            .name("inference".into())
            .spawn(move || run(&gpu, &incoming))
            .map_err(|error| eyre!("Failed to start the GPU inference thread: {error}"))?;

        Ok(Self { jobs: Some(jobs), thread: Some(thread), failed: AtomicBool::new(false) })
    }

    /// Extract, score and apply a region's ML candidates.
    ///
    /// Blocks while the GPU has the rows. Splitting that into
    /// [`Self::submit`] and [`Ticket::apply`] is what will let a caller
    /// build its next region in the meantime.
    fn score(
        &self,
        pileups: &mut [PileupMetrics],
        ml: &MachineLearning,
        score_indels: bool,
    ) -> Result<()> {
        let Some(model) = ml.model.as_ref() else {
            return Ok(());
        };
        if pileups.is_empty() {
            return Ok(());
        }

        let (rows, targets) = extract_ml_rows(pileups, ml, model, score_indels)?.split()?;
        self.submit(rows, targets)?.apply(pileups, ml.threshold)
    }

    /// Queue `rows` for scoring. Blocks only while the queue is full.
    pub fn submit(&self, rows: MlRows, targets: MlTargets) -> Result<Ticket> {
        let (scored, receiver) = bounded(1);
        self.jobs
            .as_ref()
            .ok_or_else(|| eyre!("The GPU inference thread is shutting down"))?
            .send(Job { rows, scored })
            .map_err(|_| eyre!("The GPU inference thread stopped before it could be sent work"))?;

        Ok(Ticket { targets, scored: receiver })
    }
}

/// Score `pileups` on the inference thread, if this run has one that still works.
///
/// `None` says there is no GPU stage, or it has already failed once, so the
/// caller has to score on its own thread; `Some(Err(_))` says the GPU tried
/// and failed just now, and the caller should fall back the same way. Only the
/// first failure is reported this way; after it the stage is retired for the
/// rest of the run.
pub fn score_on_gpu(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    score_indels: bool,
) -> Option<Result<()>> {
    let stage = ml.inference.as_ref()?;
    if stage.failed.load(Ordering::Relaxed) {
        return None;
    }

    let result = stage.score(pileups, ml, score_indels);
    if result.is_err() && !stage.failed.swap(true, Ordering::Relaxed) {
        warn!("GPU inference failed once; scoring the rest of the run on the CPU");
    }
    Some(result)
}

impl Ticket {
    /// Block until the scores come back, then write them into `pileups`.
    pub fn apply(self, pileups: &mut [PileupMetrics], threshold: Probability) -> Result<()> {
        let scores = self
            .scored
            .recv()
            .map_err(|_| eyre!("The GPU inference thread stopped before it returned scores"))?
            .map_err(|error| eyre!("GPU inference failed: {error:#}"))?;

        apply_ml_scores(pileups, &self.targets, &scores, threshold);
        Ok(())
    }
}

impl Drop for InferenceStage {
    fn drop(&mut self) {
        // Closing the channel is what ends `run`'s loop; the forests are then
        // dropped on the inference thread rather than on a rayon worker during
        // TLS teardown, which is what used to upset Metal.
        self.jobs = None;
        if let Some(thread) = self.thread.take()
            && thread.join().is_err()
        {
            tracing::error!("The GPU inference thread panicked");
        }
    }
}

fn run(gpu: &GpuRastairModel, incoming: &Receiver<Job>) {
    let mut tally = Tally::default();
    let started = Instant::now();

    loop {
        let waiting = Instant::now();
        let Ok(first) = incoming.recv() else { break };
        tally.idle += waiting.elapsed();

        let mut jobs = vec![first];
        // Take whatever else is already waiting, but never wait for more: a
        // timer here would add exactly the latency this thread exists to
        // remove. How much rides along is therefore set by how backed up the
        // workers are, which is the right signal.
        jobs.extend(incoming.try_iter());

        let round = Round::new(jobs);
        tally.rounds += 1;
        tally.jobs += round.replies.len() as u64;
        tally.rows += round.rows.iter().map(|(_, rows)| rows.nrows() as u64).sum::<u64>();
        trace!(jobs = round.replies.len(), "Scoring a round");

        let scoring = Instant::now();
        let scored = submit_and_collect(&round.rows, gpu).map_err(Arc::new);
        tally.busy += scoring.elapsed();

        round.reply(scored);
    }

    tally.report(started.elapsed());
}

/// What the stage did, so the next question — is it the bottleneck, and is a
/// worker's wait long enough to be worth pipelining around? — has an answer
/// that is measured rather than reasoned about.
#[derive(Default)]
struct Tally {
    rounds: u64,
    jobs: u64,
    rows: u64,
    /// Time inside `submit_and_collect`. This is the whole of what a waiting
    /// worker can be blocked on, so it bounds what pipelining could recover.
    busy: Duration,
    /// Time blocked in `recv` with nothing to do.
    idle: Duration,
}

impl Tally {
    fn report(&self, wall: Duration) {
        let per = |n: u64| n.checked_div(self.rounds);
        debug!(
            rounds = self.rounds,
            jobs = self.jobs,
            rows = self.rows,
            jobs_per_round = per(self.jobs).unwrap_or(0),
            rows_per_round = per(self.rows).unwrap_or(0),
            busy_ms = self.busy.as_millis(),
            idle_ms = self.idle.as_millis(),
            wall_ms = wall.as_millis(),
            busy_percent = percent(self.busy, wall),
            "GPU inference thread finished"
        );
    }
}

fn percent(part: Duration, whole: Duration) -> u64 {
    let whole = whole.as_micros();
    if whole == 0 {
        return 0;
    }
    u64::try_from(part.as_micros().saturating_mul(100) / whole).unwrap_or(u64::MAX)
}

/// One dispatch's worth of work: the queued jobs' rows in one set of arrays,
/// with what each contributed so the scores can be handed back apart.
struct Round {
    rows: MlRows,
    replies: Vec<(ByModel<usize>, Sender<Scored>)>,
}

impl Round {
    fn new(jobs: Vec<Job>) -> Self {
        let (mut per_job, senders): (Vec<MlRows>, Vec<_>) =
            jobs.into_iter().map(|job| (job.rows, job.scored)).unzip();
        let shares = per_job.iter().map(|rows| ByModel::from_fn(|m| rows[m].nrows()));
        let replies = shares.zip(senders).collect();

        // One job is the common case and needs no copy at all.
        let rows = match per_job.len() {
            1 => per_job.pop().unwrap_or_else(|| MlRows::from_fn(|_| Array2::zeros((0, 0)))),
            _ => concat(&per_job),
        };

        Self { rows, replies }
    }

    /// Hand each job the slice of `scored` its own rows produced.
    fn reply(self, scored: Scored) {
        let mut start = ByModel::from_fn(|_| 0usize);

        for (share, sender) in self.replies {
            let slice = scored.as_ref().map(|scores| {
                MlScores::from_fn(|m| {
                    let from = start[m];
                    start[m] = from + share[m];
                    scores[m].get(from..start[m]).unwrap_or_default().to_vec()
                })
            });
            // A worker that gave up (its region failed elsewhere) drops the
            // receiver; there is nobody left to tell, and that is fine.
            let sent = sender.send(slice.map_err(Arc::clone));
            if sent.is_err() {
                trace!("Nobody was waiting for a scored batch");
            }
        }
    }
}

fn concat(per_job: &[MlRows]) -> MlRows {
    MlRows::from_fn(|model| {
        let total = per_job.iter().map(|rows| rows[model].nrows()).sum();
        let columns = per_job.first().map_or(0, |rows| rows[model].ncols());
        let mut merged = Array2::zeros((total, columns));

        let mut at = 0;
        for rows in per_job {
            let source = &rows[model];
            merged.slice_mut(s![at..at + source.nrows(), ..]).assign(source);
            at += source.nrows();
        }
        merged
    })
}
