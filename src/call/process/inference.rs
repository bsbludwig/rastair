//! A dedicated thread for GPU inference.
//!
//! One thread owns one set of forests, and whatever batches are already queued
//! when it wakes ride along in the same dispatch. Forests per worker would
//! give the GPU `regions × models` dispatches of a few hundred rows each.
//! A worker blocks while the GPU has its rows.

use super::calc_ml::{MlRows, MlScores, apply_ml_scores, extract_ml_rows, submit_and_collect};
use crate::metrics::{
    PileupMetrics,
    ml::types::{ByModel, GpuRastairModel, MachineLearning},
};
use crate::runtime::{fault_injection, threads};
use color_eyre::eyre::{Report, Result, WrapErr as _, eyre};
use crossbeam_channel::{Receiver, Sender, bounded};
use ndarray::{Array2, s};
use seqair_types::Base;
use std::{
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    time::{Duration, Instant},
};
use tracing::{debug, error, trace, warn};

/// A dispatch failed for a whole round, so every job in it hears about it. A
/// [`Report`] is not `Clone`, and this is only ever read to log and fall back.
type Scored = Result<MlScores, Arc<Report>>;

/// The GPU-owning thread, plus the channel that feeds it.
pub struct InferenceStage {
    /// `None` only while [`Drop`] closes the channel to stop the thread.
    jobs: Option<Sender<Job>>,
    thread: Option<threads::Thread<()>>,
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

impl InferenceStage {
    /// Take ownership of `gpu` and start scoring.
    ///
    /// `workers` sizes the queue: each worker has at most one batch queued,
    /// and two per worker leaves headroom.
    pub fn spawn(gpu: GpuRastairModel, workers: usize) -> Result<Self> {
        let (jobs, incoming) = bounded(workers.saturating_mul(2).max(1));
        let server = Server { gpu, incoming, tally: Tally::default(), started: Instant::now() };
        // The run can do without the GPU: after a crash, regions fall back to
        // the CPU, like after any other GPU failure.
        let thread = threads::spawn_recovering("inference", server, Server::run, Server::recover)
            .wrap_err("Failed to start the GPU inference thread")?;

        Ok(Self { jobs: Some(jobs), thread: Some(thread), failed: AtomicBool::new(false) })
    }

    /// Extract, score and apply a region's ML candidates. Blocks only while
    /// the queue is full and while the GPU has the rows.
    fn score(
        &self,
        pileups: &mut [PileupMetrics],
        ml: &MachineLearning,
        only: impl Fn(&PileupMetrics, Base) -> bool,
    ) -> Result<()> {
        if pileups.is_empty() {
            return Ok(());
        }

        let (rows, targets) = extract_ml_rows(pileups, ml, only)?.split()?;
        let (scored, receiver) = bounded(1);
        self.jobs
            .as_ref()
            .ok_or_else(|| eyre!("The GPU inference thread is shutting down"))?
            .send(Job { rows, scored })
            .map_err(|_| eyre!("The GPU inference thread stopped before it could be sent work"))?;

        let scores = receiver
            .recv()
            .map_err(|_| eyre!("The GPU inference thread stopped before it returned scores"))?
            .map_err(|error| eyre!("GPU inference failed: {error:#}"))?;
        // A short reply fails here like a failed dispatch, so `score_on_gpu`
        // hands the region to the CPU.
        apply_ml_scores(pileups, &targets, &scores)
    }
}

/// Whether the inference thread scored a region's pileups.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[must_use = "`LeftToCaller` means the pileups still need scoring"]
pub enum GpuScoring {
    Scored,
    /// The run has no inference thread, or it failed: the caller scores
    /// them on its own.
    LeftToCaller,
}

/// Score `pileups` on the inference thread, if this run has one that still
/// works.
///
/// The first failure retires the stage for the rest of the run, and is the
/// only one logged: the regions in flight at the time all fail with it.
pub fn score_on_gpu(
    pileups: &mut [PileupMetrics],
    ml: &MachineLearning,
    only: impl Fn(&PileupMetrics, Base) -> bool,
) -> GpuScoring {
    let Some(stage) = ml.inference.as_ref() else {
        return GpuScoring::LeftToCaller;
    };
    if stage.failed.load(Ordering::Relaxed) {
        return GpuScoring::LeftToCaller;
    }

    match stage.score(pileups, ml, only) {
        Ok(()) => GpuScoring::Scored,
        Err(error) => {
            if !stage.failed.swap(true, Ordering::Relaxed) {
                warn!(
                    error = format!("{error:#}"),
                    "GPU inference failed, scoring the rest of the run on the CPU"
                );
            }
            GpuScoring::LeftToCaller
        }
    }
}

impl Drop for InferenceStage {
    fn drop(&mut self) {
        // Closing the channel is what ends `serve`'s loop (or the drain after a
        // crash); the forests are then dropped on the inference thread rather
        // than on a rayon worker during TLS teardown, which upsets Metal.
        self.jobs = None;
        if let Some(thread) = self.thread.take()
            && let Err(error) = thread.join()
        {
            error!(error = format!("{error:#}"), "The GPU inference thread failed");
        }
    }
}

/// What the inference thread owns.
struct Server {
    gpu: GpuRastairModel,
    incoming: Receiver<Job>,
    tally: Tally,
    started: Instant,
}

impl Server {
    fn run(&mut self) {
        serve(|rows| submit_and_collect(rows, &self.gpu), &self.incoming, &mut self.tally);
        self.tally.report(self.started.elapsed());
    }

    /// After a panic, `gpu` is not used again.
    fn recover(self) {
        // Jobs still queued are only dropped once both channel ends are gone,
        // and the sending end lives as long as the stage. Without draining,
        // every worker waiting on a queued job would wait forever. Dropping a
        // job drops its reply sender, so its worker gets an error and falls
        // back to the CPU like after any other GPU failure.
        self.incoming.iter().for_each(drop);
        self.tally.report(self.started.elapsed());
    }
}

fn serve(
    mut score: impl FnMut(&MlRows) -> Result<MlScores>,
    incoming: &Receiver<Job>,
    tally: &mut Tally,
) {
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
        let scored = score(&round.rows)
            .and_then(|scores| {
                fault_injection::fault_point(fault_injection::FaultPoint::GpuDispatch)?;
                Ok(scores)
            })
            .map_err(Arc::new);
        tally.busy += scoring.elapsed();

        let failed = scored.is_err();
        round.reply(scored);
        if failed {
            // biosphere leaves a forest busy after a failed collect, or with a
            // handle the error dropped uncollected, and panics on its next
            // submit. Dropping the queued jobs sends their workers to the CPU.
            incoming.iter().for_each(drop);
            break;
        }
    }
}

/// What the stage did, so whether it is the bottleneck is measured rather
/// than reasoned about.
#[derive(Default)]
struct Tally {
    rounds: u64,
    jobs: u64,
    rows: u64,
    /// Time inside `submit_and_collect`, the whole of what a waiting worker
    /// can be blocked on besides the queue.
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
    replies: Vec<Reply>,
}

/// Where a job's scores go: its rows per model, in the round's order.
struct Reply {
    rows: ByModel<usize>,
    to: Sender<Scored>,
}

impl Round {
    fn new(jobs: Vec<Job>) -> Self {
        let (mut per_job, senders): (Vec<MlRows>, Vec<_>) =
            jobs.into_iter().map(|job| (job.rows, job.scored)).unzip();
        let shares = per_job.iter().map(|rows| ByModel::from_fn(|m| rows[m].nrows()));
        let replies = shares.zip(senders).map(|(rows, to)| Reply { rows, to }).collect();

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

        for Reply { rows, to } in self.replies {
            let slice = scored.as_ref().map(|scores| {
                MlScores::from_fn(|m| {
                    let from = start[m];
                    start[m] = from + rows[m];
                    scores[m].get(from..start[m]).unwrap_or_default().to_vec()
                })
            });
            // A worker that gave up (its region failed elsewhere) drops the
            // receiver; there is nobody left to tell, and that is fine.
            let sent = to.send(slice.map_err(Arc::clone));
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

#[cfg(test)]
mod tests {
    use super::*;

    fn job() -> (Job, Receiver<Scored>) {
        let (scored, receiver) = bounded(1);
        (Job { rows: MlRows::from_fn(|_| Array2::zeros((1, 1))), scored }, receiver)
    }

    /// biosphere keeps a forest busy after a failed collect, so the GPU must
    /// not be asked again: a job queued behind a failed round is dropped, and
    /// its worker falls back to the CPU.
    #[test]
    fn a_failed_round_retires_the_gpu() {
        let (jobs, incoming) = bounded(2);
        let (first, first_reply) = job();
        let (second, second_reply) = job();
        jobs.send(first).unwrap();
        let mut queue_second = Some((jobs, second));
        let mut calls = 0;
        serve(
            |_| {
                calls += 1;
                // Queue the next job while the first is on the GPU, and close
                // the channel so `serve` ends.
                if let Some((jobs, second)) = queue_second.take() {
                    jobs.send(second).unwrap();
                }
                Err(eyre!("readback timed out"))
            },
            &incoming,
            &mut Tally::default(),
        );
        assert_eq!(calls, 1);
        assert!(first_reply.recv().unwrap().is_err());
        assert!(second_reply.recv().is_err(), "the queued job is dropped unscored");
    }
}
