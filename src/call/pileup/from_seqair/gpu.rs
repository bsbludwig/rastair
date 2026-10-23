//! The pair-HMM on the GPU (`--gpu`): a region's clusters queued into
//! launches, each run by the worker's own GPU thread while the worker prepares
//! the next.
//!
//! The calls equal the CPU path's: the kernel is bit-identical to compair's
//! CPU kernels, and a pair whose `f32` score falls under its precision floor is
//! rescored in `f64` on the CPU (`GpuScore::or_rescore`).
//!
//! - `pair_hmm::prepare` builds a cluster's haplotypes, reads and bands;
//!   `Prepared::align` runs the CPU kernel, `Prepared::push` adds the pairs to
//!   a `GpuPairs` plan and `Prepared::finish_gpu` turns a launch's scores back
//!   through the same `rows`. Keep the arms sharing everything but the kernel
//!   call.
//! - One GPU thread per worker (`Lane`, a `thread_local`): the worker sends a
//!   plan, keeps preparing, and collects later. The aligner's buffers live on
//!   that thread and grow to the largest launch once. The wait is compair's
//!   `Wait::Block`; "Buffer with 'scores' label is still mapped" on every
//!   launch after the first means a lost mapping callback there. A failed
//!   launch rebuilds the lane's aligner.
//! - A launch is at most `MAX_LAUNCH_PAIRS` pairs (a single larger cluster
//!   launches alone), so a worker holds at most two launches' reads. The
//!   per-region `Batch` is flushed after each sub-segment's cluster loop,
//!   before the next sub-segment's reads load; a batch dropped with a launch in
//!   flight drains it in `Drop`, or the next region would collect it.
//! - Every failure (no adapter, a plan too large, a lost thread, a failed
//!   launch) falls back to the CPU kernel for the clusters concerned, with a
//!   `warn!`. Without `--gpu` this module is never entered, and `cargo test`
//!   does not exercise it.
//! - `RASTAIR_LOG=rastair::call::pileup::from_seqair=debug` logs a line per
//!   launch with running totals; `grep -c "Pair-HMM GPU launch"` checks that
//!   the GPU path ran. Score every arm of a comparison with `--gpu` or every
//!   arm without it: the two are not bitwise identical.

use super::{ClusterPass, cluster::assign_scores};
use crate::{
    call::pair_hmm::{HaplotypeScratch, Prepared},
    metrics::PileupMetrics,
};
use compair::gpu::{GpuAligner, GpuContext, GpuError, GpuPairs, GpuScore, KernelOptions, Wait};
use std::{
    cell::OnceCell,
    mem,
    ops::Range,
    sync::{
        Arc, OnceLock,
        atomic::{AtomicU64, Ordering},
        mpsc,
    },
    thread,
};
use tracing::{debug, info, warn};

/// Pairs a launch holds at most. A launch costs ~1.5 ms over its pairs,
/// so smaller launches lose that; larger ones hold more prepared reads
/// in memory, twice over while one runs and the next fills. 4,096 and
/// 8,192 time the same on chr12 at 73x; 4,096 holds 120 MB less.
const MAX_LAUNCH_PAIRS: usize = 4_096;

fn context() -> Option<&'static Arc<GpuContext>> {
    static CONTEXT: OnceLock<Option<Arc<GpuContext>>> = OnceLock::new();
    CONTEXT
        .get_or_init(|| match GpuContext::new() {
            Ok(context) => {
                info!(adapter = %context.adapter_info().name, "The pair-HMM is scored on the GPU");
                Some(context)
            }
            Err(error) => {
                warn!(%error, "No GPU for the pair-HMM; it is scored on the CPU");
                None
            }
        })
        .as_ref()
}

type Launched = (GpuPairs, Result<Vec<GpuScore>, GpuError>);

/// A worker's aligner on a thread of its own: `submit` a plan, `collect`
/// it back with its scores. The aligner's buffers grow to the largest
/// launch and stay, so it outlives the region's [`Batch`].
struct Lane {
    submit: mpsc::Sender<GpuPairs>,
    collect: mpsc::Receiver<Launched>,
}

thread_local! {
    static LANE: OnceCell<Option<Lane>> = const { OnceCell::new() };
}

impl Lane {
    fn spawn(context: &Arc<GpuContext>) -> std::io::Result<Self> {
        let (submit, launches) = mpsc::channel::<GpuPairs>();
        let (results, collect) = mpsc::channel();
        // Sleep between polls rather than spin a core while a launch runs.
        let options = KernelOptions { wait: Wait::Block, ..KernelOptions::default() };
        let context = Arc::clone(context);
        let aligner = move || GpuAligner::with_options(Arc::clone(&context), options);
        let name = thread::current().name().unwrap_or("worker").to_owned();
        thread::Builder::new().name(format!("{name}-gpu")).spawn(move || {
            let mut current = aligner();
            for pairs in launches {
                let scores = current.align(&pairs);
                // A failed launch can leave the aligner's buffers mapped;
                // a fresh one costs a regrowth and keeps the lane alive.
                if scores.is_err() {
                    current = aligner();
                }
                if results.send((pairs, scores)).is_err() {
                    break;
                }
            }
        })?;
        Ok(Self { submit, collect })
    }

    fn with<R>(f: impl FnOnce(Option<&Lane>) -> R) -> R {
        LANE.with(|lane| {
            f(lane
                .get_or_init(|| {
                    Lane::spawn(context()?)
                        .inspect_err(|error| {
                            warn!(%error, "No GPU thread for the pair-HMM; it is scored on the CPU");
                        })
                        .ok()
                })
                .as_ref())
        })
    }
}

/// The region's batch, created on first use; `None` when the pair-HMM
/// is scored on the CPU.
pub(super) fn batch(slot: &mut Option<Batch>, gpu: bool) -> Option<&mut Batch> {
    if gpu && slot.is_none() && Lane::with(|lane| lane.is_some()) {
        *slot = Some(Batch::default());
    }
    slot.as_mut()
}

/// A queued cluster: its columns, its prepared state, and where its
/// scores land in the launch.
struct Pending {
    columns: Vec<usize>,
    pass: ClusterPass,
    prepared: Prepared,
    range: Range<usize>,
}

#[derive(Default)]
pub(in super::super) struct Batch {
    filling: GpuPairs,
    pending: Vec<Pending>,
    /// The clusters of the launch on the lane, if one is.
    in_flight: Option<Vec<Pending>>,
}

static LAUNCHES: AtomicU64 = AtomicU64::new(0);
static PAIRS: AtomicU64 = AtomicU64::new(0);
static RESCORED: AtomicU64 = AtomicU64::new(0);

impl Batch {
    /// Queues a cluster, launching the batch once it is full; hands the
    /// cluster back when it cannot be queued, to score on the CPU.
    pub(super) fn push(
        &mut self,
        columns: &[usize],
        pass: ClusterPass,
        prepared: Prepared,
        pileups: &mut [PileupMetrics],
        scratch: &mut HaplotypeScratch,
    ) -> Option<Prepared> {
        match prepared.push(&mut self.filling) {
            Ok(range) => {
                self.pending.push(Pending { columns: columns.to_vec(), pass, prepared, range });
                if self.filling.len() >= MAX_LAUNCH_PAIRS {
                    self.launch(pileups, scratch);
                }
                None
            }
            Err(error) => {
                warn!(%error, "A cluster did not fit the GPU launch; it is scored on the CPU");
                Some(prepared)
            }
        }
    }

    /// Scores everything queued, waiting for the GPU.
    pub(super) fn flush(&mut self, pileups: &mut [PileupMetrics], scratch: &mut HaplotypeScratch) {
        if !self.pending.is_empty() {
            self.launch(pileups, scratch);
        }
        if self.in_flight.is_some() {
            self.filling = self.collect(pileups, scratch);
        }
    }

    /// Submits the filling plan, after collecting the launch before it.
    fn launch(&mut self, pileups: &mut [PileupMetrics], scratch: &mut HaplotypeScratch) {
        let spare =
            if self.in_flight.is_some() { self.collect(pileups, scratch) } else { GpuPairs::new() };
        let plan = mem::replace(&mut self.filling, spare);
        let pending = mem::take(&mut self.pending);
        match Lane::with(|lane| lane.map(|lane| lane.submit.send(plan))) {
            Some(Ok(())) => self.in_flight = Some(pending),
            _ => {
                warn!("The pair-HMM's GPU thread is gone; its clusters are scored on the CPU");
                settle(pending, None, pileups, scratch);
            }
        }
    }

    /// Waits for the launch on the lane, assigns its scores, and returns
    /// its plan, cleared, for the next launch.
    fn collect(
        &mut self,
        pileups: &mut [PileupMetrics],
        scratch: &mut HaplotypeScratch,
    ) -> GpuPairs {
        let pending = self.in_flight.take().unwrap_or_default();
        let Some((mut plan, scores)) = Lane::with(|lane| lane?.collect.recv().ok()) else {
            warn!("The pair-HMM's GPU thread is gone; its clusters are scored on the CPU");
            settle(pending, None, pileups, scratch);
            return GpuPairs::new();
        };
        let scores = scores
            .inspect_err(|error| {
                warn!(%error, "The pair-HMM's GPU launch failed; its clusters are scored on the CPU");
            })
            .ok();
        if let Some(scores) = &scores {
            let rescored = scores.iter().filter(|score| score.trusted().is_none()).count();
            let total = |counter: &AtomicU64, n: u64| counter.fetch_add(n, Ordering::Relaxed) + n;
            debug!(
                launches = total(&LAUNCHES, 1),
                pairs = total(&PAIRS, scores.len() as u64),
                rescored = total(&RESCORED, rescored as u64),
                launch_pairs = scores.len(),
                "Pair-HMM GPU launch (running totals)"
            );
        }
        settle(pending, scores.as_deref(), pileups, scratch);
        plan.clear();
        plan
    }
}

/// A batch dropped with a launch in flight (a region that failed) takes
/// the launch's result off the lane, or the next region would collect it.
impl Drop for Batch {
    fn drop(&mut self) {
        if self.in_flight.is_some()
            && !Lane::with(|lane| lane.is_some_and(|lane| lane.collect.recv().is_ok()))
        {
            debug!("A pair-HMM launch was lost with its GPU thread");
        }
    }
}

/// Gives each cluster its scores from the launch, or scores it on the
/// CPU where the launch has none.
fn settle(
    pending: Vec<Pending>,
    scores: Option<&[GpuScore]>,
    pileups: &mut [PileupMetrics],
    scratch: &mut HaplotypeScratch,
) {
    for Pending { columns, pass, prepared, range } in pending {
        let result = match scores.and_then(|scores| scores.get(range)) {
            Some(scores) => prepared.finish_gpu(scores),
            None => prepared.align(scratch),
        };
        assign_scores(pileups, &columns, pass, result);
    }
}
