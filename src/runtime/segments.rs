//! Handing segments from the workers to the writer, in order.

use crate::utils::logging::Bug;
use color_eyre::eyre::{Report, Result, WrapErr as _, eyre};
use std::{
    fmt,
    sync::{
        Mutex, PoisonError,
        atomic::{AtomicUsize, Ordering},
    },
};
use tracing::{error, warn};

/// Segments that failed to process.
///
/// Their output is missing, but the run goes on with the other segments rather
/// than losing all of them to one bad region. It still fails in the end (see
/// [`Self::check`]), so the output is not taken for complete.
#[derive(Debug, Default)]
pub struct FailedSegments {
    count: AtomicUsize,
    /// Returned from [`Self::check`], so that a bug gets its issue link.
    first_bug: Mutex<Option<Report>>,
}

impl FailedSegments {
    /// Log `error` and count the segment as failed.
    pub fn record(&self, segment: &impl fmt::Display, error: Report) {
        error!(
            %segment,
            error = format!("{error:#}"),
            "Failed to process segment, its output is missing"
        );
        self.count.fetch_add(1, Ordering::Relaxed);
        if Bug::marks(error.as_ref()) {
            self.first_bug.lock().unwrap_or_else(PoisonError::into_inner).get_or_insert(error);
        }
    }

    pub fn check(&self, total: usize) -> Result<()> {
        let failed = match self.count.load(Ordering::Relaxed) {
            0 => return Ok(()),
            failed => SegmentsFailed { failed, total },
        };
        match self.first_bug.lock().unwrap_or_else(PoisonError::into_inner).take() {
            Some(bug) => Err(bug.wrap_err(failed)),
            None => Err(failed.into()),
        }
    }
}

#[derive(Debug, thiserror::Error)]
#[error(
    "{failed} of {total} segments failed to process (see the errors above), \
     output files are missing their records"
)]
pub struct SegmentsFailed {
    failed: usize,
    total: usize,
}

/// Of two steps that both ran, e.g. writing and then closing an output, the
/// first error is the one to report; a second one is only logged.
pub fn first_error(first: Result<()>, then: Result<()>) -> Result<()> {
    match (first, then) {
        (Err(error), Err(then)) => {
            warn!(error = format!("{then:#}"), "Another error followed");
            Err(error)
        }
        (first, then) => first.and(then),
    }
}

/// Runs `work` on each item on `pool`'s threads and hands the results to
/// `write` on the calling thread, in order (see [`ordair::in_order`]).
///
/// The workers pull `items` one at a time, so an iterator that reads its input
/// runs on them, serially. An error stops the run: no new item is started, so
/// the output is a contiguous prefix. To skip a segment that failed in `work`
/// instead, count it in [`FailedSegments`] and return no records. A panic in
/// `init`, `work` or `write` stops the run like an error, after the panic hook
/// printed its report.
pub fn process_in_order<I, S, R: Send>(
    pool: &rayon::ThreadPool,
    items: I,
    init: impl Fn() -> Result<S> + Sync,
    work: impl Fn(&mut S, I::Item) -> Result<R> + Sync,
    mut write: impl FnMut(R) -> Result<()>,
) -> Result<()>
where
    I: IntoIterator<IntoIter: Send, Item: Send + fmt::Display>,
{
    let mut written = 0_usize;
    // ordair's default window, four items per worker, lets the others carry
    // on past one slow item (deep coverage, a repeat) without every item
    // behind it piling up in memory.
    let result = ordair::in_order(items)
        .pool(pool)
        .map_init(init, |state, item| {
            let name = item.to_string();
            work(state, item).wrap_err_with(|| format!("Failed to process `{name}`"))
        })
        .try_for_each(|processed| {
            write(processed?).wrap_err("Failed to write the output")?;
            written += 1;
            Ok(())
        })
        .map_err(|panic| match panic.origin() {
            // The panic hook printed the message already
            ordair::PanicOrigin::Worker => eyre!("A worker panicked (see the crash report above)"),
            ordair::PanicOrigin::Consumer => {
                eyre!("Writing the output panicked (see the crash report above)")
            }
            _ => Report::new(panic),
        })
        .flatten();

    if result.is_err() {
        warn!(written, "Output is incomplete: stopped before all items were written");
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_failed_segments() {
        let failed = FailedSegments::default();
        assert!(failed.check(3).is_ok());
        failed.record(&"chr1:1-10", color_eyre::eyre::eyre!("boom"));
        let error = failed.check(3).expect_err("one segment failed");
        assert!(error.to_string().starts_with("1 of 3 segments failed"), "{error}");
        assert!(!Bug::marks(error.as_ref()));
    }

    #[test]
    fn a_failed_segment_keeps_the_bug() {
        use crate::utils::logging::ThisIsABug as _;
        let failed = FailedSegments::default();
        failed.record(&"chr1:1-10", color_eyre::eyre::eyre!("user error"));
        let bug = Err::<(), _>(color_eyre::eyre::eyre!("bug")).this_is_a_bug();
        failed.record(&"chr1:11-20", bug.expect_err("constructed as Err"));
        let error = failed.check(3).expect_err("two segments failed");
        assert!(Bug::marks(error.as_ref()));
        assert_eq!(
            format!("{error:#}"),
            "2 of 3 segments failed to process (see the errors above), output files are missing their records: bug"
        );
    }
}
