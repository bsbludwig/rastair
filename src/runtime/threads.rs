//! Panics: who reports them, and which ones fail a run.
//!
//! A panic's crash report is made in the panic hook, while the panicking
//! thread's spans are still entered, so the report and its issue link show the
//! region it happened in. During a run, it is printed when the run ends, after
//! the warnings about the outputs it left behind. Whoever handles the panic
//! only has to stop the run: ordair returns a panic in a worker or in the writer
//! as an error, and [`PanicCheck`] fails a run in which any other thread
//! panicked. The one exception is a thread started by [`spawn_recovering`],
//! which the run can do without: its panic is logged, and a recovery runs.

use crate::utils::logging::{BUG_MESSAGE, ThisIsABug as _, panic_report};
use color_eyre::eyre::{Context as _, Result, ensure, eyre};
use std::{
    cell::Cell,
    panic::{self, AssertUnwindSafe},
    sync::{
        Mutex, Once, PoisonError,
        atomic::{AtomicUsize, Ordering},
    },
    thread::{self, JoinHandle},
};
use tracing::error;

thread_local! {
    /// Set on a thread started by [`spawn_recovering`].
    static RECOVERS: Cell<bool> = const { Cell::new(false) };
}

/// Panics on threads that don't recover, see [`PanicCheck`].
static PANICS: AtomicUsize = AtomicUsize::new(0);

/// Crash reports held back while a [`PanicCheck`] is alive; `None` otherwise.
static DEFERRED: Mutex<Option<Vec<String>>> = Mutex::new(None);

fn deferred() -> std::sync::MutexGuard<'static, Option<Vec<String>>> {
    DEFERRED.lock().unwrap_or_else(PoisonError::into_inner)
}

fn register_panic_hook() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let previous = panic::take_hook();
        panic::set_hook(Box::new(move |info| {
            if RECOVERS.with(Cell::get) {
                error!(
                    panic = info.payload_as_str(),
                    location = info.location().map(tracing::field::display),
                    note = BUG_MESSAGE,
                    "A thread panicked and recovered"
                );
            } else {
                PANICS.fetch_add(1, Ordering::Relaxed);
                if deferred().is_none() {
                    return previous(info);
                }
                // Rendered outside the lock: a run that ends meanwhile gets the
                // report printed right away instead.
                let report = panic_report(info);
                match deferred().as_mut() {
                    Some(reports) => reports.push(report),
                    None => return eprintln!("{report}"),
                }
                // The report waits for the run to end, which can take a while
                // when a thread nobody joins panicked and the run carries on.
                error!(
                    panic = info.payload_as_str(),
                    location = info.location().map(tracing::field::display),
                    "A thread panicked, its crash report follows when the run ends"
                );
            }
        }));
    });
}

/// Fails a run in which a thread panicked, e.g. one started by a library (for
/// compression) whose panic nobody gets as an error. Its output would
/// otherwise be committed as complete.
///
/// While it is alive, crash reports are held back, and printed when it is
/// dropped. Create it before the outputs, so it is dropped after them.
#[must_use = "only `finish` reports the panic"]
pub struct PanicCheck {
    panics_before: usize,
}

impl PanicCheck {
    pub fn start() -> Self {
        register_panic_hook();
        deferred().get_or_insert_with(Vec::new);
        Self { panics_before: PANICS.load(Ordering::Relaxed) }
    }

    pub fn finish(&self) -> Result<()> {
        ensure!(
            PANICS.load(Ordering::Relaxed) == self.panics_before,
            "A thread panicked (see above), output files may be incomplete"
        );
        Ok(())
    }
}

impl Drop for PanicCheck {
    fn drop(&mut self) {
        // Blank lines set the report apart from the logs and the error around it
        for report in deferred().take().into_iter().flatten() {
            eprintln!("\n{report}\n");
        }
    }
}

/// A thread started by [`spawn_recovering`].
#[must_use = "join the thread to get its result"]
pub struct Thread<T> {
    name: String,
    handle: JoinHandle<T>,
}

impl<T> Thread<T> {
    /// The thread's result, or an error if its recovery panicked too.
    pub fn join(self) -> Result<T> {
        self.handle
            .join()
            .map_err(|_| eyre!("Thread `{}` panicked while recovering from a panic", self.name))
            .this_is_a_bug()
    }
}

/// Start a thread the run can do without. A panic in `work` is logged and does
/// not stop the run; `recover` then runs on the same thread and provides the
/// thread's result.
///
/// `recover` gets back the `state` that `work` operated on, possibly left
/// half-updated by the panic: use it to clean up, e.g. to release whoever
/// waits on this thread.
pub fn spawn_recovering<S, T>(
    name: &str,
    state: S,
    work: impl FnOnce(&mut S) -> T + Send + 'static,
    recover: impl FnOnce(S) -> T + Send + 'static,
) -> Result<Thread<T>>
where
    S: Send + 'static,
    T: Send + 'static,
{
    register_panic_hook();
    let handle = thread::Builder::new()
        .name(name.to_string())
        .spawn(move || {
            RECOVERS.with(|recovers| recovers.set(true));
            let mut state = state;
            // `AssertUnwindSafe`: after a panic, `recover` only uses `state`
            // to clean up, never to continue the interrupted work.
            match panic::catch_unwind(AssertUnwindSafe(|| work(&mut state))) {
                Ok(value) => value,
                Err(_) => {
                    RECOVERS.with(|recovers| recovers.set(false));
                    recover(state)
                }
            }
        })
        .wrap_err_with(|| format!("Failed to start thread `{name}`"))?;
    Ok(Thread { name: name.to_string(), handle })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovering_thread_gets_state_back() {
        let thread = spawn_recovering(
            "test-recovers",
            vec![1],
            |state| -> usize {
                state.push(2);
                panic!("boom after partial update")
            },
            |state| state.len(),
        )
        .expect("spawn");
        assert_eq!(thread.join().expect("recovered"), 2);
    }
}
