//! Panics: who reports them, and which ones fail a run.
//!
//! A panic prints color-eyre's crash report from the panic hook, while the
//! panicking thread's spans are still entered, so the report and its issue
//! link show the region it happened in. Whoever handles it afterwards only
//! has to stop the run: ordair returns a panic in a worker or in the writer
//! as an error, and [`PanicCheck`] fails a run in which any other thread
//! panicked. The one exception is a thread started by [`spawn_recovering`],
//! which the run can do without: its panic is logged, and a recovery runs.

use crate::utils::logging::{BUG_MESSAGE, ThisIsABug as _};
use color_eyre::eyre::{Context as _, Result, ensure, eyre};
use std::{
    cell::Cell,
    panic::{self, AssertUnwindSafe},
    sync::{
        Once,
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

fn register_panic_hook() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let previous = panic::take_hook();
        panic::set_hook(Box::new(move |info| {
            if RECOVERS.with(Cell::get) {
                error!(panic = %info, note = BUG_MESSAGE, "A thread panicked and recovered");
            } else {
                PANICS.fetch_add(1, Ordering::Relaxed);
                previous(info);
            }
        }));
    });
}

/// Fails a run in which a thread panicked, e.g. one started by a library (for
/// compression) whose panic nobody gets as an error. Its output would
/// otherwise be committed as complete.
#[must_use = "only `finish` reports the panic"]
pub struct PanicCheck {
    panics_before: usize,
}

impl PanicCheck {
    pub fn start() -> Self {
        register_panic_hook();
        Self { panics_before: PANICS.load(Ordering::Relaxed) }
    }

    pub fn finish(self) -> Result<()> {
        ensure!(
            PANICS.load(Ordering::Relaxed) == self.panics_before,
            "A thread panicked (see above), output files may be incomplete"
        );
        Ok(())
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
