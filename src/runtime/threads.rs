//! Threads whose panics become errors instead of crash reports.
//!
//! [`catch_panic`] runs work on threads we don't create (rayon workers), and
//! [`spawn_recovering`] starts a named thread the run can do without. A panic
//! in either is caught on that thread and turned into a [`Panicked`] error
//! with the thread name, message and location. The panic hook only records
//! these details instead of printing color-eyre's crash report, so each crash
//! is reported once, by whoever handles it:
//! - [`catch_panic`] returns the error, which stops the run like any other.
//! - [`spawn_recovering`] logs it, runs the given recovery on the same thread
//!   (the GPU falls back to the CPU), and the run continues.
//!
//! Panics anywhere else (e.g. the main thread) get color-eyre's usual crash
//! report, and fail a running pipeline at its end (see [`PanicCheck`]).

use crate::utils::logging::{BUG_MESSAGE, ThisIsABug as _, backtrace, bug};
use color_eyre::eyre::{Context as _, Report, Result, ensure, eyre};
use std::{
    cell::{Cell, RefCell},
    panic::{self, AssertUnwindSafe, PanicHookInfo},
    sync::{
        Once,
        atomic::{AtomicUsize, Ordering},
    },
    thread::{self, JoinHandle},
};
use tracing::{debug, error};

/// A panic caught on a thread, as an error.
#[derive(Debug, Clone, thiserror::Error)]
#[error("Thread `{thread}` panicked at {location}: {message}")]
pub struct Panicked {
    thread: String,
    message: String,
    location: String,
    /// Where the panic happened, see [`backtrace`]. Only the panic hook can
    /// capture this: by the time the panic is caught, the stack has unwound.
    backtrace: Option<String>,
}

impl Panicked {
    fn from_hook(info: &PanicHookInfo<'_>) -> Self {
        Self {
            thread: current_thread_name(),
            message: info.payload_as_str().unwrap_or("unknown panic payload").to_string(),
            location: info
                .location()
                .map_or_else(|| "unknown location".into(), ToString::to_string),
            backtrace: backtrace(),
        }
    }

    fn into_report(mut self) -> Report {
        let backtrace = self.backtrace.take();
        bug(self, backtrace)
    }
}

fn current_thread_name() -> String {
    thread::current().name().unwrap_or("unnamed").to_string()
}

thread_local! {
    /// Set while running inside [`catch`].
    static CATCHING: Cell<bool> = const { Cell::new(false) };
    /// Recorded by the panic hook, picked up by [`catch`] after unwinding.
    static CAUGHT: RefCell<Option<Report>> = const { RefCell::new(None) };
}

/// Panics that no [`catch`] saw, see [`PanicCheck`].
static UNCAUGHT_PANICS: AtomicUsize = AtomicUsize::new(0);

/// The hook runs on the panicking thread before unwinding, so this is where
/// the location and the spans (captured when the report is created, for the
/// issue link) are known.
fn register_panic_hook() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let previous = panic::take_hook();
        panic::set_hook(Box::new(move |info| {
            if CATCHING.with(Cell::get) {
                let report = Panicked::from_hook(info).into_report();
                CAUGHT.with(|caught| caught.replace(Some(report)));
            } else {
                UNCAUGHT_PANICS.fetch_add(1, Ordering::Relaxed);
                previous(info);
            }
        }));
    });
}

/// Fails a run in which a thread panicked without anyone getting an error for
/// it: a thread started by a library (e.g. for compression) rather than
/// through this module. Its output would otherwise be committed as complete.
#[must_use = "only `finish` reports the panic"]
pub struct PanicCheck {
    uncaught_before: usize,
}

impl PanicCheck {
    pub fn start() -> Self {
        register_panic_hook();
        Self { uncaught_before: UNCAUGHT_PANICS.load(Ordering::Relaxed) }
    }

    pub fn finish(self) -> Result<()> {
        ensure!(
            UNCAUGHT_PANICS.load(Ordering::Relaxed) == self.uncaught_before,
            "A thread panicked (see above), output files may be incomplete"
        );
        Ok(())
    }
}

/// `AssertUnwindSafe` is fine for all callers: after a panic, whatever `work`
/// captured is only used to clean up, never to continue the interrupted work.
fn catch<T>(work: impl FnOnce() -> T) -> Result<T> {
    register_panic_hook();
    let outer = CATCHING.with(|catching| catching.replace(true));
    let result = panic::catch_unwind(AssertUnwindSafe(work));
    CATCHING.with(|catching| catching.set(outer));
    result.map_err(|_| {
        // Only empty if something replaced our hook
        CAUGHT
            .with(RefCell::take)
            .unwrap_or_else(|| bug(eyre!("Thread `{}` panicked", current_thread_name()), None))
    })
}

/// Run `work` on the current thread. A panic is returned as an error.
///
/// For work on threads we don't spawn ourselves, i.e. rayon workers: catching
/// per segment on the worker (rather than around `install`) keeps the error
/// next to the segment it happened in.
pub fn catch_panic<T>(work: impl FnOnce() -> Result<T>) -> Result<T> {
    catch(work).this_is_a_bug()?
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
    let handle = thread::Builder::new()
        .name(name.to_string())
        .spawn(move || {
            let mut state = state;
            match catch(|| work(&mut state)) {
                Ok(value) => value,
                Err(report) => {
                    error!(
                        panic = format!("{report:#}"),
                        note = BUG_MESSAGE,
                        "A thread panicked and recovered"
                    );
                    debug!("{report:?}");
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
    fn catch_panic_turns_panics_into_errors() {
        let caught = catch_panic::<()>(|| panic!("boom in {}", "worker"));
        let message = format!("{:#}", caught.expect_err("panic should become an error"));
        assert!(message.contains("boom in worker"), "{message}");
        assert!(
            message.contains("src/runtime/threads.rs"),
            "location should be included: {message}"
        );

        assert_eq!(catch_panic(|| Ok(7)).expect("no panic"), 7);
    }

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
