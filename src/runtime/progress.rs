use super::threads::Thread;
use anstyle_progress::TermProgress;
use jiff::Timestamp;
use std::{
    io::{IsTerminal as _, Write as _},
    ops::Add,
    sync::{
        Arc, Once,
        atomic::{AtomicBool, AtomicU8, Ordering},
        mpsc,
    },
    time::{Duration, Instant},
};
use tracing::{info, warn};

/// Flag set by the SIGINFO (ctrl+t on macOS) / SIGUSR1 (Linux) signal handler
/// to request a progress report from the writer thread.
static PRINT_REQUESTED: AtomicBool = AtomicBool::new(false);

extern "C" fn on_progress_signal(_: libc::c_int) {
    PRINT_REQUESTED.store(true, Ordering::Relaxed);
}

/// Register the OS signal handler that sets [`PRINT_REQUESTED`] on
/// SIGINFO (macOS ctrl+t) or SIGUSR1 (Linux).
///
/// Safe to call multiple times; registration happens only once per process.
pub fn register_signal_handler() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        #[cfg(target_os = "macos")]
        let sig = libc::SIGINFO;
        #[cfg(target_os = "linux")]
        let sig = libc::SIGUSR1;

        // SAFETY: We only write to an atomic bool in the handler — async-signal-safe.
        #[cfg(any(target_os = "macos", target_os = "linux"))]
        unsafe {
            let mut sa: libc::sigaction = std::mem::zeroed();
            sa.sa_sigaction = on_progress_signal as *const () as libc::sighandler_t;
            libc::sigemptyset(&mut sa.sa_mask);
            sa.sa_flags = libc::SA_RESTART;
            if libc::sigaction(sig, &sa, std::ptr::null_mut()) != 0 {
                tracing::warn!(error = ?std::io::Error::last_os_error(), "Failed to register the progress signal handler");
            }
        }
    });
}

// Require enough completed segments and elapsed time before the first
// automatic ETA print, so the estimate is based on meaningful data.
const MIN_CALIBRATION_SEGMENTS: usize = 50;
const MIN_CALIBRATION_TIME: Duration = Duration::from_secs(30);

/// Tracks segment completion in the writer thread and logs ETA estimates.
pub struct ProgressTracker {
    enabled: bool,
    total: usize,
    completed: usize,
    calibrated: bool,
    start: Instant,
    terminal: Option<TerminalProgress>,
}

impl ProgressTracker {
    /// Create a progress tracker, disabled if the "CI" environment variable is set (to avoid spamming CI logs).
    pub fn new(total_segments: usize) -> Self {
        let enabled = std::env::var("CI").err() == Some(std::env::VarError::NotPresent);
        let terminal = TerminalProgress::detect(enabled);
        Self {
            total: total_segments,
            completed: 0,
            calibrated: false,
            start: Instant::now(),
            enabled,
            terminal,
        }
    }

    /// Call after each segment has been fully written
    pub fn segment_done(&mut self) {
        if !self.enabled {
            return;
        }

        self.completed += 1;

        if let Some(terminal) = &self.terminal {
            terminal.update(percent_complete(self.completed, self.total));
        }

        let signal_requested = PRINT_REQUESTED
            .compare_exchange(true, false, Ordering::Relaxed, Ordering::Relaxed)
            .is_ok();

        if signal_requested {
            self.log();
        } else if !self.calibrated
            && self.completed >= MIN_CALIBRATION_SEGMENTS
            && self.start.elapsed() >= MIN_CALIBRATION_TIME
        {
            self.calibrated = true;
            self.log_one_off_estimate();
        }
    }

    fn log(&self) {
        let Estimate { percent, eta, done } = self.estimate();

        info!(
            percent = %format!("{percent:.1}%"),
            time_left = %format_duration(eta),
            done_at = %done,
            "{}/{} segments",
            self.completed,
            self.total,
        );
    }

    fn log_one_off_estimate(&self) {
        let Estimate { eta, done, .. } = self.estimate();

        info!(
            time_left = %format_duration(eta),
            done_at = %done,
            "Runtime estimate",
        );
    }

    fn estimate(&self) -> Estimate {
        let elapsed = self.start.elapsed();
        let pct = self.completed as f64 / self.total as f64 * 100.0;
        let remaining = self.total - self.completed;
        let secs_per_segment = elapsed.as_secs_f64() / self.completed as f64;
        let eta = Duration::from_secs_f64(secs_per_segment * remaining as f64);
        let done = Timestamp::now().add(eta);

        Estimate { percent: pct, eta, done }
    }
}

/// Whole percent, clamped to `0..=100` as required by the OSC 9;4 protocol.
fn percent_complete(completed: usize, total: usize) -> u8 {
    let percent = completed.saturating_mul(100).checked_div(total).unwrap_or(100);
    u8::try_from(percent.min(100)).unwrap_or(100)
}

/// Progress indicator in the terminal tab/taskbar via the `ConEmu` OSC 9;4 escape
/// sequence (e.g. Windows Terminal, Ghostty, `WezTerm`, iTerm2, Konsole).
///
/// Only used when stderr is a terminal known to support it: other terminals
/// interpret bare OSC 9 as a desktop notification, and in batch jobs (Slurm,
/// Nextflow) the sequence would only end up as noise in log files.
///
/// The state is re-sent from a background thread every [`KEEPALIVE_INTERVAL`]
/// because some terminals (Ghostty) reset it after ~15 s without an update, and
/// segments can take longer than that to arrive at the in-order writer. This also
/// means a killed process (e.g. ctrl+c, which skips `Drop`) doesn't leave a stale
/// indicator behind in those terminals.
struct TerminalProgress {
    percent: Arc<AtomicU8>,
    stop: Option<mpsc::Sender<()>>,
    keepalive: Option<Thread<bool>>,
}

const KEEPALIVE_INTERVAL: Duration = Duration::from_secs(1);

impl TerminalProgress {
    fn detect(enabled: bool) -> Option<Self> {
        if !enabled || !anstyle_progress::supports_term_progress(std::io::stderr().is_terminal()) {
            return None;
        }

        let percent = Arc::new(AtomicU8::new(0));
        let (stop, stop_requested) = mpsc::channel::<()>();
        // The indicator is cosmetic, so losing it must not stop the run. The
        // thread's result says whether stderr is still writable, which a panic
        // doesn't rule out, so still try to remove the indicator afterwards.
        let keepalive = super::threads::spawn_recovering(
            "term-progress",
            (Arc::clone(&percent), stop_requested),
            |(percent, stop_requested)| loop {
                let current = percent.load(Ordering::Relaxed);
                if !write_to_stderr(TermProgress::start().percent(current)) {
                    return false;
                }
                match stop_requested.recv_timeout(KEEPALIVE_INTERVAL) {
                    Err(mpsc::RecvTimeoutError::Timeout) => {}
                    Ok(()) | Err(mpsc::RecvTimeoutError::Disconnected) => return true,
                }
            },
            |_| true,
        );

        match keepalive {
            Ok(keepalive) => Some(Self { percent, stop: Some(stop), keepalive: Some(keepalive) }),
            Err(error) => {
                warn!(error = format!("{error:#}"), "Terminal progress indicator disabled");
                None
            }
        }
    }

    fn update(&self, percent: u8) {
        self.percent.store(percent, Ordering::Relaxed);
    }
}

impl Drop for TerminalProgress {
    fn drop(&mut self) {
        // Dropping the sender wakes the keepalive thread immediately
        drop(self.stop.take());
        let still_writable = match self.keepalive.take().map(Thread::join) {
            Some(Ok(writable)) => writable,
            Some(Err(error)) => {
                warn!(error = format!("{error:#}"), "Terminal progress indicator failed");
                false
            }
            None => false,
        };
        if still_writable {
            write_to_stderr(TermProgress::remove());
        }
    }
}

fn write_to_stderr(progress: TermProgress) -> bool {
    let mut stderr = std::io::stderr().lock();
    match write!(stderr, "{progress}").and_then(|()| stderr.flush()) {
        Ok(()) => true,
        Err(error) => {
            warn!(?error, "Failed to write terminal progress indicator, disabling it");
            false
        }
    }
}

struct Estimate {
    percent: f64,
    eta: Duration,
    done: Timestamp,
}

fn format_duration(d: Duration) -> String {
    let total_secs = d.as_secs();
    let hours = total_secs / 3600;
    let mins = (total_secs % 3600) / 60;
    let secs = total_secs % 60;
    if hours > 0 {
        format!("{hours}h {mins:02}m {secs:02}s")
    } else if mins > 0 {
        format!("{mins}m {secs:02}s")
    } else {
        format!("{secs}s")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percent_is_clamped_and_handles_empty_totals() {
        assert_eq!(percent_complete(0, 200), 0);
        assert_eq!(percent_complete(1, 200), 0);
        assert_eq!(percent_complete(2, 200), 1);
        assert_eq!(percent_complete(200, 200), 100);
        assert_eq!(percent_complete(201, 200), 100);
        assert_eq!(percent_complete(0, 0), 100);
        assert_eq!(percent_complete(usize::MAX, 3), 100);
    }

    #[test]
    fn osc_sequences() {
        assert_eq!(TermProgress::start().percent(42).to_string(), "\x1b]9;4;1;42\x1b\\");
        assert_eq!(TermProgress::remove().to_string(), "\x1b]9;4;0;\x1b\\");
    }
}
