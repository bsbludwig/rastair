use jiff::Timestamp;
use std::{
    ops::Add,
    sync::{
        Once,
        atomic::{AtomicBool, Ordering},
    },
    time::{Duration, Instant},
};
use tracing::info;

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
        // Clear any stale flag left over from a previous call() invocation.
        PRINT_REQUESTED.store(false, Ordering::Relaxed);

        #[cfg(target_os = "macos")]
        let sig = libc::SIGINFO;
        #[cfg(target_os = "linux")]
        let sig = libc::SIGUSR1;

        // SAFETY: We only write to an atomic bool in the handler — async-signal-safe.
        #[cfg(any(target_os = "macos", target_os = "linux"))]
        unsafe {
            let mut sa: libc::sigaction = std::mem::zeroed();
            sa.sa_sigaction = on_progress_signal as libc::sighandler_t;
            libc::sigemptyset(&mut sa.sa_mask);
            sa.sa_flags = libc::SA_RESTART;
            libc::sigaction(sig, &sa, std::ptr::null_mut());
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
}

impl ProgressTracker {
    /// Create a progress tracker, disabled if the "CI" environment variable is set (to avoid spamming CI logs).
    pub fn new(total_segments: usize) -> Self {
        let enabled = std::env::var("CI").err() == Some(std::env::VarError::NotPresent);
        Self {
            total: total_segments,
            completed: 0,
            calibrated: false,
            start: Instant::now(),
            enabled,
        }
    }

    /// Call after each segment has been fully written
    pub fn segment_done(&mut self) {
        if !self.enabled {
            return;
        }

        self.completed += 1;

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
