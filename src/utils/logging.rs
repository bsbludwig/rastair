use color_eyre::{ErrorKind, Section as _, SectionExt as _, eyre::Report};
use std::{
    error::Error as StdError,
    fmt,
    panic::PanicHookInfo,
    sync::{
        OnceLock,
        atomic::{AtomicBool, Ordering},
    },
};
use tracing_subscriber::{EnvFilter, layer::SubscriberExt as _};

pub static LOG_VAR: &str = "RASTAIR_LOG";
pub static BUG_MESSAGE: &str = "This is a bug in Rastair.";
const ISSUE_URL: &str = "https://github.com/bsbludwig/rastair/issues/new";

/// Setup logging and error handling
///
/// To be called once at the start of the program.
pub fn setup_logging(verbose: bool) {
    setup_tracing(verbose);
    setup_eyre(verbose);
}

fn setup_tracing(verbose: bool) {
    let subscriber = {
        let default_log_settings = if verbose { "info,rastair=debug" } else { "warn,rastair=info" };
        let mut env_filter = EnvFilter::new(default_log_settings);
        if let Ok(env) = std::env::var(LOG_VAR) {
            for directive in env.split(',') {
                if directive.is_empty() {
                    continue;
                }
                match directive.parse() {
                    Ok(parsed_directive) => {
                        env_filter = env_filter.add_directive(parsed_directive);
                    }
                    Err(error) => {
                        eprintln!("Warning: Invalid log directive `{directive}`: {error:#}");
                    }
                }
            }
        }

        tracing_subscriber::Registry::default()
            .with(tracing_error::ErrorLayer::default())
            .with(env_filter)
            .with(
                tracing_subscriber::fmt::Layer::default()
                    .with_target(true)
                    .with_thread_names(verbose)
                    // .with_span_events(FmtSpan::CLOSE) // maybe enable with flag
                    .with_writer(std::io::stderr),
            )
    };
    if let Err(error) = tracing::subscriber::set_global_default(subscriber) {
        eprintln!("Failed to register logging: {error:#}");
    }
}

fn setup_eyre(verbose: bool) {
    // color-eyre puts a backtrace into the issue link, which is then ~30 kB,
    // far over the ~8 kB GitHub accepts. So color-eyre must capture none, and
    // only `backtrace()` decides, from what the user asked for at startup.
    let wants_backtrace =
        verbose || std::env::var_os("RUST_BACKTRACE").is_some_and(|value| value != "0");
    WANTS_BACKTRACE.store(wants_backtrace, Ordering::Relaxed);
    // SAFETY: This is set at the very start of the program
    unsafe { std::env::set_var("RUST_BACKTRACE", "0") };
    if std::env::var_os("RUST_LIB_BACKTRACE").is_none() {
        // SAFETY: as above
        unsafe { std::env::set_var("RUST_LIB_BACKTRACE", "0") };
    }
    let hooks = hook_builder(verbose).try_into_hooks().and_then(|(panic_hook, eyre_hook)| {
        eyre_hook.install()?;
        PANIC_HOOK
            .set(panic_hook)
            .map_err(|_| color_eyre::eyre::eyre!("The panic hook was already set up"))?;
        std::panic::set_hook(Box::new(|info| eprintln!("{}", panic_report(info))));
        Ok(())
    });
    if let Err(error) = hooks.note("Seeing this error message is somewhat ironic, we know") {
        eprintln!("Failed to register panic handler: {error:#}");
    }
}

static PANIC_HOOK: OnceLock<color_eyre::config::PanicHook> = OnceLock::new();

/// The crash report for a panic, with a backtrace if one was asked for.
///
/// Call it from the panic hook: the report's span trace shows the spans
/// entered on the current thread.
pub fn panic_report(info: &PanicHookInfo<'_>) -> String {
    let report = match PANIC_HOOK.get() {
        Some(hook) => hook.panic_report(info).to_string(),
        None => info.to_string(),
    };
    match backtrace() {
        Some(backtrace) => format!("{report}\n\nBacktrace:\n{backtrace}"),
        None => report,
    }
}

fn hook_builder(verbose: bool) -> color_eyre::config::HookBuilder {
    color_eyre::config::HookBuilder::default()
        .panic_section(BUG_MESSAGE)
        .issue_url(ISSUE_URL)
        .add_issue_metadata("version", env!("CARGO_PKG_VERSION"))
        .add_issue_metadata(
            "pileup backend",
            if cfg!(feature = "experimental-seqair") { "seqair" } else { "htslib" },
        )
        // Only defects get a link: a user error like "BAM not found" is not
        // worth an issue, and its message is full of the user's paths.
        .issue_filter(|kind| match kind {
            ErrorKind::NonRecoverable(_) => true,
            ErrorKind::Recoverable(error) => Bug::marks(error),
        })
        .capture_span_trace_by_default(true)
        .display_env_section(verbose)
        .display_location_section(verbose)
        .theme(if std::env::var("NO_COLOR").is_ok() {
            color_eyre::config::Theme::new()
        } else {
            color_eyre::config::Theme::dark()
        })
}

pub trait ThisIsABug<T> {
    /// Note to user that this is a bug in the program not an expected error,
    /// and offer a pre-filled issue link for it.
    ///
    /// The link shows the spans entered when the report is created here, and
    /// sections added before this call are lost, so call it early.
    fn this_is_a_bug(self) -> Result<T, Report>;
}

impl<T, E> ThisIsABug<T> for Result<T, E>
where
    E: Into<Report>,
{
    // Not `map_err`: a closure would end the `#[track_caller]` chain, and the
    // report's location would be in here instead of at the caller.
    #[track_caller]
    fn this_is_a_bug(self) -> Result<T, Report> {
        match self {
            Ok(value) => Ok(value),
            Err(error) => Err(bug(error, backtrace())),
        }
    }
}

/// See [`ThisIsABug`]. The backtrace is shown in the report but, unlike
/// color-eyre's own, not put into the issue link.
#[track_caller]
fn bug(error: impl Into<Report>, backtrace: Option<String>) -> Report {
    let report = error.into();
    if Bug::marks(report.as_ref()) {
        // Wrapping it again would lose its spans and sections
        return report;
    }
    let report = Report::new(Bug(report.into())).note(BUG_MESSAGE);
    match backtrace {
        Some(backtrace) => report.section(backtrace.header("Backtrace:")),
        None => report,
    }
}

/// Whether `RUST_BACKTRACE` (or `--verbose`) asked for a backtrace before
/// `setup_eyre` turned the variable off.
static WANTS_BACKTRACE: AtomicBool = AtomicBool::new(false);

fn backtrace() -> Option<String> {
    WANTS_BACKTRACE
        .load(Ordering::Relaxed)
        .then(|| std::backtrace::Backtrace::force_capture().to_string())
}

/// Marks an error as a defect in Rastair, see [`ThisIsABug`]. Transparent:
/// the error chain reads the same with it as without it.
#[derive(Debug)]
pub(crate) struct Bug(Box<dyn StdError + Send + Sync + 'static>);

impl fmt::Display for Bug {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl StdError for Bug {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        self.0.source()
    }
}

impl Bug {
    /// Whether `error`, or anything it wraps, went through [`ThisIsABug`].
    pub(crate) fn marks(error: &(dyn StdError + 'static)) -> bool {
        std::iter::successors(Some(error), |&error| error.source()).any(|error| error.is::<Self>())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bug_location_is_the_caller() {
        let line = line!() + 1;
        let error = Err::<(), _>(color_eyre::eyre::eyre!("boom")).this_is_a_bug();
        let report = format!("{:?}", error.expect_err("constructed as Err"));
        assert!(report.contains(&format!("{}:{line}", file!())), "{report}");
    }
}
