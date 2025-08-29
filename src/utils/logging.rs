use color_eyre::{Section as _, eyre::Report};
use tracing_subscriber::{EnvFilter, layer::SubscriberExt as _};

pub static LOG_VAR: &str = "RASTAIR_LOG";
pub static BUG_MESSAGE: &str = "This is a bug in Rastair, please report it at <https://bitbucket.org/bsblabludwig/rastair/issues/new>";

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
    if verbose {
        // SAFETY: This is set at the very start of the program
        unsafe { std::env::set_var("RUST_BACKTRACE", "1") };
    }
    let hook = color_eyre::config::HookBuilder::default()
        .panic_section(BUG_MESSAGE)
        .display_env_section(verbose)
        .display_location_section(verbose)
        .theme(if std::env::var("NO_COLOR").is_ok() {
            color_eyre::config::Theme::new()
        } else {
            color_eyre::config::Theme::dark()
        })
        .install()
        .note("Seeing this error message is somewhat ironic, we know");
    if let Err(error) = hook {
        eprintln!("Failed to register panic handler: {error:#}");
    }
}

pub trait ThisIsABug<T> {
    fn this_is_a_bug(self) -> Result<T, Report>;
}

impl<T, E> ThisIsABug<T> for Result<T, E>
where
    E: Into<Report>,
{
    fn this_is_a_bug(self) -> Result<T, Report> {
        use color_eyre::Help;

        self.map_err(|error| error.into()).map_err(|report| report.note(BUG_MESSAGE))
    }
}
