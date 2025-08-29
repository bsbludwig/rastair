use color_eyre::eyre::Report;
use tracing_subscriber::{EnvFilter, layer::SubscriberExt as _};

pub fn setup_tracing(verbose: bool) {
    let subscriber = {
        let default_log_settings =
            if verbose { "info,rastair2=debug" } else { "warn,rastair2=info" };
        let mut env_filter = EnvFilter::new(default_log_settings);
        if let Ok(env) = std::env::var("RASTAIR_LOG") {
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

pub static BUG_MESSAGE: &str = "This is a bug in Rastair, please report it at <https://bitbucket.org/bsblabludwig/rastair/issues/new>";

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
