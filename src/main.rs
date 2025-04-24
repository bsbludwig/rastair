use clap::Parser as _;
use color_eyre::eyre::Result;
use rastair2::call::{CallParams, call};
use tracing::debug;
use tracing_subscriber::{fmt::format::FmtSpan, layer::SubscriberExt as _};

#[derive(Debug, clap::Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Subcommand,
}

#[derive(Debug, clap::Subcommand)]
enum Subcommand {
    Call(CallParams),
}

fn main() -> Result<()> {
    color_eyre::install()?;
    let subscriber =
        tracing_subscriber::Registry::default().with(tracing_error::ErrorLayer::default()).with(
            tracing_subscriber::fmt::Layer::default()
                .with_target(true)
                .with_span_events(FmtSpan::CLOSE)
                .with_writer(std::io::stderr),
        );
    tracing::subscriber::set_global_default(subscriber)?;

    let args = Cli::parse();

    match args.command {
        Subcommand::Call(params) => {
            debug!(?params, "Running call command");
            call(&params)?;
        }
    }

    Ok(())
}

// bitflags! {
//     #[derive(Debug, Clone, Copy)]
//     struct Flags: u16 {
//         const IS_PAIRED = 0x1;
//         const IS_PROPERLY_PAIRED = 0x2;
//         const IS_UNMAPPED = 0x4;
//         const MATE_IS_UNMAPPED = 0x8;
//         const IS_REVERSE_STRAND = 0x10;
//         const MATE_IS_REVERSE_STRAND = 0x20;
//         const IS_FIRST_IN_PAIR = 0x40;
//         const IS_SECOND_IN_PAIR = 0x80;
//         const IS_NOT_PRIMARY = 0x100;
//         const IS_FAILED = 0x200;
//         const IS_DUPLICATE = 0x400;
//         const IS_SUPPLEMENTAL = 0x800;
//     }
// }
