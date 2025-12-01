//! Convert three separate model files into a single combined `RastairModel` file

use biosphere::RandomForest;
use clap::Parser;
use color_eyre::{Result, eyre::Context};
use rastair::setup_logging;
use std::{fs, path::PathBuf};
use tracing::info;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct RastairModel {
    pub cpg: RandomForest,
    pub denovo: RandomForest,
    pub others: RandomForest,
}

#[derive(Parser)]
#[command(name = "combine_models")]
#[command(about = "Combine three separate model files into one RastairModel file")]
struct Args {
    /// Path to the CpG model file
    #[arg(long)]
    cpg: PathBuf,

    /// Path to the denovo CpG model file
    #[arg(long)]
    denovo: PathBuf,

    /// Path to the others model file
    #[arg(long)]
    others: PathBuf,

    /// Output path for the combined model file
    #[arg(long)]
    output: PathBuf,
}

fn main() -> Result<()> {
    setup_logging(true);
    let args = Args::parse();

    info!("Loading CpG model from: {}", args.cpg.display());
    let cpg = load_model(&args.cpg)?;

    info!("Loading denovo CpG model from: {}", args.denovo.display());
    let denovo = load_model(&args.denovo)?;

    info!("Loading others model from: {}", args.others.display());
    let others = load_model(&args.others)?;

    let combined = RastairModel { cpg, denovo, others };

    info!("Saving combined model to: {}", args.output.display());
    save_combined_model(&combined, &args.output)?;

    info!("✓ Successfully created combined model file");

    Ok(())
}

fn load_model(path: &PathBuf) -> Result<RandomForest> {
    let file = fs::read(path)
        .wrap_err_with(|| format!("Failed to read model file: {}", path.display()))?;

    let decompress = lz4::Decoder::new(&file[..]).wrap_err("Failed to create LZ4 decoder")?;

    rmp_serde::from_read(decompress).wrap_err("Failed to deserialize random forest")
}

fn save_combined_model(model: &RastairModel, path: &PathBuf) -> Result<()> {
    let file = fs::File::create(path)
        .wrap_err_with(|| format!("Failed to create output file: {}", path.display()))?;

    let mut encoder = lz4::EncoderBuilder::new()
        .level(16)
        .build(file)
        .wrap_err("Failed to create LZ4 encoder")?;

    rmp_serde::encode::write(&mut encoder, model).wrap_err("Failed to serialize combined model")?;

    let (_output, result) = encoder.finish();
    result.wrap_err("Failed to finish LZ4 encoding")?;

    Ok(())
}
