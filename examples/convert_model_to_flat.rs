//! Convert a combined `RastairModel` file to a flat `RastairFlatModel` file.
//!
//! The flat model uses BFS-linearised f32 forests for faster CPU inference.

use biosphere::{FlatForest, RandomForest};
use clap::Parser;
use color_eyre::{
    Result,
    eyre::{Context, ensure},
};
use rastair::{
    metrics::ml::types::{MlFeatureSet, PlattScaling, RastairFlatModel},
    setup_logging,
};
use std::{
    fs,
    io::Read,
    path::{Path, PathBuf},
};
use tracing::{debug, info};

#[derive(Parser)]
#[command(name = "convert_model_to_flat")]
#[command(about = "Convert a combined RastairModel file to a flat RastairFlatModel file")]
struct Args {
    /// Path to the input combined model file (.rf.mpk.lz4)
    #[arg(long)]
    input: PathBuf,

    /// Output path for the flat model file
    #[arg(long)]
    output: PathBuf,
}

fn main() -> Result<()> {
    setup_logging(true);
    let args = Args::parse();

    info!("Loading combined model from: {}", args.input.display());
    let combined =
        load_combined_model(Some(&args.input), &[]).wrap_err("Failed to load combined model")?;

    let feature_nums = combined.feature_set.get_calculator().feature_num();

    info!("Converting forests to flat representation");
    let flat_cpg = FlatForest::from_forest(&combined.cpg, feature_nums.cpg);
    let flat_denovo = FlatForest::from_forest(&combined.denovo, feature_nums.denovo_cpg);
    let flat_others = FlatForest::from_forest(&combined.others, feature_nums.others);

    let flat_model = RastairFlatModel {
        cpg: flat_cpg,
        denovo: flat_denovo,
        others: flat_others,
        cpg_platt: combined.cpg_platt,
        denovo_platt: combined.denovo_platt,
        others_platt: combined.others_platt,
        feature_set: combined.feature_set,
    };

    info!("Saving flat model to: {}", args.output.display());
    save_flat_model(&flat_model, &args.output)?;

    info!(path=?args.output, "Successfully wrote flat model file");

    Ok(())
}

pub fn load_combined_model(path: Option<&Path>, built_in: &[u8]) -> Result<RastairModel> {
    fn load_model_from_file(path: &Path) -> Result<RastairModel> {
        ensure!(path.exists(), "Model file does not exist");
        let file = fs::read(path).wrap_err("Failed to read model file")?;
        load_combined(&file[..]).wrap_err("Failed to load combined model")
    }

    fn load_combined(reader: impl Read) -> Result<RastairModel> {
        let decompress = lz4::Decoder::new(reader).wrap_err("Failed to create LZ4 decoder")?;
        rmp_serde::from_read(decompress).wrap_err("Failed to deserialize combined model")
    }

    if let Some(path) = path {
        let model = load_model_from_file(path)
            .wrap_err_with(|| format!("Failed to load combined model from {path:?}"));
        debug!(?path, "Loaded combined model from file");
        model
    } else {
        let model = load_combined(built_in);
        debug!("Loaded built-in combined model");
        model
    }
}

/// Combined model file containing all three random forest models
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct RastairModel {
    pub cpg: RandomForest,
    pub denovo: RandomForest,
    pub others: RandomForest,
    #[serde(default)]
    pub cpg_platt: PlattScaling,
    #[serde(default)]
    pub denovo_platt: PlattScaling,
    #[serde(default)]
    pub others_platt: PlattScaling,
    #[serde(default)]
    pub feature_set: MlFeatureSet,
}

fn save_flat_model(model: &RastairFlatModel, path: &PathBuf) -> Result<()> {
    let file = fs::File::create(path)
        .wrap_err_with(|| format!("Failed to create output file: {}", path.display()))?;

    let mut encoder = lz4::EncoderBuilder::new()
        .level(16)
        .build(file)
        .wrap_err("Failed to create LZ4 encoder")?;

    rmp_serde::encode::write(&mut encoder, model).wrap_err("Failed to serialize flat model")?;

    let (_output, result) = encoder.finish();
    result.wrap_err("Failed to finish LZ4 encoding")?;

    Ok(())
}
