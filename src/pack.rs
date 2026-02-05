use crate::{
    metrics::ml::types::{MlFeatureSet, RastairModel},
    utils::cli,
};
use biosphere::RandomForest;
use clio::ClioPath;
use color_eyre::eyre::{Context, Result, bail, ensure};
use lz4::{Decoder, EncoderBuilder};
use std::{fs, path::Path};
use tracing::instrument;

#[derive(Debug, clap::Args)]
pub struct PackModelParams {
    /// CpG random forest model (.mpk.lz4)
    #[arg(help_heading = cli::sections::INPUT, value_parser=clap::value_parser!(ClioPath).exists().is_file())]
    pub cpg: ClioPath,

    /// De-novo CpG random forest model (.mpk.lz4)
    #[arg(help_heading = cli::sections::INPUT, value_parser=clap::value_parser!(ClioPath).exists().is_file())]
    pub denovo: ClioPath,

    /// Other positions random forest model (.mpk.lz4)
    #[arg(help_heading = cli::sections::INPUT, value_parser=clap::value_parser!(ClioPath).exists().is_file())]
    pub others: ClioPath,

    /// Output combined model file (.mpk.lz4)
    #[arg(short = 'o', long = "output", default_value = "./models/rastair_combined.rf.mpk.lz4")]
    #[arg(help_heading = cli::sections::OUTPUT, value_hint=clap::ValueHint::FilePath)]
    pub output: ClioPath,

    /// Feature set used to train the input models
    #[arg(long, default_value_t = MlFeatureSet::Standard)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub ml_features: MlFeatureSet,

    /// Platt scaling parameter A for CpG model
    #[arg(long, default_value_t = 1.0)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub cpg_platt_a: f64,

    /// Platt scaling parameter B for CpG model
    #[arg(long, default_value_t = 0.0)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub cpg_platt_b: f64,

    /// Platt scaling parameter A for de-novo CpG model
    #[arg(long, default_value_t = 1.0)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub denovo_platt_a: f64,

    /// Platt scaling parameter B for de-novo CpG model
    #[arg(long, default_value_t = 0.0)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub denovo_platt_b: f64,

    /// Platt scaling parameter A for others model
    #[arg(long, default_value_t = 1.0)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub others_platt_a: f64,

    /// Platt scaling parameter B for others model
    #[arg(long, default_value_t = 0.0)]
    #[arg(help_heading = cli::sections::PROCESSING)]
    pub others_platt_b: f64,
}

#[instrument(level = "info", skip_all)]
pub fn pack_models(params: &PackModelParams) -> Result<()> {
    let cpg = load_random_forest(params.cpg.path()).wrap_err("Failed to load CpG model")?;
    let denovo =
        load_random_forest(params.denovo.path()).wrap_err("Failed to load de-novo CpG model")?;
    let others =
        load_random_forest(params.others.path()).wrap_err("Failed to load others model")?;

    let combined = RastairModel {
        cpg,
        denovo,
        others,
        cpg_platt: crate::metrics::ml::types::PlattScaling {
            a: params.cpg_platt_a,
            b: params.cpg_platt_b,
        },
        denovo_platt: crate::metrics::ml::types::PlattScaling {
            a: params.denovo_platt_a,
            b: params.denovo_platt_b,
        },
        others_platt: crate::metrics::ml::types::PlattScaling {
            a: params.others_platt_a,
            b: params.others_platt_b,
        },
        feature_set: params.ml_features,
    };

    params.output.parent().and_then(|p| std::fs::create_dir_all(p).ok());

    serialize_model(&combined, params.output.clone())
        .wrap_err_with(|| format!("Failed to serialize model to {}", params.output.display()))?;

    Ok(())
}

fn load_random_forest(path: &Path) -> Result<RandomForest> {
    ensure!(path.exists(), "Model file does not exist");
    let bytes = fs::read(path).wrap_err("Failed to read model file")?;

    match decode_random_forest(&bytes) {
        Ok(model) => Ok(model),
        Err(err) => {
            if decode_combined_model(&bytes).is_ok() {
                bail!(
                    "Expected a single RandomForest model, but {path:?} appears to be a combined Rastair model"
                );
            }
            Err(err).wrap_err_with(|| format!("Failed to deserialize RandomForest from {path:?}"))
        }
    }
}

fn decode_random_forest(bytes: &[u8]) -> Result<RandomForest> {
    let decoder = Decoder::new(bytes).wrap_err("Failed to create LZ4 decoder")?;
    rmp_serde::from_read(decoder).wrap_err("Failed to deserialize RandomForest model")
}

fn decode_combined_model(bytes: &[u8]) -> Result<RastairModel> {
    let decoder = Decoder::new(bytes).wrap_err("Failed to create LZ4 decoder")?;
    rmp_serde::from_read(decoder).wrap_err("Failed to deserialize combined model")
}

/// Serialize a model to disk with LZ4 compression
fn serialize_model(model: &RastairModel, path: ClioPath) -> Result<()> {
    let file = path.create().wrap_err("Failed to create output file for model serialization")?;
    let mut encoder =
        EncoderBuilder::new().level(16).build(file).wrap_err("Failed to create LZ4 encoder")?;

    rmp_serde::encode::write(&mut encoder, &model).wrap_err("Failed to serialize model")?;

    let (_output, result) = encoder.finish();
    result.wrap_err("Failed to finalize LZ4 compression")?;

    Ok(())
}
