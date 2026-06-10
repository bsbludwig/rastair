//! Splice freshly-trained insertion/deletion forests into an existing
//! (pre-indel) model file while leaving its cpg/denovo/others forests untouched.
//!
//! This is a script that is only useful while our indel support is
//! experimental.
//!
//! ```
//!   cargo run --example splice_indel_models -- \
//!       --old models/rastair_default.rff.mpk.lz4 \
//!       --new models/rastair_indels.rff.mpk.lz4 \
//!       --output models/rastair_merged.rff.mpk.lz4
//! ```

use std::{fs, path::PathBuf};

use biosphere::FlatForest;
use clap::Parser;
use color_eyre::eyre::{Context as _, Result, ensure};
use lz4::{Decoder, EncoderBuilder};
use rastair::metrics::ml::types::{MlFeatureSet, PlattScaling, RastairFlatModel};

// Matches main 3c4f8515ab9ef863789e37c5c3da0f5cbb699472
#[derive(serde::Deserialize)]
struct OldGrouped {
    cpg: FlatForest,
    denovo: FlatForest,
    others: FlatForest,
    #[serde(default)]
    cpg_platt: PlattScaling,
    #[serde(default)]
    denovo_platt: PlattScaling,
    #[serde(default)]
    others_platt: PlattScaling,
    #[serde(default)]
    feature_set: MlFeatureSet,
}

#[derive(Parser)]
#[command(about = "Merge indel forests from a new model into an old (pre-indel) model file")]
struct Args {
    /// Pre-indel model file to take cpg/denovo/others forests from
    #[arg(long)]
    old: PathBuf,

    /// Newly trained model file to take insertion/deletion forests from
    #[arg(long)]
    new: PathBuf,

    /// Where to write the merged model
    #[arg(short = 'o', long)]
    output: PathBuf,
}

fn read_decompressed(path: &PathBuf) -> Result<Vec<u8>> {
    use std::io::Read as _;
    let bytes = fs::read(path).wrap_err_with(|| format!("Failed to read {}", path.display()))?;
    let mut decoder = Decoder::new(&bytes[..]).wrap_err("Failed to create LZ4 decoder")?;
    let mut out = Vec::new();
    decoder.read_to_end(&mut out).wrap_err("Failed to decompress LZ4 stream")?;
    Ok(out)
}

/// Load the pre-indel model, auto-detecting grouped vs interleaved field order.
fn load_old(path: &PathBuf) -> Result<OldGrouped> {
    let raw = read_decompressed(path)?;
    rmp_serde::from_slice(&raw).wrap_err("could not load old model")
}

fn load_new(path: &PathBuf) -> Result<RastairFlatModel> {
    let raw = read_decompressed(path)?;
    rmp_serde::from_slice(&raw)
        .wrap_err_with(|| format!("Failed to deserialize new model from {}", path.display()))
}

fn write_model(model: &RastairFlatModel, path: &PathBuf) -> Result<()> {
    let file =
        fs::File::create(path).wrap_err_with(|| format!("Failed to create {}", path.display()))?;
    let mut encoder =
        EncoderBuilder::new().level(16).build(file).wrap_err("Failed to create LZ4 encoder")?;
    rmp_serde::encode::write(&mut encoder, model).wrap_err("Failed to serialize merged model")?;
    let (_file, result) = encoder.finish();
    result.wrap_err("Failed to finalize LZ4 compression")?;
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    rastair::setup_logging(false);

    let old = load_old(&args.old)?;
    let new = load_new(&args.new)?;

    // All five forests are evaluated with one shared feature calculator, so the
    // forests we keep and the ones we splice in must agree on the feature set.
    ensure!(
        old.feature_set.to_string() == new.feature_set.to_string(),
        "Feature set mismatch: old model is `{}`, new model is `{}`. \
         The cpg/denovo/others and insertion/deletion forests would be scored with \
         incompatible features. Retrain both with the same --ml-features.",
        old.feature_set,
        new.feature_set,
    );

    let merged = RastairFlatModel {
        cpg: old.cpg,
        cpg_platt: old.cpg_platt,
        denovo: old.denovo,
        denovo_platt: old.denovo_platt,
        others: old.others,
        others_platt: old.others_platt,
        insertion: new.insertion,
        insertion_platt: new.insertion_platt,
        deletion: new.deletion,
        deletion_platt: new.deletion_platt,
        feature_set: old.feature_set,
    };

    write_model(&merged, &args.output)?;

    tracing::info!(
        "Merged model written to {}\n  cpg/denovo/others: {}\n  \
         insertion/deletion: {}\n  feature set: {}",
        args.output.display(),
        args.old.display(),
        args.new.display(),
        merged.feature_set,
    );

    Ok(())
}
