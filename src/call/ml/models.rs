use biosphere::RandomForest;
use color_eyre::{
    Result,
    eyre::{Context, ensure},
};
use std::{fmt, fs, io::Read, path::Path};
use tracing::{debug, instrument};

use crate::vcf::Record;

pub struct MachineLearning {
    pub threshold: f64,
    pub cpg: Option<Box<RandomForest>>,
    pub denovo_cpg: Option<Box<RandomForest>>,
    pub others: Option<Box<RandomForest>>,
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub enum MlResult {
    None,
    Prediction { prediction: f64, threshold: f64 },
}

impl MlResult {
    pub fn pass(&self) -> bool {
        match self {
            MlResult::None => false,
            MlResult::Prediction { prediction, threshold } => prediction >= threshold,
        }
    }
}

impl fmt::Debug for MlResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MlResult::None => f.debug_tuple("MlResult::None").finish(),
            MlResult::Prediction { prediction, threshold } => f
                .debug_tuple(if self.pass() { "MlResult::PASS" } else { "MlResult::FAIL" })
                .field(prediction)
                .finish(),
        }
    }
}

impl std::ops::Deref for MlResult {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        match self {
            MlResult::None => &0.0,
            MlResult::Prediction { prediction, .. } => prediction,
        }
    }
}

impl MachineLearning {
    pub fn disabled() -> Self {
        Self { threshold: 1., cpg: None, denovo_cpg: None, others: None }
    }

    pub fn predict(
        &self,
        record: &Record,
        before: Option<&Record>,
        after: Option<&Record>,
    ) -> MlResult {
        let (model, features) = if *record.info.in_cp_g {
            (self.cpg.as_ref(), super::cpg::params_from_record(record, before, after))
        } else if *record.info.de_novo_cp_g_candidate {
            (self.denovo_cpg.as_ref(), super::denovo_cpg::params_from_record(record, before, after))
        } else {
            (self.others.as_ref(), super::others::params_from_record(record, before, after))
        };

        let Some(model) = self.cpg.as_ref() else {
            return MlResult::None;
        };

        let features = super::cpg::params_from_record(record, before, after);
        let prediction = model.predict(&features.view());
        match prediction.get(0).copied() {
            Some(p) => MlResult::Prediction { prediction: p, threshold: self.threshold },
            None => MlResult::None,
        }
    }
}

#[instrument(level = "debug", skip_all)]
pub fn load_model(path: Option<&Path>, built_in: &[u8]) -> Result<RandomForest> {
    let Some(path) = path else {
        return load_rf(built_in);
    };
    load_model_from_file(path)
        .wrap_err_with(|| format!("Failed to load model from: {}", path.display()))
}

fn load_model_from_file(path: &Path) -> Result<RandomForest> {
    ensure!(path.exists(), "Model file does not exist");
    let file = fs::read(path).wrap_err("Failed to read model file")?;
    load_rf(&file[..]).wrap_err("Failed to load model")
}

fn load_rf(reader: impl Read) -> Result<RandomForest> {
    let decompress = lz4::Decoder::new(reader).wrap_err("Failed to create LZ4 decoder")?;
    rmp_serde::from_read(decompress).wrap_err("Failed to deserialize random forest")
}
