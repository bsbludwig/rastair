use crate::vcf::{DeNovoCpGCandidate, InCpG, Record};
use biosphere::RandomForest;
use color_eyre::{
    Result,
    eyre::{Context, ensure},
};
use ndarray::{Array1, Axis};
use rastair2_types::Base;
use smallvec::SmallVec;
use smol_str::SmolStr;
use std::{fmt, fs, io::Read, path::Path};
use tracing::{debug, instrument, warn};

pub struct MachineLearning {
    pub disabled: bool,
    pub threshold: f64,
    pub cpg: Option<Box<RandomForest>>,
    pub denovo_cpg: Option<Box<RandomForest>>,
    pub others: Option<Box<RandomForest>>,
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub enum MlResult {
    None,
    Predictions(SmallVec<Prediction, 1>),
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Prediction {
    pub model: MlModel,
    pub allele: Base,
    pub prediction: f64,
    pub threshold: f64,
    #[serde(skip)]
    pub features: Array1<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MlModel {
    Cpg,
    DenovoCpg,
    Others,
}

impl Prediction {
    pub fn pass(&self) -> bool {
        self.prediction >= self.threshold
    }

    pub fn empty() -> Self {
        Self {
            model: MlModel::Others,
            allele: Base::Unknown,
            prediction: 0.,
            threshold: 1.,
            features: Array1::from_elem(0, 0.),
        }
    }
}

impl fmt::Debug for Prediction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple(if self.pass() { "MlResult::PASS" } else { "MlResult::FAIL" })
            .field(&self.prediction)
            .field(&self.model)
            .finish()
    }
}

impl std::ops::Deref for Prediction {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.prediction
    }
}

impl MachineLearning {
    pub fn disabled() -> Self {
        Self { disabled: true, threshold: 1., cpg: None, denovo_cpg: None, others: None }
    }

    pub fn predict(
        &self,
        record: &Record,
        before: Option<&Record>,
        after: Option<&Record>,
    ) -> MlResult {
        if self.disabled {
            return MlResult::None;
        }

        // check all alts for this record and calculate predictions for each based on which model fits
        let mut predictions = SmallVec::new();
        for alt in &record.main.alt {
            let (name, model, features) = if let InCpG::C = record.info.in_cp_g
                && alt == "T"
            {
                (
                    MlModel::Cpg,
                    self.cpg.as_ref(),
                    super::cpg::params_from_record(record, before, after),
                )
            } else if let InCpG::G = record.info.in_cp_g
                && alt == "A"
            {
                (
                    MlModel::Cpg,
                    self.cpg.as_ref(),
                    super::cpg::params_from_record(record, before, after),
                )
            } else if let DeNovoCpGCandidate::Candidate { ref_base, alt_base, alt_index } =
                record.info.de_novo_cp_g_candidate
                && *alt == alt_base
            {
                (
                    MlModel::DenovoCpg,
                    self.denovo_cpg.as_ref(),
                    super::denovo_cpg::params_from_record(record, before, after),
                )
            } else {
                (
                    MlModel::Others,
                    self.others.as_ref(),
                    super::others::params_from_record(record, before, after),
                )
            };

            let Some(model) = model else {
                warn!(model=?name, "No model found");
                predictions.push(Prediction::empty());
                continue;
            };
            let prediction = model.predict(&features.view());
            match prediction.get(0).copied() {
                Some(p) => predictions.push(Prediction {
                    prediction: p,
                    threshold: self.threshold,
                    allele: alt.parse().unwrap_or(Base::Unknown),
                    features: features.row(0).to_owned(),
                    model: name,
                }),
                None => {
                    warn!(model=?name, "No predictions");
                    predictions.push(Prediction::empty())
                }
            }
        }

        MlResult::Predictions(predictions)
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
