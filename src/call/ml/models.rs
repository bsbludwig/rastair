use biosphere::RandomForest;
use color_eyre::{Result, eyre::Context};
use std::{fmt, io::Read};
use tracing::{debug, instrument};

use crate::vcf::Record;

pub struct MachineLearning {
    threshold: f64,
    cpg: Option<Box<RandomForest>>,
    denovo_cpg: Option<Box<RandomForest>>,
    others: Option<Box<RandomForest>>,
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

    pub fn with_threshold(threshold: f64) -> Self {
        Self {
            threshold,
            cpg: Some(Box::new(
                load_rf(&include_bytes!("../../../models/BS_RF_800-2_CpG.rf.mpk.lz4")[..])
                    .expect("Failed to load CpG RF model"),
            )),
            denovo_cpg: Some(Box::new(
                load_rf(&include_bytes!("../../../models/BS_RF_800-2_denovo.rf.mpk.lz4")[..])
                    .expect("Failed to load DeNovo CpG RF model"),
            )),
            others: Some(Box::new(
                load_rf(&include_bytes!("../../../models/BS_RF_800-2_other.rf.mpk.lz4")[..])
                    .expect("Failed to load Others RF model"),
            )),
        }
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
fn load_rf(reader: impl Read) -> Result<RandomForest> {
    let decompress = lz4::Decoder::new(reader).wrap_err("Failed to create LZ4 decoder")?;
    rmp_serde::from_read(decompress).wrap_err("Failed to deserialize random forest")
}
