use crate::utils::cli;
use better_default::Default;
use seqair_types::{BaseQuality, MappingQuality};

#[derive(Debug, Clone, Default, clap::Args, serde::Serialize, serde::Deserialize)]
#[serde(from = "QualityBytes", into = "QualityBytes")]
pub struct QualityFilterParams {
    /// Minimum mapping quality to consider a read
    #[arg(short = 'q', long, default_value = "1", value_parser = mapping_quality)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(MappingQuality::from_byte(1))]
    pub min_mapq: MappingQuality,
    /// Minimum base quality to consider a base
    #[arg(short = 'Q', long, default_value = "10", value_parser = base_quality)]
    #[arg(help_heading = cli::sections::FILTER)]
    #[default(BaseQuality::from_byte(10))]
    pub min_baseq: BaseQuality,
}

impl QualityFilterParams {
    /// Whether a base of quality `baseq` counts. An unavailable base quality
    /// never does, and a minimum of 255 is above every available one.
    #[must_use]
    pub fn keeps_base(&self, baseq: BaseQuality) -> bool {
        baseq.get().is_some_and(|q| q >= self.min_baseq.as_byte())
    }
}

fn mapping_quality(s: &str) -> Result<MappingQuality, std::num::ParseIntError> {
    s.parse().map(MappingQuality::from_byte)
}

fn base_quality(s: &str) -> Result<BaseQuality, std::num::ParseIntError> {
    s.parse().map(BaseQuality::from_byte)
}

/// The thresholds as the run's config records them: their bytes.
#[derive(serde::Serialize, serde::Deserialize)]
struct QualityBytes {
    min_mapq: u8,
    min_baseq: u8,
}

impl From<QualityBytes> for QualityFilterParams {
    fn from(bytes: QualityBytes) -> Self {
        Self {
            min_mapq: MappingQuality::from_byte(bytes.min_mapq),
            min_baseq: BaseQuality::from_byte(bytes.min_baseq),
        }
    }
}

impl From<QualityFilterParams> for QualityBytes {
    fn from(params: QualityFilterParams) -> Self {
        Self { min_mapq: params.min_mapq.as_byte(), min_baseq: params.min_baseq.as_byte() }
    }
}
