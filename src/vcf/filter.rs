//! The closed set of VCF FILTER codes rastair can emit.
//!
//! Filters are modelled as a `#[repr(u8)]` enum rather than carried around as
//! strings, so storage on [`PileupMetrics`](crate::metrics::PileupMetrics) is a
//! single byte per filter. This enum is the single source of truth: the VCF
//! header is registered from [`RastairFilter::ALL`] (in declaration order) and
//! each filter resolves to its [`FilterId`](seqair::vcf::FilterId) by name at
//! write time (see [`crate::vcf::schema`]).

/// A VCF FILTER code. `PASS` is implicit (an empty filter set), so it is not a
/// variant here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[repr(u8)]
pub enum RastairFilter {
    /// `lowDp`
    LowDp,
    /// `dnCpG_lowDp`
    DnCpgLowDp,
    /// `dnCpG_bq`
    DnCpgBq,
    /// `dnCpG_mapq`
    DnCpgMapq,
    /// `dnCpG_vaf`
    DnCpgVaf,
    /// `dnCpG_adj`
    DnCpgAdj,
    /// `m_vaf`
    MVaf,
    /// `m_bq_ratio`
    MBqRatio,
    /// `m_pos`
    MPos,
    /// `m_highDp`
    MHighDp,
    /// `pre_ml`
    PreMl,
    /// `low_ml_score`
    LowMlScore,
    /// `indel_strand`
    IndelStrand,
}

impl RastairFilter {
    /// Number of filter variants. The discriminants are `0..COUNT`, so a
    /// `RastairFilter` doubles as an index into a `[T; COUNT]` lookup table.
    pub const COUNT: usize = Self::ALL.len();

    /// Every filter, in registration (header / dictionary) order.
    pub const ALL: [RastairFilter; 13] = [
        Self::LowDp,
        Self::DnCpgLowDp,
        Self::DnCpgBq,
        Self::DnCpgMapq,
        Self::DnCpgVaf,
        Self::DnCpgAdj,
        Self::MVaf,
        Self::MBqRatio,
        Self::MPos,
        Self::MHighDp,
        Self::PreMl,
        Self::LowMlScore,
        Self::IndelStrand,
    ];

    /// The filter's VCF ID, as it appears in the header and FILTER column.
    pub const fn name(self) -> &'static str {
        match self {
            Self::LowDp => "lowDp",
            Self::DnCpgLowDp => "dnCpG_lowDp",
            Self::DnCpgBq => "dnCpG_bq",
            Self::DnCpgMapq => "dnCpG_mapq",
            Self::DnCpgVaf => "dnCpG_vaf",
            Self::DnCpgAdj => "dnCpG_adj",
            Self::MVaf => "m_vaf",
            Self::MBqRatio => "m_bq_ratio",
            Self::MPos => "m_pos",
            Self::MHighDp => "m_highDp",
            Self::PreMl => "pre_ml",
            Self::LowMlScore => "low_ml_score",
            Self::IndelStrand => "indel_strand",
        }
    }

    /// The filter's header description.
    pub const fn description(self) -> &'static str {
        match self {
            Self::LowDp => "Low read depth",
            Self::DnCpgLowDp => "Low read depth for de-novo CpG candidate",
            Self::DnCpgBq => "Low base quality for de-novo CpG candidate",
            Self::DnCpgMapq => "Low mapping quality for de-novo CpG candidate",
            Self::DnCpgVaf => "Low variant allele frequency for de-novo CpG candidate",
            Self::DnCpgAdj => {
                "Included as adjacent position for de-novo CpG candidate, but other position did not pass filters"
            }
            Self::MVaf => "Low variant allele frequency for methylation candidate",
            Self::MBqRatio => "Low quality ratio for methylation candidate",
            Self::MPos => "Alt allele evidence from read edges for methylation candidate",
            Self::MHighDp => "Excessive coverage for methylation candidate",
            Self::PreMl => "Low amount of usable evidence, skipping ML",
            Self::LowMlScore => "Machine Learning module prediction below threshold",
            Self::IndelStrand => "Indel allele supported on only one bisulfite strand",
        }
    }
}
