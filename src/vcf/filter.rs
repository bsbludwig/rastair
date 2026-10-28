//! The closed set of VCF FILTER codes rastair can emit.
//!
//! Filters are modelled as a `#[repr(u8)]` enum rather than carried around as
//! strings. A record carries at most one: an SNV alt's
//! [`Rejection`](crate::metrics::Rejection), or an indel the pair-HMM does
//! not call. This enum is the single source of truth: the VCF header is registered from
//! [`RastairFilter::ALL`] (in declaration order) and each filter resolves to
//! its [`FilterId`](seqair::vcf::FilterId) by name at write time (see
//! [`crate::vcf::schema`]).
//!
//! `RastairFilter as usize` indexes the resolved-`FilterId` table in
//! [`crate::vcf::schema`], so the discriminants are `0..COUNT`.
//!
//! Rejected records are only written under `--all`, and every committed VCF
//! snapshot is all `PASS`; a few unit tests assert a FILTER code. When
//! changing filter emission, look at the column directly:
//!
//! ```bash
//! cargo run -- call --fasta-file=tests/data/test.fasta.gz tests/data/test.bam --all \
//!   | grep -v '^#' | awk -F'\t' '$7!="PASS"{print $7}' | sort | uniq -c | sort -rn
//! ```
//!
//! An SNV alt's code is decided once, in `set_alt_calls`.

/// A VCF FILTER code. `PASS` is a record without one, so it is not a variant
/// here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum RastairFilter {
    /// `pre_ml`
    PreMl,
    /// `m_highDp`
    MHighDp,
    /// `low_ml_score`
    LowMlScore,
    /// `low_indel_qual`
    LowIndelQual,
    /// `no_genotype_copy`
    NoGenotypeCopy,
}

impl RastairFilter {
    /// Number of filter variants. The discriminants are `0..COUNT`, so a
    /// `RastairFilter` doubles as an index into a `[T; COUNT]` lookup table.
    pub const COUNT: usize = Self::ALL.len();

    /// Every filter, in registration (header / dictionary) order.
    pub const ALL: [RastairFilter; 5] =
        [Self::PreMl, Self::MHighDp, Self::LowMlScore, Self::LowIndelQual, Self::NoGenotypeCopy];

    /// The filter's VCF ID, as it appears in the header and FILTER column.
    pub const fn name(self) -> &'static str {
        match self {
            Self::PreMl => "pre_ml",
            Self::MHighDp => "m_highDp",
            Self::LowMlScore => "low_ml_score",
            Self::LowIndelQual => "low_indel_qual",
            Self::NoGenotypeCopy => "no_genotype_copy",
        }
    }

    /// The filter's header description.
    pub const fn description(self) -> &'static str {
        match self {
            Self::PreMl => "Low amount of usable evidence, skipping ML",
            Self::MHighDp => "Excessive coverage for methylation candidate",
            Self::LowMlScore => "Machine Learning module prediction below threshold",
            Self::LowIndelQual => "Indel allele the pair-HMM does not call present",
            Self::NoGenotypeCopy => {
                "SNV allele the model passes but the pair-HMM's genotype holds no copy of"
            }
        }
    }
}
