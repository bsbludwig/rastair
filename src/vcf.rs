//! VCF/BCF schema and direct record emission.
//!
//! The schema (every INFO/FORMAT/FILTER field) is defined once in [`schema`]
//! and encoded straight from [`PileupMetrics`](crate::metrics::PileupMetrics)
//! by [`emit`] using seqair's typestate writer — there is no intermediate
//! record struct. The domain field types (`InCpG`, `Methylated`, …) live in
//! their own modules and are read directly during emission.
//!
//! See
//! <https://github.com/samtools/hts-specs/blob/0d7f8774658f7cee0a4540b0682174e460726432/VCFv4.5.tex>
//! for the VCF spec.

use seqair_types::SmolStr;

mod cpg;
pub use cpg::InCpG;
mod denovo_cpg;
pub use denovo_cpg::DeNovoCpGCandidate;
mod methylation;
pub use methylation::{CpgBeta, CpgOrigin, Methylated, MethylationAltDepth, MethylationDepth};
mod filter;
pub use filter::RastairFilter;

pub mod schema;
pub use schema::{FieldConfig, FormatFieldId, InfoFieldId, Schema, register};

mod emit;
pub use emit::emit_pileup;

pub use crate::metrics::MethylationEvidenceStrandInfo;
pub use crate::utils::SequenceContext;

/// A contig (chromosome) entry for the VCF header.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize)]
pub struct Contig {
    /// Name of the contig (chromosome).
    pub name: SmolStr,
    /// Length of the contig.
    pub length: u64,
}
