use crate::metrics::PileupMetrics;
use rastair_vcf::Contig;
use seqair_types::SmolStr;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MpkHeader {
    pub rastair_version: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MpkVcfHeader {
    pub contigs: Vec<Contig>,
    pub samples: Vec<SmolStr>,
    pub metadata: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[expect(clippy::large_enum_variant, reason = "headers only occur once")]
pub enum MpkEntry<'r> {
    Header(MpkHeader),
    VcfHeader(MpkVcfHeader),
    Record(Cow<'r, PileupMetrics>),
}
