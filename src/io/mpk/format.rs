use crate::metrics::PileupMetrics;
use rastair_vcf::Contig;
use seqair_types::SmolStr;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

/// Serialization format version of the `.mpk` payload.
///
/// Bump this whenever the serialized shape of anything reachable from
/// [`MpkEntry`] changes — a field added, removed or reordered on [`PileupMetrics`]
/// or any of its components. `rmp_serde::encode::write` emits structs as
/// *positional* arrays, so such a change shifts every field after it and
/// `#[serde(default)]` cannot recover: a file written by another version would
/// otherwise decode into neighbouring fields rather than fail. The reader
/// compares this value and refuses mismatches.
pub const MPK_FORMAT_VERSION: u32 = 2;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MpkHeader {
    pub rastair_version: String,
    pub format_version: u32,
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
