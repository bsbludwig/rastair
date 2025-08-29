use crate::vcf;
use rastair_vcf::Contig;
use serde::{Deserialize, Serialize};
use smol_str::SmolStr;
use std::borrow::Cow;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MpkHeader {
    pub rastair2_version: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MpkVcfHeader {
    pub contigs: Vec<Contig>,
    pub samples: Vec<SmolStr>,
    pub metadata: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(clippy::large_enum_variant)] // all but two entries in a file are `Record`s
pub enum MpkEntry<'r> {
    Header(MpkHeader),
    VcfHeader(MpkVcfHeader),
    Record(Cow<'r, vcf::Record>),
}
