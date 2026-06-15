//! seqair-native VCF/BCF schema: the single source of truth for every INFO,
//! FORMAT and FILTER field rastair emits.
//!
//! Each field's header definition lives once as a `const *_DEF`. [`register`]
//! walks those defs through seqair's typestate header builder and returns a
//! [`Schema`] holding the resolved, pre-indexed keys plus contig/filter lookup
//! tables. [`crate::vcf::emit`] drives the encoder directly from
//! [`PileupMetrics`](crate::metrics::PileupMetrics) using these keys — there is
//! no intermediate record struct.

use std::io::Write;

use color_eyre::{Result, eyre::Context as _};
use rustc_hash::FxHashMap;
use seqair::vcf::record_encoder::{
    Arr, FieldDescription, FilterFieldDef, FormatFieldDef, InfoFieldDef, Scalar, Str,
};
use seqair::vcf::{
    ContigId, FilterId, FormatFloat, FormatFloats, FormatGt, FormatInt, FormatInts, InfoFlag,
    InfoFloat, InfoFloats, InfoInt, InfoInts, InfoString, Number, ValueType, VcfHeader,
};
use seqair_types::SmolStr;

use crate::vcf::{Contig, RastairFilter};

// ── INFO field definitions ──────────────────────────────────────────────
// Field order here is the VCF header order; it is cosmetic (seqair resolves
// keys by name), so it only needs to stay readable.

const AD_DEF: InfoFieldDef<Arr<i32>> = InfoFieldDef::new(
    "AD",
    Number::ReferenceAlternateBases,
    ValueType::Integer,
    "Total read depth for each allele",
);
const BQ_DEF: InfoFieldDef<Scalar<f32>> =
    InfoFieldDef::new("BQ", Number::Count(1), ValueType::Float, "RMS base quality");
const DP_DEF: InfoFieldDef<Scalar<i32>> =
    InfoFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Combined depth across samples");
const MQ_DEF: InfoFieldDef<Scalar<f32>> =
    InfoFieldDef::new("MQ", Number::Count(1), ValueType::Float, "RMS mapping quality");
const MQ0_DEF: InfoFieldDef<Scalar<i32>> =
    InfoFieldDef::new("MQ0", Number::Count(1), ValueType::Integer, "Number of MAPQ == 0 reads");
const NS_DEF: InfoFieldDef<Scalar<i32>> =
    InfoFieldDef::new("NS", Number::Count(1), ValueType::Integer, "Number of samples with data");
const AS_SB_OT_DEF: InfoFieldDef<Arr<i32>> = InfoFieldDef::new(
    "AS_SB_OT",
    Number::ReferenceAlternateBases,
    ValueType::Integer,
    "OT counts per allele",
);
const AS_SB_OB_DEF: InfoFieldDef<Arr<i32>> = InfoFieldDef::new(
    "AS_SB_OB",
    Number::ReferenceAlternateBases,
    ValueType::Integer,
    "OB counts per allele",
);
const SC5_DEF: InfoFieldDef<Str> = InfoFieldDef::new(
    "SC5",
    Number::Count(1),
    ValueType::String,
    "5-base sequence context centered on the variant position",
);
const AF_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "AF",
    Number::AlternateBases,
    ValueType::Float,
    "Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)",
);
const ABQ_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "ABQ",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "RMS Base quality per allele",
);
const AMQ_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "AMQ",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "RMS Map quality per allele",
);
const AS_SS_BQ_OT_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "AS_SS_BQ_OT",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "Strand-specific RMS of base quality per allele on the original top strand",
);
const AS_SS_BQ_OB_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "AS_SS_BQ_OB",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "Strand-specific RMS of base quality per allele on the original bottom strand",
);
const AS_SS_MQ_OT_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "AS_SS_MQ_OT",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "Strand-specific RMS of mapping quality per allele on the original top strand",
);
const AS_SS_MQ_OB_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "AS_SS_MQ_OB",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "Strand-specific RMS of mapping quality per allele on the original bottom strand",
);
const PIR_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "PIR",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "RMS of relative position in read",
);
const ENT100_DEF: InfoFieldDef<Scalar<f32>> = InfoFieldDef::new(
    "ENT100",
    Number::Count(1),
    ValueType::Float,
    "Shannon entropy of 100bp sequence context around variant position. Value range (0..2)",
);
const NAB_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "NAB",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "RMS of number of aligned bases",
);
const NOI_DEF: InfoFieldDef<Arr<f32>> = InfoFieldDef::new(
    "NOI",
    Number::ReferenceAlternateBases,
    ValueType::Float,
    "RMS of number of indels",
);
const M5MC_STRANDS_DEF: InfoFieldDef<Arr<i32>> = InfoFieldDef::new(
    "M5mC_Strands",
    Number::Count(4),
    ValueType::Integer,
    "Number of reads that are evidence for unmodified, modified, no SNP, snp",
);
const CPG_DEF: InfoFieldDef<seqair::vcf::Flag> =
    InfoFieldDef::new("CPG", Number::Count(0), ValueType::Flag, "Is this a CpG site?");
const CPGNOVO_DEF: InfoFieldDef<seqair::vcf::Flag> = InfoFieldDef::new(
    "CPGnovo",
    Number::Count(0),
    ValueType::Flag,
    "De-novo CPG candidate: Could the alt alleles create a new CpG site?",
);

// ── FORMAT field definitions ────────────────────────────────────────────

const GT_DEF: FormatFieldDef<seqair::vcf::Gt> =
    FormatFieldDef::new("GT", Number::Count(1), ValueType::String, "Genotype");
const GL_DEF: FormatFieldDef<Scalar<f32>> = FormatFieldDef::new(
    "GL",
    Number::Genotypes,
    ValueType::Float,
    "Genotype likelihoods, Phred-scaled",
);
const GC_DEF: FormatFieldDef<Scalar<f32>> = FormatFieldDef::new(
    "GC",
    Number::Genotypes,
    ValueType::Float,
    "Genotype confidence, Phred-scaled",
);
const SAMPLE_DP_DEF: FormatFieldDef<Scalar<i32>> =
    FormatFieldDef::new("DP", Number::Count(1), ValueType::Integer, "Read depth");
const M5MC_DEF: FormatFieldDef<Arr<f32>> = FormatFieldDef::new(
    "M5mC",
    Number::BaseModification,
    ValueType::Float,
    "Methylation level at CpG sites",
);
const DPM5MC_DEF: FormatFieldDef<Arr<i32>> = FormatFieldDef::new(
    "DPM5mC",
    Number::BaseModification,
    ValueType::Integer,
    "Total read depth for 5-methylcytosine detection",
);
const ADM5MC_DEF: FormatFieldDef<Arr<i32>> = FormatFieldDef::new(
    "ADM5mC",
    Number::BaseModification,
    ValueType::Integer,
    "Read depth supporting 5-methylcytosine modification",
);
const ML_DEF: FormatFieldDef<Arr<f32>> = FormatFieldDef::new(
    "ML",
    Number::AlternateBases,
    ValueType::Float,
    "Prediction of methylation/variant likelihood by rastair's machine learning model",
);

// ── FILTER definitions ──────────────────────────────────────────────────
// The FILTER set is owned by `RastairFilter` (the type stored on
// `PileupMetrics`); the header is registered from `RastairFilter::ALL` so the
// names/descriptions can never drift from the enum. PASS is registered
// automatically by seqair.

/// Every INFO field definition, in header order (both columns of dual-column
/// fields appear separately). Drives doc generation.
const INFO_DEFS: &[&dyn FieldDescription] = &[
    &AD_DEF,
    &BQ_DEF,
    &DP_DEF,
    &MQ_DEF,
    &MQ0_DEF,
    &NS_DEF,
    &AS_SB_OT_DEF,
    &AS_SB_OB_DEF,
    &SC5_DEF,
    &AF_DEF,
    &ABQ_DEF,
    &AMQ_DEF,
    &AS_SS_BQ_OT_DEF,
    &AS_SS_BQ_OB_DEF,
    &AS_SS_MQ_OT_DEF,
    &AS_SS_MQ_OB_DEF,
    &PIR_DEF,
    &ENT100_DEF,
    &NAB_DEF,
    &NOI_DEF,
    &M5MC_STRANDS_DEF,
    &CPG_DEF,
    &CPGNOVO_DEF,
];

/// Every FORMAT field definition, in header order. Drives doc generation.
const FORMAT_DEFS: &[&dyn FieldDescription] =
    &[&GT_DEF, &GL_DEF, &GC_DEF, &SAMPLE_DP_DEF, &M5MC_DEF, &DPM5MC_DEF, &ADM5MC_DEF, &ML_DEF];

/// Write the VCF field reference documentation as markdown.
pub fn write_vcf_docs<W: Write>(w: &mut W) -> std::io::Result<()> {
    writeln!(w, "# VCF output fields\n")?;

    writeln!(w, "## FILTER\n")?;
    writeln!(w, "| ID | Description |")?;
    writeln!(w, "|----|-------------|")?;
    writeln!(w, "| PASS | All filters passed |")?;
    for f in RastairFilter::ALL {
        writeln!(w, "| {} | {} |", f.name(), f.description())?;
    }

    let table = |w: &mut W, title: &str, defs: &[&dyn FieldDescription]| -> std::io::Result<()> {
        writeln!(w, "\n## {title}\n")?;
        writeln!(w, "| ID | Number | Type | Description |")?;
        writeln!(w, "|----|--------|------|-------------|")?;
        for d in defs {
            writeln!(
                w,
                "| {} | {} | {} | {} |",
                d.name(),
                d.number().as_str(),
                d.value_type().as_str(),
                d.description()
            )?;
        }
        Ok(())
    };
    table(w, "INFO", INFO_DEFS)?;
    table(w, "FORMAT", FORMAT_DEFS)?;
    Ok(())
}

/// All INFO field IDs, including both columns of dual-column fields. Used to
/// validate `--vcf-info-fields` CLI input.
pub const ALL_INFO_IDS: &[&str] = &[
    "AD",
    "BQ",
    "DP",
    "MQ",
    "MQ0",
    "NS",
    "AS_SB",
    "SC5",
    "AF",
    "ABQ",
    "AMQ",
    "AS_SS_BQ",
    "AS_SS_MQ",
    "PIR",
    "ENT100",
    "NAB",
    "NOI",
    "M5mC_Strands",
    "CPG",
    "CPGnovo",
];

/// All FORMAT field IDs. Used to validate `--vcf-format-fields` CLI input.
pub const ALL_FORMAT_IDS: &[&str] = &["GT", "GL", "GC", "DP", "M5mC", "DPM5mC", "ADM5mC", "ML"];

// ── Resolved keys ───────────────────────────────────────────────────────

/// All INFO keys, resolved against the header (BCF dict indices pre-computed).
pub(crate) struct InfoKeys {
    pub ad: InfoInts,
    pub bq: InfoFloat,
    pub dp: InfoInt,
    pub mq: InfoFloat,
    pub mq0: InfoInt,
    pub ns: InfoInt,
    pub as_sb_ot: InfoInts,
    pub as_sb_ob: InfoInts,
    pub sc5: InfoString,
    pub af: InfoFloats,
    pub abq: InfoFloats,
    pub amq: InfoFloats,
    pub as_ss_bq_ot: InfoFloats,
    pub as_ss_bq_ob: InfoFloats,
    pub as_ss_mq_ot: InfoFloats,
    pub as_ss_mq_ob: InfoFloats,
    pub pir: InfoFloats,
    pub ent100: InfoFloat,
    pub nab: InfoFloats,
    pub noi: InfoFloats,
    pub m5mc_strands: InfoInts,
    pub cpg: InfoFlag,
    pub cpgnovo: InfoFlag,
}

/// All FORMAT keys, resolved against the header.
pub(crate) struct FormatKeys {
    pub gt: FormatGt,
    pub gl: FormatFloat,
    pub gc: FormatFloat,
    pub dp: FormatInt,
    pub m5mc: FormatFloats,
    pub dpm5mc: FormatInts,
    pub adm5mc: FormatInts,
    pub ml: FormatFloats,
}

/// Resolved VCF schema: typed keys plus contig and filter lookup tables.
pub struct Schema {
    pub(crate) info: InfoKeys,
    pub(crate) format: FormatKeys,
    contigs: FxHashMap<SmolStr, ContigId>,
    /// Pre-resolved `FilterId`s indexed by `RastairFilter as usize`, so a
    /// filter resolves to its handle by a single array index instead of a
    /// per-record string hash.
    filters: [FilterId; RastairFilter::COUNT],
}

impl Schema {
    /// Resolve a contig name to its pre-indexed handle.
    pub(crate) fn contig(&self, name: &str) -> Option<&ContigId> {
        self.contigs.get(name)
    }

    /// Resolve a [`RastairFilter`] to its pre-indexed BCF handle.
    #[expect(
        clippy::indexing_slicing,
        reason = "RastairFilter discriminants are 0..COUNT, matching the array length"
    )]
    pub(crate) fn filter(&self, filter: RastairFilter) -> &FilterId {
        &self.filters[filter as usize]
    }
}

/// Build the VCF header and resolve all field/contig/filter keys.
///
/// Threads the typestate header builder through every phase
/// (contigs → filters → INFO → FORMAT → samples), registering each field once.
pub fn register(
    contigs: &[Contig],
    samples: &[SmolStr],
    metadata: &[String],
) -> Result<(VcfHeader, Schema)> {
    let mut builder = VcfHeader::builder();
    for line in metadata {
        builder.add_other_line(SmolStr::from(line.as_str()));
    }

    let mut contig_ids = FxHashMap::default();
    for contig in contigs {
        let id = builder
            .register_contig(
                contig.name.clone(),
                seqair::vcf::ContigDef { length: Some(contig.length) },
            )
            .wrap_err_with(|| format!("Failed to register contig {}", contig.name))?;
        contig_ids.insert(contig.name.clone(), id);
    }

    let mut builder = builder.filters();
    // Registered in `RastairFilter::ALL` order, which is discriminant order, so
    // the resulting Vec indexes line up with `RastairFilter as usize`.
    let mut filter_ids: Vec<FilterId> = Vec::with_capacity(RastairFilter::COUNT);
    for filter in RastairFilter::ALL {
        let def = FilterFieldDef::new(filter.name(), filter.description());
        filter_ids.push(builder.register_filter(&def).wrap_err("Failed to register filter")?);
    }
    let filters: [FilterId; RastairFilter::COUNT] =
        filter_ids.try_into().map_err(|_| color_eyre::eyre::eyre!("filter count mismatch"))?;

    let mut builder = builder.infos();
    let info = InfoKeys {
        ad: builder.register_info(&AD_DEF)?,
        bq: builder.register_info(&BQ_DEF)?,
        dp: builder.register_info(&DP_DEF)?,
        mq: builder.register_info(&MQ_DEF)?,
        mq0: builder.register_info(&MQ0_DEF)?,
        ns: builder.register_info(&NS_DEF)?,
        as_sb_ot: builder.register_info(&AS_SB_OT_DEF)?,
        as_sb_ob: builder.register_info(&AS_SB_OB_DEF)?,
        sc5: builder.register_info(&SC5_DEF)?,
        af: builder.register_info(&AF_DEF)?,
        abq: builder.register_info(&ABQ_DEF)?,
        amq: builder.register_info(&AMQ_DEF)?,
        as_ss_bq_ot: builder.register_info(&AS_SS_BQ_OT_DEF)?,
        as_ss_bq_ob: builder.register_info(&AS_SS_BQ_OB_DEF)?,
        as_ss_mq_ot: builder.register_info(&AS_SS_MQ_OT_DEF)?,
        as_ss_mq_ob: builder.register_info(&AS_SS_MQ_OB_DEF)?,
        pir: builder.register_info(&PIR_DEF)?,
        ent100: builder.register_info(&ENT100_DEF)?,
        nab: builder.register_info(&NAB_DEF)?,
        noi: builder.register_info(&NOI_DEF)?,
        m5mc_strands: builder.register_info(&M5MC_STRANDS_DEF)?,
        cpg: builder.register_info(&CPG_DEF)?,
        cpgnovo: builder.register_info(&CPGNOVO_DEF)?,
    };

    let mut builder = builder.formats();
    let format = FormatKeys {
        gt: builder.register_format(&GT_DEF)?,
        gl: builder.register_format(&GL_DEF)?,
        gc: builder.register_format(&GC_DEF)?,
        dp: builder.register_format(&SAMPLE_DP_DEF)?,
        m5mc: builder.register_format(&M5MC_DEF)?,
        dpm5mc: builder.register_format(&DPM5MC_DEF)?,
        adm5mc: builder.register_format(&ADM5MC_DEF)?,
        ml: builder.register_format(&ML_DEF)?,
    };

    let mut builder = builder.samples();
    for sample in samples {
        builder.add_sample(sample.clone()).wrap_err("Failed to add sample")?;
    }

    let header = builder.build().wrap_err("Failed to build VCF header")?;
    Ok((header, Schema { info, format, contigs: contig_ids, filters }))
}

// ── Field selection (CLI `--vcf-info-fields` / `--vcf-format-fields`) ────

/// Which INFO/FORMAT fields to write. Defaults to a minimal set; the rest can
/// be enabled per field via the CLI.
#[derive(Debug, Clone)]
pub struct FieldConfig {
    pub info: InfoSelection,
    pub format: FormatSelection,
}

/// Per-field INFO selection flags.
#[derive(Debug, Clone)]
pub struct InfoSelection {
    pub ad: bool,
    pub bq: bool,
    pub dp: bool,
    pub mq: bool,
    pub mq0: bool,
    pub ns: bool,
    pub as_sb: bool,
    pub sc5: bool,
    pub af: bool,
    pub abq: bool,
    pub amq: bool,
    pub as_ss_bq: bool,
    pub as_ss_mq: bool,
    pub pir: bool,
    pub ent100: bool,
    pub nab: bool,
    pub noi: bool,
    pub m5mc_strands: bool,
    pub cpg: bool,
    pub cpgnovo: bool,
}

/// Per-field FORMAT selection flags.
#[derive(Debug, Clone)]
pub struct FormatSelection {
    pub gt: bool,
    pub gl: bool,
    pub gc: bool,
    pub dp: bool,
    pub m5mc: bool,
    pub dpm5mc: bool,
    pub adm5mc: bool,
    pub ml: bool,
}

impl Default for FieldConfig {
    fn default() -> Self {
        // Defaults mirror the previous macro schema's `default` markers.
        Self {
            info: InfoSelection {
                ad: true,
                bq: true,
                dp: true,
                mq: true,
                mq0: false,
                ns: false,
                as_sb: false,
                sc5: false,
                af: false,
                abq: false,
                amq: false,
                as_ss_bq: false,
                as_ss_mq: false,
                pir: false,
                ent100: false,
                nab: false,
                noi: false,
                m5mc_strands: true,
                cpg: true,
                cpgnovo: true,
            },
            format: FormatSelection {
                gt: true,
                gl: true,
                gc: true,
                dp: true,
                m5mc: true,
                dpm5mc: true,
                adm5mc: true,
                ml: true,
            },
        }
    }
}

impl FieldConfig {
    /// Enable every field.
    pub fn with_all_fields(mut self) -> Self {
        let i = &mut self.info;
        for f in [
            &mut i.ad,
            &mut i.bq,
            &mut i.dp,
            &mut i.mq,
            &mut i.mq0,
            &mut i.ns,
            &mut i.as_sb,
            &mut i.sc5,
            &mut i.af,
            &mut i.abq,
            &mut i.amq,
            &mut i.as_ss_bq,
            &mut i.as_ss_mq,
            &mut i.pir,
            &mut i.ent100,
            &mut i.nab,
            &mut i.noi,
            &mut i.m5mc_strands,
            &mut i.cpg,
            &mut i.cpgnovo,
        ] {
            *f = true;
        }
        let m = &mut self.format;
        for f in [
            &mut m.gt,
            &mut m.gl,
            &mut m.gc,
            &mut m.dp,
            &mut m.m5mc,
            &mut m.dpm5mc,
            &mut m.adm5mc,
            &mut m.ml,
        ] {
            *f = true;
        }
        self
    }

    /// Enable the given additional INFO/FORMAT fields (by VCF ID) on top of the
    /// defaults.
    pub fn with_field_ids(
        mut self,
        info_fields: &[InfoFieldId],
        format_fields: &[FormatFieldId],
    ) -> Self {
        for id in info_fields {
            self.enable_info(id.0);
        }
        for id in format_fields {
            self.enable_format(id.0);
        }
        self
    }

    fn enable_info(&mut self, id: &str) {
        let i = &mut self.info;
        match id {
            "AD" => i.ad = true,
            "BQ" => i.bq = true,
            "DP" => i.dp = true,
            "MQ" => i.mq = true,
            "MQ0" => i.mq0 = true,
            "NS" => i.ns = true,
            "AS_SB" => i.as_sb = true,
            "SC5" => i.sc5 = true,
            "AF" => i.af = true,
            "ABQ" => i.abq = true,
            "AMQ" => i.amq = true,
            "AS_SS_BQ" => i.as_ss_bq = true,
            "AS_SS_MQ" => i.as_ss_mq = true,
            "PIR" => i.pir = true,
            "ENT100" => i.ent100 = true,
            "NAB" => i.nab = true,
            "NOI" => i.noi = true,
            "M5mC_Strands" => i.m5mc_strands = true,
            "CPG" => i.cpg = true,
            "CPGnovo" => i.cpgnovo = true,
            _ => {}
        }
    }

    fn enable_format(&mut self, id: &str) {
        let m = &mut self.format;
        match id {
            "GT" => m.gt = true,
            "GL" => m.gl = true,
            "GC" => m.gc = true,
            "DP" => m.dp = true,
            "M5mC" => m.m5mc = true,
            "DPM5mC" => m.dpm5mc = true,
            "ADM5mC" => m.adm5mc = true,
            "ML" => m.ml = true,
            _ => {}
        }
    }
}

/// A validated INFO field ID for the CLI (`--vcf-info-fields`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct InfoFieldId(pub &'static str);

impl InfoFieldId {
    pub const ALL_IDS: &'static [&'static str] = ALL_INFO_IDS;
}

impl std::str::FromStr for InfoFieldId {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ALL_INFO_IDS.iter().find(|id| **id == s).map(|id| InfoFieldId(id)).ok_or_else(|| {
            format!("Unknown INFO field: '{s}'. Available: {}", ALL_INFO_IDS.join(", "))
        })
    }
}

/// A validated FORMAT field ID for the CLI (`--vcf-format-fields`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FormatFieldId(pub &'static str);

impl FormatFieldId {
    pub const ALL_IDS: &'static [&'static str] = ALL_FORMAT_IDS;
}

impl std::str::FromStr for FormatFieldId {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ALL_FORMAT_IDS.iter().find(|id| **id == s).map(|id| FormatFieldId(id)).ok_or_else(|| {
            format!("Unknown FORMAT field: '{s}'. Available: {}", ALL_FORMAT_IDS.join(", "))
        })
    }
}
