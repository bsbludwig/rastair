//! Direct VCF/BCF emission from [`PileupMetrics`].
//!
//! Per emitted line we build typed [`Alleles`] from the pileup's `Base`s and
//! drive seqair's typestate encoder straight from the metrics — no
//! intermediate record struct, no string round-trips.

use std::{
    io::Write,
    num::{NonZeroU8, NonZeroU16},
};

use color_eyre::{Result, eyre::Context as _, eyre::ContextCompat as _, eyre::ensure};
use enumset::EnumSet;
use seqair::vcf::{Alleles, ContigId, Genotype as SeqGenotype, Ready, Writer};
use seqair_types::{Base, Phred, Pos1, Probability, SmallVec, SmolStr};

use crate::{
    call::{
        RecordFilters,
        pair_hmm::EventIdx,
        phasing::{self, GtAllele, PhaseCall, SiteKind},
        pileup::indels::{IndelAllele, Scored},
        variant_calling::{
            ErrorModel, GenotypeTag,
            indel_calling::{IndelCall, IndelRecords},
        },
    },
    metrics::{AlleleMetrics, Alt, AltCall, PileupMetrics},
    utils::IntoF64 as _,
    vcf::{
        CpgBeta, InCpG, Methylated, MethylationAltDepth, MethylationDepth, RastairFilter,
        schema::{FieldConfig, Schema},
    },
};

/// Encode every VCF record this pileup produces into `writer`.
///
/// Mirrors the previous `to_vcf_records` + `VcfRecordSet::to_vec` selection:
/// a main record (real variants, or reference-only for CpG tracking), optional
/// rejected records (methylation evidence / read errors, only with `--vcf-all`),
/// and any indel records.
#[expect(clippy::too_many_arguments, reason = "self-contained per-pileup encode entry point")]
pub fn emit_pileup<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    ml_threshold: Option<Probability>,
    error_model: &ErrorModel,
    record_filter: &RecordFilters,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    for alt in &pileup.alts {
        ensure!(
            !matches!(alt.call, AltCall::Uncalled),
            "Alt {} at position {} is Uncalled - this should not happen",
            alt.base,
            pileup.pos
        );
    }

    let mut real_variants: SmallVec<&Alt, 2> = SmallVec::new();
    let mut methylation_evidence: SmallVec<&Alt, 2> = SmallVec::new();
    let mut read_errors: SmallVec<&Alt, 2> = SmallVec::new();
    for alt in &pileup.alts {
        match alt.call {
            AltCall::RealVariant => real_variants.push(alt),
            AltCall::MethylationEvidenceOnly { .. } => methylation_evidence.push(alt),
            AltCall::ReadError => read_errors.push(alt),
            AltCall::Uncalled => unreachable!("checked above"),
        }
    }

    // Line selection.
    let t = &pileup.tags;
    let cpg = t.cpg || t.denovo_cpg || t.denovo_cpg_partner;
    // A reference-only record (no real variants) is only emitted when it is a
    // CpG/de-novo-CpG; otherwise a covered non-CpG would carry M5mC values
    // without the CPG/CPGnovo tags set.
    let main_is_ref_only = real_variants.is_empty();
    let want_main = (!main_is_ref_only || cpg)
        && match (record_filter.vcf_all, record_filter.cpgs_only) {
            (false, false) => t.covered,
            (false, true) => t.covered && cpg,
            (true, false) => true,
            (true, true) => cpg,
        };
    let emit_rejected = match (record_filter.vcf_all, record_filter.cpgs_only) {
        (false, _) => false,
        (true, false) => true,
        (true, true) => cpg,
    };

    if want_main {
        emit_main_record(pileup, schema, contig, config, &real_variants, error_model, writer)
            .wrap_err("Failed to emit main VCF record")?;
    }

    if emit_rejected {
        for alt in methylation_evidence.iter().chain(read_errors.iter()) {
            emit_rejected_record(pileup, schema, contig, config, alt, ml_threshold, writer)
                .wrap_err("Failed to emit rejected VCF record")?;
        }
    }

    let Some(scored) = pileup.indel_data.as_ref().and_then(|d| d.scored.as_ref()) else {
        return Ok(());
    };
    match &scored.records {
        IndelRecords::Compound(calls) => {
            emit_compound_het_record(pileup, schema, contig, config, scored, calls, writer)
                .wrap_err("Failed to emit compound-het indel VCF record")?;
        }
        IndelRecords::Separate(calls) => {
            for call in calls.iter().filter(|call| record_filter.vcf_all || call.is_called()) {
                emit_indel_record(pileup, schema, contig, config, scored, call, writer)
                    .wrap_err("Failed to emit indel VCF record")?;
            }
        }
    }

    Ok(())
}

/// Emit a compound heterozygote as one multi-allelic record.
///
/// Both alleles share an anchor, so REF spans the longest deletion among them
/// and each ALT is written against that span.
fn emit_compound_het_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    scored: &Scored,
    [first, second]: &[IndelCall; 2],
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let anchor = pileup.reference_base;
    let render = |bases: &[Base]| -> String { bases.iter().map(|b| b.as_str()).collect() };
    let deleted: &[Base] = [first, second]
        .into_iter()
        .filter_map(|c| match &c.allele {
            IndelAllele::Deletion(bases) => Some(&bases[..]),
            IndelAllele::Insertion(_) => None,
        })
        .max_by_key(|b| b.len())
        .unwrap_or(&[]);
    let alt_of = |call: &IndelCall| -> SmolStr {
        match &call.allele {
            IndelAllele::Deletion(bases) => {
                format!("{anchor}{}", render(deleted.get(bases.len()..).unwrap_or(&[]))).into()
            }
            IndelAllele::Insertion(bases) => {
                format!("{anchor}{}{}", render(bases), render(deleted)).into()
            }
        }
    };
    let mut alts: SmallVec<SmolStr, 2> = SmallVec::new();
    alts.push(alt_of(first));
    alts.push(alt_of(second));
    let alleles = Alleles::complex(format!("{anchor}{}", render(deleted)).into(), alts);

    let enc = writer.begin_record(
        contig,
        pos1(pileup)?,
        &alleles,
        Some(first.quality.as_int() as f32),
    )?;
    let mut enc = enc.filter_pass();

    let depth = scored.haplotypes.depth();
    if config.info.dp {
        schema.info.dp.encode(&mut enc, i32::try_from(depth).unwrap_or(i32::MAX))?;
    }
    if config.info.ad {
        schema.info.ad.encode(&mut enc, &indel_ad(scored, &[first.event, second.event]))?;
    }

    let mut enc = enc.begin_samples();
    if config.format.gt {
        encode_indel_genotype(&mut enc, schema, pileup, first.genotype, true)?;
    }
    encode_genotype_quality(&mut enc, schema, config, scored)?;
    if config.format.dp {
        schema.format.dp.encode(&mut enc, &[i32::try_from(depth).unwrap_or(i32::MAX)])?;
    }
    enc.emit()?;
    Ok(())
}

/// An indel record's `GT`, and its `PS` where phased.
///
/// Phasing made a site of the column's one passing heterozygous indel record,
/// so the column's indel phase call belongs to whichever record that is; a
/// filtered record under `--all` stays unphased, and so does a homozygous one.
fn encode_indel_genotype(
    enc: &mut seqair::vcf::RecordEncoder<'_, seqair::vcf::WithSamples>,
    schema: &Schema,
    pileup: &PileupMetrics,
    tag: GenotypeTag,
    passes: bool,
) -> Result<()> {
    let phase = phasing::call_of(pileup, SiteKind::Indel).filter(|_| passes);
    let first = phase.map(|call| match call.first {
        GtAllele::Ref => 0,
        GtAllele::Alt(n) => NonZeroU16::from(n).get(),
    });
    let genotype = to_seqair_gt(tag, first);
    let phased = genotype.phased.contains(&true);
    schema.format.gt.encode(enc, &[genotype])?;
    if let Some(call) = phase.filter(|_| phased) {
        schema.format.ps.encode(enc, &[call.set.as_i32()])?;
        schema.format.pq.encode(enc, &[i32::from(call.quality)])?;
    }
    Ok(())
}

/// The pair-HMM's genotype quality as `GC`, the field SNV records carry
/// theirs in.
fn encode_genotype_quality(
    enc: &mut seqair::vcf::RecordEncoder<'_, seqair::vcf::WithSamples>,
    schema: &Schema,
    config: &FieldConfig,
    scored: &Scored,
) -> Result<()> {
    if config.format.gc {
        schema.format.gc.encode(enc, &[scored.haplotypes.verdict.quality.as_int() as f32])?;
    }
    Ok(())
}

fn pos1(pileup: &PileupMetrics) -> Result<Pos1> {
    Pos1::new(pileup.pos.saturating_add(1))
        .wrap_err_with(|| format!("Invalid 1-based position from {}", pileup.pos))
}

/// The genotype, phased around `first` (the VCF allele haplotype 1 carries)
/// when that names one of a heterozygote's two alleles.
fn to_seqair_gt(tag: GenotypeTag, first: Option<u16>) -> SeqGenotype {
    let allele = |n: NonZeroU8| NonZeroU16::from(n).get();
    let (lower, higher) = match tag {
        GenotypeTag::HomRef => return SeqGenotype::unphased(0, 0),
        GenotypeTag::HomAlt(n) => return SeqGenotype::unphased(allele(n), allele(n)),
        GenotypeTag::RefHet(n) => (0, allele(n)),
        GenotypeTag::AltHet(m, n) => (allele(m), allele(n)),
    };
    match first {
        Some(hap1) if hap1 == lower => SeqGenotype::phased_diploid(lower, higher),
        Some(hap1) if hap1 == higher => SeqGenotype::phased_diploid(higher, lower),
        _ => SeqGenotype::unphased(lower, higher),
    }
}

fn emit_main_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    real_variants: &[&Alt],
    error_model: &ErrorModel,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let ref_base = pileup.reference_base;
    let alleles = if real_variants.is_empty() {
        Alleles::reference(ref_base)
    } else {
        let alts: SmallVec<Base, 2> = real_variants.iter().map(|a| a.base).collect();
        Alleles::snv_multi(ref_base, &alts).wrap_err("Failed to build SNV alleles")?
    };

    let depth = pileup.pos_metrics.depth.max(1).f();
    let ml_qual: Option<i32> = if real_variants.is_empty() {
        Some(Phred::from_phred(99_u8).as_int())
    } else {
        real_variants
            .iter()
            .filter_map(|alt| alt.filters.ml)
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|p| Phred::from(p.inverted()).as_int())
    };
    let qual = match ml_qual {
        Some(q) => q,
        None => Phred::from(
            Probability::new(*error_model.error_rate() / depth)
                .wrap_err("Failed to compute QUAL from error rate")?,
        )
        .as_int(),
    };

    let main_alts: &[&Alt] = if real_variants.is_empty() { &[] } else { real_variants };

    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, Some(qual as f32))?;
    let mut enc = enc.filter_pass();
    encode_info(&mut enc, schema, config, pileup, main_alts)?;
    let mut enc = enc.begin_samples();
    encode_format(&mut enc, schema, config, pileup, main_alts)?;
    enc.emit()?;
    Ok(())
}

fn emit_rejected_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    alt: &Alt,
    ml_threshold: Option<Probability>,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let alleles = Alleles::snv(pileup.reference_base, alt.base)
        .wrap_err("Failed to build rejected SNV alleles")?;
    let qual = alt.filters.ml.map(|ml| Phred::from(ml.inverted()).as_int() as f32);

    // Filters: low_ml_score (if below threshold) + position + alt filters.
    let mut filters = pileup.pos_filters.as_set() | alt.filters.filters.as_set();
    if alt.filters.ml < ml_threshold {
        filters.insert(RastairFilter::LowMlScore);
    }

    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, qual)?;
    let mut enc = encode_filters(enc, schema, filters)?;
    let alts = [alt];
    encode_info(&mut enc, schema, config, pileup, &alts)?;
    let mut enc = enc.begin_samples();
    encode_format(&mut enc, schema, config, pileup, &alts)?;
    enc.emit()?;
    Ok(())
}

fn emit_indel_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    scored: &Scored,
    call: &IndelCall,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let anchor = pileup.reference_base;
    let alleles = match &call.allele {
        IndelAllele::Insertion(bases) => Alleles::insertion(anchor, bases),
        IndelAllele::Deletion(bases) => Alleles::deletion(anchor, bases),
    }
    .wrap_err("Failed to build indel alleles")?;

    let qual = call.quality.as_int() as f32;
    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, Some(qual))?;
    // A pair-HMM call is rejected exactly when its genotype leaves the allele out.
    let mut enc = if call.is_called() {
        enc.filter_pass()
    } else {
        encode_filters(enc, schema, RastairFilter::LowIndelQual.into())?
    };

    let depth = scored.haplotypes.depth();
    if config.info.dp {
        schema.info.dp.encode(&mut enc, i32::try_from(depth).unwrap_or(i32::MAX))?;
    }
    if config.info.ad {
        schema.info.ad.encode(&mut enc, &indel_ad(scored, &[call.event]))?;
    }

    let mut enc = enc.begin_samples();
    if config.format.gt {
        encode_indel_genotype(&mut enc, schema, pileup, call.genotype, call.is_called())?;
    }
    encode_genotype_quality(&mut enc, schema, config, scored)?;
    if config.format.dp {
        schema.format.dp.encode(&mut enc, &[i32::try_from(depth).unwrap_or(i32::MAX)])?;
    }
    enc.emit()?;
    Ok(())
}

/// An indel record's `AD`: the reads the pair-HMM assigns to the reference and
/// to each of the record's alleles. Reads assigned to an allele another record
/// writes are in `DP` but in no `AD` slot here.
fn indel_ad(scored: &Scored, events: &[EventIdx]) -> SmallVec<i32, 3> {
    scored
        .haplotypes
        .support_over(events)
        .into_iter()
        .map(|count| i32::try_from(count).unwrap_or(i32::MAX))
        .collect()
}

/// Resolve a filter set into the encoder's `Filtered` state. Empty (or all
/// unknown) → PASS; otherwise the resolved fail filters.
fn encode_filters<'a>(
    enc: seqair::vcf::RecordEncoder<'a, seqair::vcf::Begun>,
    schema: &Schema,
    filters: EnumSet<RastairFilter>,
) -> Result<seqair::vcf::RecordEncoder<'a, seqair::vcf::Filtered>> {
    if filters.is_empty() {
        return Ok(enc.filter_pass());
    }
    Ok(enc.filter_fail(filters.iter().map(|f| schema.filter(f))))
}

/// Per-allele metric refs: reference allele first, then the given alts.
fn ref_alts<'a>(pileup: &'a PileupMetrics, alts: &[&'a Alt]) -> SmallVec<&'a AlleleMetrics, 3> {
    let mut xs: SmallVec<&AlleleMetrics, 3> = SmallVec::new();
    xs.push(&pileup.ref_metrics);
    for alt in alts {
        xs.push(&alt.metrics);
    }
    xs
}

#[expect(clippy::cast_possible_truncation, reason = "VCF float fields are f32")]
fn encode_info(
    enc: &mut seqair::vcf::RecordEncoder<'_, seqair::vcf::Filtered>,
    schema: &Schema,
    config: &FieldConfig,
    pileup: &PileupMetrics,
    alts: &[&Alt],
) -> Result<()> {
    let i = &schema.info;
    let c = &config.info;
    let pm = &pileup.pos_metrics;
    let ra = ref_alts(pileup, alts);
    let only_alts = ra.get(1..).unwrap_or(&[]);

    if c.ad {
        let v: SmallVec<i32, 3> =
            ra.iter().map(|m| i32::try_from(m.depth).unwrap_or(i32::MAX)).collect();
        i.ad.encode(enc, &v)?;
    }
    if c.bq {
        i.bq.encode(enc, pm.baseq.f() as f32)?;
    }
    if c.dp {
        i.dp.encode(enc, i32::try_from(pm.depth).unwrap_or(i32::MAX))?;
    }
    if c.mq {
        i.mq.encode(enc, pm.mapq.f() as f32)?;
    }
    if c.mq0 {
        i.mq0.encode(enc, i32::try_from(pm.mapq0).unwrap_or(i32::MAX))?;
    }
    if c.ns {
        i.ns.encode(enc, 1)?;
    }
    if c.as_sb {
        let ot: SmallVec<i32, 3> =
            ra.iter().map(|m| i32::try_from(m.strand_count.ot).unwrap_or(i32::MAX)).collect();
        let ob: SmallVec<i32, 3> =
            ra.iter().map(|m| i32::try_from(m.strand_count.ob).unwrap_or(i32::MAX)).collect();
        i.as_sb_ot.encode(enc, &ot)?;
        i.as_sb_ob.encode(enc, &ob)?;
    }
    if c.sc5 {
        i.sc5.encode(enc, pileup.context.as_vcf_str().as_str());
    }
    if c.af {
        let v: SmallVec<f32, 2> = only_alts.iter().map(|m| m.allele_frequency.f() as f32).collect();
        // Number=A: omit entirely for reference-only sites (no ALT alleles).
        if !v.is_empty() {
            i.af.encode(enc, &v)?;
        }
    }
    if c.abq {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.baseq.f() as f32).collect();
        i.abq.encode(enc, &v)?;
    }
    if c.amq {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.mapq.f() as f32).collect();
        i.amq.encode(enc, &v)?;
    }
    if c.as_ss_bq {
        let ot: SmallVec<f32, 3> = ra.iter().map(|m| *m.baseq_s.ot as f32).collect();
        let ob: SmallVec<f32, 3> = ra.iter().map(|m| *m.baseq_s.ob as f32).collect();
        i.as_ss_bq_ot.encode(enc, &ot)?;
        i.as_ss_bq_ob.encode(enc, &ob)?;
    }
    if c.as_ss_mq {
        let ot: SmallVec<f32, 3> = ra.iter().map(|m| *m.mapq_s.ot as f32).collect();
        let ob: SmallVec<f32, 3> = ra.iter().map(|m| *m.mapq_s.ob as f32).collect();
        i.as_ss_mq_ot.encode(enc, &ot)?;
        i.as_ss_mq_ob.encode(enc, &ob)?;
    }
    if c.pir {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.position_in_read.f() as f32).collect();
        i.pir.encode(enc, &v)?;
    }
    if c.ent100 {
        i.ent100.encode(enc, pm.region_entropy as f32)?;
    }
    if c.nab {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.num_aligned_bases.f() as f32).collect();
        i.nab.encode(enc, &v)?;
    }
    if c.noi {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.num_indels.f() as f32).collect();
        i.noi.encode(enc, &v)?;
    }
    if c.m5mc_strands {
        let s = pm.extended.methylation_strand_info;
        let v = [
            i32::try_from(s.unmod).unwrap_or(i32::MAX),
            i32::try_from(s.modified).unwrap_or(i32::MAX),
            i32::try_from(s.no_snp).unwrap_or(i32::MAX),
            i32::try_from(s.snp).unwrap_or(i32::MAX),
        ];
        i.m5mc_strands.encode(enc, &v)?;
    }
    if c.cpg && pm.cpg != InCpG::No {
        i.cpg.encode(enc);
    }
    if c.cpgnovo {
        let tags = &pileup.tags;
        if tags.denovo_cpg || tags.denovo_cpg_partner {
            i.cpgnovo.encode(enc);
        }
    }
    Ok(())
}

fn encode_format(
    enc: &mut seqair::vcf::RecordEncoder<'_, seqair::vcf::WithSamples>,
    schema: &Schema,
    config: &FieldConfig,
    pileup: &PileupMetrics,
    main_alts: &[&Alt],
) -> Result<()> {
    let f = &schema.format;
    let c = &config.format;

    let genotype = compute_genotype(pileup, main_alts);

    if c.gt {
        f.gt.encode(enc, &[genotype.genotype])?;
        // An unphased record omits `PS` and `PQ` rather than carrying empty values.
        if let Some(call) = genotype.phase {
            f.ps.encode(enc, &[call.set.as_i32()])?;
            f.pq.encode(enc, &[i32::from(call.quality)])?;
        }
    }
    if c.gl {
        f.gl.encode(enc, &[genotype.likelihood])?;
    }
    if c.gc {
        f.gc.encode(enc, &[genotype.confidence])?;
    }
    if c.dp {
        f.dp.encode(enc, &[i32::try_from(pileup.pos_metrics.depth).unwrap_or(i32::MAX)])?;
    }
    // Only emit the methylation fields when the position actually has a CpG
    // context: htslib complains when it reads a zero-length FORMAT field.
    let methylated = effective_methylation(pileup);
    if c.m5mc {
        #[expect(clippy::cast_possible_truncation, reason = "VCF float fields are f32")]
        let values: SmallVec<f32, 2> =
            methylated.ordered_values(|b| b.beta.f() as f32).into_iter().flatten().collect();
        if !values.is_empty() {
            f.m5mc.encode_single_sample(enc, &values)?;
        }
    }
    if c.dpm5mc {
        let values = counts(&MethylationDepth::from(&methylated).0);
        if !values.is_empty() {
            f.dpm5mc.encode_single_sample(enc, &values)?;
        }
    }
    if c.adm5mc {
        let values = counts(&MethylationAltDepth::from(&methylated).0);
        if !values.is_empty() {
            f.adm5mc.encode_single_sample(enc, &values)?;
        }
    }
    if c.ml {
        let has_ml = main_alts.iter().any(|alt| alt.filters.ml.is_some());
        if has_ml {
            #[expect(clippy::cast_possible_truncation, reason = "VCF float fields are f32")]
            let values: SmallVec<f32, 2> =
                main_alts.iter().map(|alt| *alt.filters.ml.unwrap_or_default() as f32).collect();
            if !values.is_empty() {
                f.ml.encode_single_sample(enc, &values)?;
            }
        }
    }
    Ok(())
}

/// Phred-scaled GL/GC value as a single float (matches the previous output of
/// one value per record).
fn phred_value(p: Probability) -> f32 {
    Phred::from(p).as_int() as f32
}

/// A record's genotype fields.
struct GenotypeFields {
    genotype: SeqGenotype,
    likelihood: f32,
    confidence: f32,
    /// Set exactly where `genotype` is phased.
    phase: Option<PhaseCall>,
}

/// The genotype plus its GL/GC values and phase call.
///
/// The estimate is the one the calling pipeline stored, not a fresh one: an
/// estimate re-derived here would silently ignore `--error-model` and so
/// contradict the BED output and the methylation calls, which read the stored
/// value. The phase call names alleles the same way, so one remap serves both.
fn compute_genotype(pileup: &PileupMetrics, main_alts: &[&Alt]) -> GenotypeFields {
    let Some(estimated) = pileup.pos_metrics.extended.genotype else {
        return GenotypeFields {
            genotype: SeqGenotype::unphased(0, 0),
            likelihood: phred_value(Probability::ZERO),
            confidence: phred_value(Probability::ZERO),
            phase: None,
        };
    };

    // Remap genotype allele indices from `pileup.alts` positions to the VCF ALT
    // order (only real variants appear in the main record).
    let mut self_to_vcf: SmallVec<Option<usize>, 2> = SmallVec::new();
    self_to_vcf.resize(pileup.alts.len(), None);
    for (vcf_idx, main_alt) in main_alts.iter().enumerate() {
        for (self_idx, self_alt) in pileup.alts.iter().enumerate() {
            if std::ptr::eq(*main_alt, self_alt) {
                if let Some(slot) = self_to_vcf.get_mut(self_idx) {
                    *slot = Some(vcf_idx + 1);
                }
                break;
            }
        }
    }
    let map = |n: NonZeroU8| -> Option<NonZeroU8> {
        let self_idx = usize::from(n.get()).checked_sub(1)?;
        let vcf_idx = self_to_vcf.get(self_idx).and_then(|x| *x)?;
        NonZeroU8::new(u8::try_from(vcf_idx).unwrap_or(n.get()))
    };

    let remapped = match estimated.genotype {
        GenotypeTag::HomRef => GenotypeTag::HomRef,
        GenotypeTag::RefHet(n) => map(n).map_or(GenotypeTag::HomRef, GenotypeTag::ref_het),
        GenotypeTag::HomAlt(n) => map(n).map_or(GenotypeTag::HomRef, GenotypeTag::hom_alt),
        GenotypeTag::AltHet(m, n) => match (map(m), map(n)) {
            (Some(vm), Some(vn)) => GenotypeTag::alt_het(vm, vn),
            (Some(v), None) | (None, Some(v)) => GenotypeTag::ref_het(v),
            (None, None) => GenotypeTag::HomRef,
        },
    };

    // A phase call whose allele does not survive the remap cannot be written
    // as a genotype either, and the record stays unphased.
    let phase = phasing::call_of(pileup, SiteKind::Snv);
    let first = phase.and_then(|call| match call.first {
        GtAllele::Ref => Some(0),
        GtAllele::Alt(n) => map(n).map(|vcf| NonZeroU16::from(vcf).get()),
    });
    let genotype = to_seqair_gt(remapped, first);
    let phased = genotype.phased.contains(&true);

    GenotypeFields {
        genotype,
        likelihood: phred_value(estimated.likelihood),
        confidence: phred_value(estimated.confidence),
        phase: phase.filter(|_| phased),
    }
}

/// Read counts in canonical CpG order, flattened for encoding.
fn counts(values: &SmallVec<Option<u32>, 2>) -> SmallVec<i32, 2> {
    values.iter().flatten().map(|&n| i32::try_from(n).unwrap_or(i32::MAX)).collect()
}

/// The methylation values to write for this position: real measurements when
/// present, or zero-filled entries when the position sits in a CpG context that
/// produced no evidence. Empty means there is no CpG context at all, and the
/// M5mC/DPM5mC/ADM5mC fields render as the missing value `.`.
///
/// `origins()` only reports a de-novo partner as a CpG once its partner's alt
/// was actually called (`other_pos_in_denovo_passes`), so a rejected de-novo
/// candidate never yields a zero `M5mC` without the `CPG`/`CPGnovo` tags set.
fn effective_methylation(pileup: &PileupMetrics) -> Methylated {
    let observed = &pileup.pos_metrics.methylated;
    if !observed.is_empty() {
        return observed.clone();
    }

    Methylated(
        crate::metrics::methylation::origins(pileup)
            .into_iter()
            .map(|origin| CpgBeta { origin, beta: Probability::ZERO, mod_count: 0, total_count: 0 })
            .collect(),
    )
}
