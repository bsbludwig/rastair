//! Direct VCF/BCF emission from [`PileupMetrics`].
//!
//! This replaces the old `to_vcf_records` → intermediate `Record` → htslib
//! pipeline: per emitted line we build typed [`Alleles`] from the pileup's
//! `Base`s and drive seqair's typestate encoder straight from the metrics — no
//! intermediate record struct, no string round-trips.

use std::io::Write;

use color_eyre::{Result, eyre::Context as _, eyre::ContextCompat as _, eyre::ensure};
use seqair::vcf::{Alleles, Genotype as SeqGenotype, Ready, Writer};
use seqair_types::{Base, Phred, Pos1, Probability, SmallVec};

use crate::{
    call::{
        RecordFilters,
        pileup::indels::IndelAllele,
        variant_calling::{ErrorModel, GenotypeTag},
    },
    metrics::{AlleleMetrics, Alt, AltCall, PileupMetrics},
    utils::IntoF64 as _,
    vcf::{
        CpgBeta, CpgOrigin, InCpG, Methylated, MethylationAltDepth, MethylationDepth,
        RastairFilter,
        schema::{FieldConfig, Schema},
    },
};

/// Encode every VCF record this pileup produces into `writer`.
///
/// Mirrors the previous `to_vcf_records` + `VcfRecordSet::to_vec` selection:
/// a main record (real variants, or reference-only for CpG tracking), optional
/// rejected records (methylation evidence / read errors, only with `--vcf-all`),
/// and any indel records.
pub fn emit_pileup<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
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
            pileup.pileup.pos
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

    // Line selection, ported from `VcfRecordSet::to_vec`.
    let t = &pileup.tags;
    let cpg = t.cpg || t.denovo_cpg || t.denovo_cpg_partner;
    let (emit_main, emit_rejected) = match (record_filter.vcf_all, record_filter.cpgs_only) {
        (false, false) => (t.covered, false),
        (false, true) => (t.covered && cpg, false),
        (true, false) => (true, true),
        (true, true) => (cpg, cpg),
    };

    if emit_main {
        emit_main_record(pileup, schema, config, &real_variants, ml_threshold, error_model, writer)
            .wrap_err("Failed to emit main VCF record")?;
    }

    if emit_rejected {
        for alt in methylation_evidence.iter().chain(read_errors.iter()) {
            emit_rejected_record(pileup, schema, config, alt, ml_threshold, writer)
                .wrap_err("Failed to emit rejected VCF record")?;
        }
    }

    for call in &pileup.indel_calls {
        emit_indel_record(pileup, schema, config, call, ml_threshold, writer)
            .wrap_err("Failed to emit indel VCF record")?;
    }

    Ok(())
}

fn pos1(pileup: &PileupMetrics) -> Result<Pos1> {
    Pos1::new(pileup.pileup.pos.saturating_add(1))
        .wrap_err_with(|| format!("Invalid 1-based position from {}", pileup.pileup.pos))
}

fn to_seqair_gt(tag: GenotypeTag) -> SeqGenotype {
    let allele = |n: std::num::NonZeroU8| u16::from(n.get());
    match tag {
        GenotypeTag::HomRef => SeqGenotype::unphased(0, 0),
        GenotypeTag::RefHet(n) => SeqGenotype::unphased(0, allele(n)),
        GenotypeTag::AltHet(m, n) => SeqGenotype::unphased(allele(m), allele(n)),
        GenotypeTag::HomAlt(n) => SeqGenotype::unphased(allele(n), allele(n)),
    }
}

#[expect(clippy::too_many_arguments, reason = "self-contained per-record encode")]
fn emit_main_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    config: &FieldConfig,
    real_variants: &[&Alt],
    ml_threshold: Option<Probability>,
    error_model: &ErrorModel,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let ref_base = pileup.pileup.reference_base;
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

    let contig = schema
        .contig(&pileup.pileup.contig())
        .wrap_err_with(|| format!("Contig {} not in header", pileup.pileup.contig()))?;
    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, Some(qual as f32))?;
    let mut enc = enc.filter_pass();
    encode_info(&mut enc, schema, config, pileup, main_alts)?;
    let mut enc = enc.begin_samples();
    encode_format(&mut enc, schema, config, pileup, main_alts, ml_threshold)?;
    enc.emit()?;
    Ok(())
}

fn emit_rejected_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    config: &FieldConfig,
    alt: &Alt,
    ml_threshold: Option<Probability>,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let alleles = Alleles::snv(pileup.pileup.reference_base, alt.base)
        .wrap_err("Failed to build rejected SNV alleles")?;
    let qual = alt.filters.ml.map(|ml| Phred::from(ml.inverted()).as_int() as f32);

    // Filters: low_ml_score (if below threshold) + position + alt filters.
    let mut filters: SmallVec<RastairFilter, 8> = SmallVec::new();
    if alt.filters.ml < ml_threshold {
        filters.push(RastairFilter::LowMlScore);
    }
    filters.extend(pileup.pos_filters.iter().copied());
    filters.extend(alt.filters.filters.iter().copied());

    let contig = schema
        .contig(&pileup.pileup.contig())
        .wrap_err_with(|| format!("Contig {} not in header", pileup.pileup.contig()))?;
    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, qual)?;
    let mut enc = encode_filters(enc, schema, &filters)?;
    let alts = [alt];
    encode_info(&mut enc, schema, config, pileup, &alts)?;
    let mut enc = enc.begin_samples();
    encode_format(&mut enc, schema, config, pileup, &alts, ml_threshold)?;
    enc.emit()?;
    Ok(())
}

#[expect(clippy::cast_possible_truncation, reason = "VCF float fields are f32")]
fn emit_indel_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    config: &FieldConfig,
    call: &crate::call::variant_calling::indel_calling::IndelCall,
    ml_threshold: Option<Probability>,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let anchor = pileup.pileup.reference_base;
    let alleles = match &call.allele {
        IndelAllele::Insertion(bases) => Alleles::insertion(anchor, bases),
        IndelAllele::Deletion(bases) => Alleles::deletion(anchor, bases),
    }
    .wrap_err("Failed to build indel alleles")?;

    let qual = call
        .ml
        .map(|ml| Phred::from(ml.inverted()).as_int())
        .unwrap_or_else(|| call.quality.as_int()) as f32;

    let ml_below = call.ml.zip(ml_threshold).is_some_and(|(ml, threshold)| ml < threshold);

    let contig = schema
        .contig(&pileup.pileup.contig())
        .wrap_err_with(|| format!("Contig {} not in header", pileup.pileup.contig()))?;
    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, Some(qual))?;
    let mut enc = if ml_below {
        encode_filters(enc, schema, &[RastairFilter::LowMlScore])?
    } else {
        enc.filter_pass()
    };

    // Minimal INFO: combined depth + per-allele depth.
    if config.info.dp {
        schema.info.dp.encode(&mut enc, i32::try_from(call.depth).unwrap_or(i32::MAX));
    }
    if config.info.ad {
        let ad = [
            i32::try_from(call.depth.saturating_sub(call.alt_count)).unwrap_or(i32::MAX),
            i32::try_from(call.alt_count).unwrap_or(i32::MAX),
        ];
        schema.info.ad.encode(&mut enc, &ad);
    }

    let mut enc = enc.begin_samples();
    if config.format.gt {
        schema.format.gt.encode(&mut enc, &[to_seqair_gt(call.genotype)])?;
    }
    if config.format.dp {
        schema.format.dp.encode(&mut enc, &[i32::try_from(call.depth).unwrap_or(i32::MAX)])?;
    }
    if config.format.ml
        && let Some(ml) = call.ml
    {
        schema.format.ml.encode(&mut enc, &[&[*ml as f32]])?;
    }
    enc.emit()?;
    Ok(())
}

/// Resolve a filter set into the encoder's `Filtered` state. Empty (or all
/// unknown) → PASS; otherwise the resolved fail filters.
fn encode_filters<'a>(
    enc: seqair::vcf::RecordEncoder<'a, seqair::vcf::Begun>,
    schema: &Schema,
    filters: &[RastairFilter],
) -> Result<seqair::vcf::RecordEncoder<'a, seqair::vcf::Filtered>> {
    if filters.is_empty() {
        return Ok(enc.filter_pass());
    }
    let ids: SmallVec<&seqair::vcf::FilterId, 8> =
        filters.iter().filter_map(|f| schema.filter(f.name())).collect();
    if ids.is_empty() { Ok(enc.filter_pass()) } else { Ok(enc.filter_fail(ids.iter().copied())) }
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
        i.ad.encode(enc, &v);
    }
    if c.bq {
        i.bq.encode(enc, pm.baseq.f() as f32);
    }
    if c.dp {
        i.dp.encode(enc, i32::try_from(pm.depth).unwrap_or(i32::MAX));
    }
    if c.mq {
        i.mq.encode(enc, pm.mapq.f() as f32);
    }
    if c.mq0 {
        i.mq0.encode(enc, i32::try_from(pm.mapq0).unwrap_or(i32::MAX));
    }
    if c.ns {
        i.ns.encode(enc, 1);
    }
    if c.as_sb {
        let ot: SmallVec<i32, 3> =
            ra.iter().map(|m| i32::try_from(m.strand_count.ot).unwrap_or(i32::MAX)).collect();
        let ob: SmallVec<i32, 3> =
            ra.iter().map(|m| i32::try_from(m.strand_count.ob).unwrap_or(i32::MAX)).collect();
        i.as_sb_ot.encode(enc, &ot);
        i.as_sb_ob.encode(enc, &ob);
    }
    if c.sc5 {
        i.sc5.encode(enc, pileup.pileup.context.as_vcf_str().as_str());
    }
    if c.af {
        let v: SmallVec<f32, 2> = only_alts.iter().map(|m| m.allele_frequency.f() as f32).collect();
        // Number=A: omit entirely for reference-only sites (no ALT alleles).
        if !v.is_empty() {
            i.af.encode(enc, &v);
        }
    }
    if c.abq {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.baseq.f() as f32).collect();
        i.abq.encode(enc, &v);
    }
    if c.amq {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.mapq.f() as f32).collect();
        i.amq.encode(enc, &v);
    }
    if c.as_ss_bq {
        let ot: SmallVec<f32, 3> = ra.iter().map(|m| *m.baseq_s.ot as f32).collect();
        let ob: SmallVec<f32, 3> = ra.iter().map(|m| *m.baseq_s.ob as f32).collect();
        i.as_ss_bq_ot.encode(enc, &ot);
        i.as_ss_bq_ob.encode(enc, &ob);
    }
    if c.as_ss_mq {
        let ot: SmallVec<f32, 3> = ra.iter().map(|m| *m.mapq_s.ot as f32).collect();
        let ob: SmallVec<f32, 3> = ra.iter().map(|m| *m.mapq_s.ob as f32).collect();
        i.as_ss_mq_ot.encode(enc, &ot);
        i.as_ss_mq_ob.encode(enc, &ob);
    }
    if c.pir {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.position_in_read.f() as f32).collect();
        i.pir.encode(enc, &v);
    }
    if c.ent100 {
        i.ent100.encode(enc, pm.region_entropy as f32);
    }
    if c.nab {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.num_aligned_bases.f() as f32).collect();
        i.nab.encode(enc, &v);
    }
    if c.noi {
        let v: SmallVec<f32, 3> = ra.iter().map(|m| m.num_indels.f() as f32).collect();
        i.noi.encode(enc, &v);
    }
    if c.m5mc_strands {
        let s = pm.extended.methylation_strand_info;
        let v = [
            i32::try_from(s.unmod).unwrap_or(i32::MAX),
            i32::try_from(s.modified).unwrap_or(i32::MAX),
            i32::try_from(s.no_snp).unwrap_or(i32::MAX),
            i32::try_from(s.snp).unwrap_or(i32::MAX),
        ];
        i.m5mc_strands.encode(enc, &v);
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
    ml_threshold: Option<Probability>,
) -> Result<()> {
    let f = &schema.format;
    let c = &config.format;

    let (genotype, gl, gc) = compute_genotype(pileup, main_alts, ml_threshold);

    if c.gt {
        f.gt.encode(enc, &[genotype])?;
    }
    if c.gl {
        f.gl.encode(enc, &[gl])?;
    }
    if c.gc {
        f.gc.encode(enc, &[gc])?;
    }
    if c.dp {
        f.dp.encode(enc, &[i32::try_from(pileup.pileup.reads.len()).unwrap_or(i32::MAX)])?;
    }
    if let Some(methylated) = effective_methylation(pileup) {
        if c.m5mc {
            #[expect(clippy::cast_possible_truncation, reason = "VCF float fields are f32")]
            let values: SmallVec<f32, 2> =
                methylated.ordered_values(|b| b.beta.f() as f32).into_iter().flatten().collect();
            f.m5mc.encode(enc, &[values.as_slice()])?;
        }
        if c.dpm5mc {
            let values = counts(&MethylationDepth::from(&methylated).0);
            f.dpm5mc.encode(enc, &[values.as_slice()])?;
        }
        if c.adm5mc {
            let values = counts(&MethylationAltDepth::from(&methylated).0);
            f.adm5mc.encode(enc, &[values.as_slice()])?;
        }
    }
    if c.ml {
        let has_ml = main_alts.iter().any(|alt| alt.filters.ml.is_some());
        if has_ml {
            #[expect(clippy::cast_possible_truncation, reason = "VCF float fields are f32")]
            let values: SmallVec<f32, 2> =
                main_alts.iter().map(|alt| *alt.filters.ml.unwrap_or_default() as f32).collect();
            if !values.is_empty() {
                f.ml.encode(enc, &[values.as_slice()])?;
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

/// The genotype plus its GL/GC values, ported from the old `build_format`.
fn compute_genotype(
    pileup: &PileupMetrics,
    main_alts: &[&Alt],
    ml_threshold: Option<Probability>,
) -> (SeqGenotype, f32, f32) {
    use std::num::NonZeroU8;

    let Some(estimated) = pileup.estimate_genotype(ml_threshold, ErrorModel::default()) else {
        return (
            SeqGenotype::unphased(0, 0),
            phred_value(Probability::ZERO),
            phred_value(Probability::ZERO),
        );
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

    (to_seqair_gt(remapped), phred_value(estimated.likelihood), phred_value(estimated.confidence))
}

/// Read counts in canonical CpG order, flattened for encoding.
fn counts(values: &SmallVec<Option<u32>, 2>) -> SmallVec<i32, 2> {
    values.iter().flatten().map(|&n| i32::try_from(n).unwrap_or(i32::MAX)).collect()
}

/// The methylation values to write for this position: real measurements when
/// present, or zero-filled entries when the position sits in a CpG context that
/// produced no evidence. `None` means there is no CpG context at all and the
/// M5mC/DPM5mC/ADM5mC fields are omitted.
fn effective_methylation(pileup: &PileupMetrics) -> Option<Methylated> {
    let observed = &pileup.pos_metrics.methylated;
    if !observed.is_empty() {
        return Some(observed.clone());
    }

    // Original CpG context: the position is a reference CpG, or it is the
    // matching partner of a de-novo CpG.
    let is_orig_cpg =
        matches!(pileup.pos_metrics.cpg, InCpG::C | InCpG::G) || *pileup.pos_metrics.denovo_adj;
    // `FormsDenovo` already encodes the adjacency check, so it is not repeated
    // here. All call types count, so that rejected records written with
    // `--vcf-all` also receive a beta value.
    let is_denovo_cpg = pileup.alts.iter().any(|a| *a.metrics.denovo);

    let zero = CpgBeta {
        origin: CpgOrigin::Original,
        beta: Probability::ZERO,
        mod_count: 0,
        total_count: 0,
    };
    let mut betas = SmallVec::new();
    if is_orig_cpg {
        betas.push(zero);
    }
    if is_denovo_cpg {
        betas.push(CpgBeta { origin: CpgOrigin::DeNovo, ..zero });
    }
    (!betas.is_empty()).then(|| Methylated(betas))
}
