//! Direct VCF/BCF emission from [`PileupMetrics`].
//!
//! Per emitted line we build typed [`Alleles`] from the pileup's `Base`s and
//! drive seqair's typestate encoder straight from the metrics — no
//! intermediate record struct, no string round-trips.
//!
//! `QUAL` is Phred-scaled, from whatever decided the record: for an SNV
//! record, called or rejected, the forest's probability that its alt is not
//! real (a multi-allelic record takes its best-scored alt's); for an indel
//! record the pair-HMM's margin, a likelihood ratio and not a probability
//! (`Margin::qual`); and for a reference-only record its `GQ`.
//! A record whose decider gave no number has `QUAL` missing.
//!
//! `GQ` and `PL` come from one vector, the record's genotype likelihoods over
//! its own alleles (`pair_hmm::GenotypeLikelihoods`): `GQ` is the chance the
//! written genotype is wrong under a flat prior over all of them, `0/0`
//! included. A `0/0` from the binomial has `GQ` and no `PL`; a column with no
//! reads to weigh is `./.` with neither.
//!
//! `DP` and `AD` are the decider's reads by the same rule: an SNV record
//! counts the column's reads per base, an indel record the reads the pair-HMM
//! tells apart at its site, each for the allele its likeliest haplotype
//! carries (`HaplotypeScores::reads_over`). The two are not interchanged:
//! a window holds reads with no base at the anchor column, so an indel's `AD`
//! against the column's depth could sum to more than `DP`.

use std::{io::Write, num::NonZeroU8};

use color_eyre::{Result, eyre::Context as _};
use seqair::vcf::{Alleles, ContigId, Genotype as SeqGenotype, Ready, Writer};
use seqair_types::{Base, Phred, Pos1, SmallVec, SmolStr};

use crate::{
    call::{
        RecordFilters,
        pair_hmm::EventIdx,
        phasing::{self, PhaseCall, SiteKind},
        pileup::indels::{IndelAllele, Scored},
        variant_calling::{
            GenotypeTag, GtAllele,
            indel_calling::{BiallelicRecord, IndelCall, IndelRecords, RecordScore},
        },
    },
    metrics::{AlleleMetrics, Alt, AltVerdict, MethylationEvidenceStrandInfo, PileupMetrics},
    utils::IntoF64 as _,
    vcf::{
        RastairFilter,
        schema::{FieldConfig, Schema},
    },
};

/// Encode every VCF record a pileup the record filters keep
/// ([`RecordFilters::matches`]) produces into `writer`.
///
/// A main record (the called alts, or reference-only at a CpG), one rejected
/// record per alt that is not called (only with `--all`), and the indel
/// records (none in CpG-only output).
pub fn emit_pileup<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    record_filter: &RecordFilters,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let indices = (1..=u8::MAX).filter_map(NonZeroU8::new);
    let (real_variants, rejected): (SmallVec<IndexedAlt<'_>, 2>, SmallVec<IndexedAlt<'_>, 2>) =
        indices.zip(&pileup.alts).partition(|(_, alt)| alt.called());

    // A reference-only record (no real variants) is only emitted when it is a
    // CpG/de-novo-CpG; otherwise a covered non-CpG would carry M5mC values
    // without the CPG/CPGnovo tags set.
    if !real_variants.is_empty() || pileup.in_cpg() {
        emit_main_record(pileup, schema, contig, config, &real_variants, writer)
            .wrap_err("Failed to emit main VCF record")?;
    }

    if record_filter.vcf_all {
        for &alt in &rejected {
            emit_rejected_record(pileup, schema, contig, config, alt, writer)
                .wrap_err("Failed to emit rejected VCF record")?;
        }
    }

    // An indel anchored on a CpG's base is no CpG record.
    if record_filter.cpgs_only {
        return Ok(());
    }
    let Some(scored) = pileup.indel_data.as_ref().and_then(|d| d.scored.as_ref()) else {
        return Ok(());
    };
    match &scored.records {
        IndelRecords::Compound { calls, score } => {
            emit_compound_het_record(
                pileup,
                schema,
                contig,
                config,
                scored,
                (calls, score),
                writer,
            )
            .wrap_err("Failed to emit compound-het indel VCF record")?;
        }
        IndelRecords::Separate(records) => {
            for record in records
                .iter()
                .filter(|record| record_filter.vcf_all || record.call.genotype.carries_an_alt())
            {
                emit_indel_record(pileup, schema, contig, config, scored, record, writer)
                    .wrap_err("Failed to emit indel VCF record")?;
            }
        }
    }

    Ok(())
}

/// An alt with its 1-based index in the column's alts, as genotypes number it.
type IndexedAlt<'a> = (NonZeroU8, &'a Alt);

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
    ([first, second], score): (&[IndelCall; 2], &RecordScore),
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
        Some(score.quality.as_int() as f32),
    )?;
    let enc = enc.filter_pass();
    let events = [first.event, second.event];
    encode_indel_fields(enc, schema, config, pileup, scored, &events, (first.genotype, score))
}

/// What an indel record and a compound one share: the pair-HMM's reads over
/// `events`, and the record's genotype and score.
fn encode_indel_fields(
    mut enc: seqair::vcf::RecordEncoder<'_, seqair::vcf::Filtered>,
    schema: &Schema,
    config: &FieldConfig,
    pileup: &PileupMetrics,
    scored: &Scored,
    events: &[EventIdx],
    (genotype, score): (GenotypeTag, &RecordScore),
) -> Result<()> {
    // Reads assigned to an allele another record writes are in `DP` but in
    // no `AD` slot here.
    let reads = scored.scores.reads_over(events);
    let depth = i32::try_from(reads.depth).unwrap_or(i32::MAX);
    if config.info.dp {
        schema.info.dp.encode(&mut enc, depth)?;
    }
    if config.info.ad {
        let ad: SmallVec<i32, 3> =
            reads.ad.iter().map(|&count| i32::try_from(count).unwrap_or(i32::MAX)).collect();
        schema.info.ad.encode(&mut enc, &ad)?;
    }

    let mut enc = enc.begin_samples();
    if config.format.gt {
        encode_indel_genotype(&mut enc, schema, pileup, genotype)?;
    }
    encode_genotype_quality(&mut enc, schema, config, genotype, score)?;
    if config.format.dp {
        schema.format.dp.encode(&mut enc, &[depth])?;
    }
    enc.emit()?;
    Ok(())
}

/// An indel record's `GT`, and its `PS` where phased.
///
/// Phasing made a site of the column's one passing heterozygous indel record,
/// so the column's indel phase call belongs to whichever record that is; a
/// filtered (`0/0`) record under `--all` stays unphased, and so does a
/// homozygous one.
fn encode_indel_genotype(
    enc: &mut seqair::vcf::RecordEncoder<'_, seqair::vcf::WithSamples>,
    schema: &Schema,
    pileup: &PileupMetrics,
    tag: GenotypeTag,
) -> Result<()> {
    let phase = phasing::call_of(pileup, SiteKind::Indel).filter(|_| tag.carries_an_alt());
    let genotype = to_seqair_gt(tag, phase.map(|call| call.first));
    let phased = genotype.phased.contains(&true);
    schema.format.gt.encode(enc, &[genotype])?;
    if let Some(call) = phase.filter(|_| phased) {
        schema.format.ps.encode(enc, &[call.set.as_i32()])?;
        schema.format.pq.encode(enc, &[i32::from(call.quality)])?;
    }
    Ok(())
}

/// An indel call's `GQ` and `PL`.
fn encode_genotype_quality(
    enc: &mut seqair::vcf::RecordEncoder<'_, seqair::vcf::WithSamples>,
    schema: &Schema,
    config: &FieldConfig,
    genotype: GenotypeTag,
    score: &RecordScore,
) -> Result<()> {
    if config.format.gq {
        schema.format.gq.encode(enc, &[score.genotype_quality(genotype).as_int()])?;
    }
    if config.format.pl {
        schema.format.pl.encode_single_sample(enc, &score.likelihoods.phred())?;
    }
    Ok(())
}

fn pos1(pileup: &PileupMetrics) -> Result<Pos1> {
    pileup
        .pos
        .to_one_based()
        .wrap_err_with(|| format!("Invalid 1-based position from {}", pileup.pos))
}

/// The genotype, phased around `first` (the VCF allele haplotype 1 carries)
/// when that names one of a heterozygote's two alleles.
fn to_seqair_gt(tag: GenotypeTag, first: Option<GtAllele>) -> SeqGenotype {
    let [lower, higher] = tag.alleles();
    let number = |allele: GtAllele| u16::from(allele.number());
    match first {
        _ if lower == higher => SeqGenotype::unphased(number(lower), number(higher)),
        Some(hap1) if hap1 == lower => SeqGenotype::phased_diploid(number(lower), number(higher)),
        Some(hap1) if hap1 == higher => SeqGenotype::phased_diploid(number(higher), number(lower)),
        _ => SeqGenotype::unphased(number(lower), number(higher)),
    }
}

fn emit_main_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    real_variants: &[IndexedAlt<'_>],
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let ref_base = pileup.reference_base;
    let (indices, real_variants): (SmallVec<NonZeroU8, 2>, SmallVec<&Alt, 2>) =
        real_variants.iter().copied().unzip();
    let alleles = if real_variants.is_empty() {
        Alleles::reference(ref_base)
    } else {
        let alts: SmallVec<Base, 2> = real_variants.iter().map(|a| a.base).collect();
        Alleles::snv_multi(ref_base, &alts).wrap_err("Failed to build SNV alleles")?
    };

    let genotype = compute_genotype(pileup, &indices);
    let qual = if real_variants.is_empty() {
        genotype.quality.map(|quality| quality.as_int() as f32)
    } else {
        real_variants.iter().filter_map(|alt| snv_qual(alt)).max_by(f32::total_cmp)
    };

    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, qual)?;
    let mut enc = enc.filter_pass();
    encode_info(&mut enc, schema, config, pileup, &real_variants)?;
    let mut enc = enc.begin_samples();
    encode_format(&mut enc, schema, config, pileup, &real_variants, genotype)?;
    enc.emit()?;
    Ok(())
}

fn emit_rejected_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    (index, alt): IndexedAlt<'_>,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let alleles = Alleles::snv(pileup.reference_base, alt.base)
        .wrap_err("Failed to build rejected SNV alleles")?;
    let qual = snv_qual(alt);

    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, qual)?;
    let mut enc = match alt.verdict {
        AltVerdict::Rejected(rejection) => filter_fail(enc, schema, rejection.filter()),
        AltVerdict::Unscored | AltVerdict::Scored(_) | AltVerdict::Called(_) => enc.filter_pass(),
    };
    let alts = [alt];
    encode_info(&mut enc, schema, config, pileup, &alts)?;
    let mut enc = enc.begin_samples();
    let genotype = compute_genotype(pileup, &[index]);
    encode_format(&mut enc, schema, config, pileup, &alts, genotype)?;
    enc.emit()?;
    Ok(())
}

fn emit_indel_record<W: Write>(
    pileup: &PileupMetrics,
    schema: &Schema,
    contig: &ContigId,
    config: &FieldConfig,
    scored: &Scored,
    BiallelicRecord { call, score }: &BiallelicRecord,
    writer: &mut Writer<W, Ready>,
) -> Result<()> {
    let anchor = pileup.reference_base;
    let alleles = match &call.allele {
        IndelAllele::Insertion(bases) => Alleles::insertion(anchor, bases),
        IndelAllele::Deletion(bases) => Alleles::deletion(anchor, bases),
    }
    .wrap_err("Failed to build indel alleles")?;

    let qual = score.quality.as_int() as f32;
    let enc = writer.begin_record(contig, pos1(pileup)?, &alleles, Some(qual))?;
    // A pair-HMM call is rejected exactly when its genotype leaves the allele out.
    let enc = if call.genotype.carries_an_alt() {
        enc.filter_pass()
    } else {
        filter_fail(enc, schema, RastairFilter::LowIndelQual)
    };
    encode_indel_fields(enc, schema, config, pileup, scored, &[call.event], (call.genotype, score))
}

/// An SNV alt's `QUAL`: the forest's score, `None` where it did not score
/// the alt.
fn snv_qual(alt: &Alt) -> Option<f32> {
    Some(Phred::from(alt.verdict.ml()?.inverted()).as_int() as f32)
}

/// A record's one FILTER code.
fn filter_fail<'a>(
    enc: seqair::vcf::RecordEncoder<'a, seqair::vcf::Begun>,
    schema: &Schema,
    filter: RastairFilter,
) -> seqair::vcf::RecordEncoder<'a, seqair::vcf::Filtered> {
    enc.filter_fail([schema.filter(filter)])
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
        let s = MethylationEvidenceStrandInfo::from_pileup(pileup);
        let v = [
            i32::try_from(s.unmod).unwrap_or(i32::MAX),
            i32::try_from(s.modified).unwrap_or(i32::MAX),
            i32::try_from(s.no_snp).unwrap_or(i32::MAX),
            i32::try_from(s.snp).unwrap_or(i32::MAX),
        ];
        i.m5mc_strands.encode(enc, &v)?;
    }
    if c.cpg && pileup.reference_cpg().is_some() {
        i.cpg.encode(enc);
    }
    if c.cpgnovo && pileup.in_denovo_cpg() {
        i.cpgnovo.encode(enc);
    }
    Ok(())
}

fn encode_format(
    enc: &mut seqair::vcf::RecordEncoder<'_, seqair::vcf::WithSamples>,
    schema: &Schema,
    config: &FieldConfig,
    pileup: &PileupMetrics,
    main_alts: &[&Alt],
    genotype: GenotypeFields,
) -> Result<()> {
    let f = &schema.format;
    let c = &config.format;

    if c.gt {
        f.gt.encode(enc, &[genotype.genotype])?;
        // An unphased record omits `PS` and `PQ` rather than carrying empty values.
        if let Some(call) = genotype.phase {
            f.ps.encode(enc, &[call.set.as_i32()])?;
            f.pq.encode(enc, &[i32::from(call.quality)])?;
        }
    }
    if c.gq
        && let Some(quality) = genotype.quality
    {
        f.gq.encode(enc, &[quality.as_int()])?;
    }
    if c.pl
        && let Some(pl) = &genotype.pl
    {
        f.pl.encode_single_sample(enc, pl)?;
    }
    if c.dp {
        f.dp.encode(enc, &[i32::try_from(pileup.pos_metrics.depth).unwrap_or(i32::MAX)])?;
    }
    // Only emit the methylation fields when the position actually has a CpG
    // context: htslib complains when it reads a zero-length FORMAT field.
    let written = crate::metrics::methylation::written(pileup);
    if !written.is_empty() {
        let cpgs = written.iter().map(|(_, cpg)| cpg.as_ref());
        let count = |n: u32| i32::try_from(n).unwrap_or(i32::MAX);
        if c.m5mc {
            #[expect(clippy::cast_possible_truncation, reason = "VCF float fields are f32")]
            let values: SmallVec<Option<f32>, 2> =
                cpgs.clone().map(|cpg| cpg.map(|cpg| cpg.beta.f() as f32)).collect();
            f.m5mc.encode_single_sample(enc, &values)?;
        }
        if c.dpm5mc {
            let values: SmallVec<i32, 2> =
                cpgs.clone().map(|cpg| cpg.map_or(0, |cpg| count(cpg.total_count.get()))).collect();
            f.dpm5mc.encode_single_sample(enc, &values)?;
        }
        if c.adm5mc {
            let values: SmallVec<i32, 2> =
                cpgs.map(|cpg| cpg.map_or(0, |cpg| count(cpg.mod_count))).collect();
            f.adm5mc.encode_single_sample(enc, &values)?;
        }
    }
    // A record writes called alts, or one rejected alt: every alt has a
    // score, or the record's one alt has none.
    if c.ml {
        #[expect(clippy::cast_possible_truncation, reason = "VCF float fields are f32")]
        let values: Option<SmallVec<f32, 2>> =
            main_alts.iter().map(|alt| Some(*alt.verdict.ml()? as f32)).collect();
        if let Some(values) = values.filter(|values| !values.is_empty()) {
            f.ml.encode_single_sample(enc, &values)?;
        }
    }
    Ok(())
}

/// A record's genotype fields.
struct GenotypeFields {
    genotype: SeqGenotype,
    /// `GQ`; `None` where `genotype` is missing.
    quality: Option<Phred>,
    /// `PL`, in the record's ALT order.
    pl: Option<SmallVec<i32, 6>>,
    /// Set exactly where `genotype` is phased.
    phase: Option<PhaseCall>,
}

/// The genotype plus its `GQ`, `PL` and phase call of a record writing the
/// column's alts `vcf_alts`, in that order. A record without an estimate, a
/// column with no reads to weigh, is `./.`.
///
/// The estimate is the one the calling pipeline stored, not a fresh one: an
/// estimate re-derived here would silently ignore `--error-model` and so
/// contradict the BED output and the methylation calls, which read the stored
/// value. The phase call names alleles the same way, so one remap serves both.
fn compute_genotype(pileup: &PileupMetrics, vcf_alts: &[NonZeroU8]) -> GenotypeFields {
    let Some(estimated) = &pileup.pos_metrics.genotype else {
        return GenotypeFields {
            genotype: SeqGenotype::missing_diploid(),
            quality: None,
            pl: None,
            phase: None,
        };
    };

    // A genotype allele as the record's VCF allele; `None` for an alt the
    // record does not write.
    let vcf = |allele: GtAllele| match allele {
        GtAllele::Ref => Some(GtAllele::Ref),
        GtAllele::Alt(n) => {
            let vcf_idx = vcf_alts.iter().position(|&alt| alt == n)? + 1;
            Some(GtAllele::Alt(NonZeroU8::new(u8::try_from(vcf_idx).ok()?)?))
        }
    };

    // An alt the record does not write is written as the reference.
    let [a, b] = estimated.genotype.alleles().map(|allele| vcf(allele).unwrap_or(GtAllele::Ref));
    let remapped = GenotypeTag::from_alleles(a, b);

    // A phase call whose allele does not survive the remap cannot be written
    // as a genotype either, and the record stays unphased.
    let phase = phasing::call_of(pileup, SiteKind::Snv);
    let genotype = to_seqair_gt(remapped, phase.and_then(|call| vcf(call.first)));
    let phased = genotype.phased.contains(&true);

    let pl = estimated
        .source
        .likelihoods()
        .and_then(|likelihoods| likelihoods.over(vcf_alts))
        .map(|likelihoods| likelihoods.phred());

    // A genotype that names an alt this record does not write (a filtered
    // record beside a called alt) is written without its alt, and its `GQ`
    // does not describe what is written.
    let written_as_estimated =
        estimated.genotype.alleles().into_iter().all(|allele| vcf(allele).is_some());

    GenotypeFields {
        genotype,
        quality: written_as_estimated.then_some(estimated.quality),
        pl,
        phase: phase.filter(|_| phased),
    }
}
