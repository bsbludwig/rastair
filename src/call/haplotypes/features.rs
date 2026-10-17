use super::{
    ExpansionParams, Fragment, Indel, LikelihoodParams, Reference, Site, likelihood::log_add,
};
use crate::metrics::AlleleMetrics;
use crate::utils::IntoF32 as _;
use crate::{
    call::{
        phasing::{ColumnPhase, PhaseMask},
        pileup::indels::IndelAllele,
    },
    metrics::{MetricsForAlt, PileupMetrics, ml::features::HaplotypeEvidence},
    sequence::Segment,
};
use color_eyre::eyre::{OptionExt as _, Result};
use rayon::prelude::*;
use seqair_types::{Base, Strand};
use std::collections::{BTreeMap, BTreeSet};

/// A column's molecular SNV evidence: the forest's features for each
/// alternate base and the diploid posterior its genotype comes from. Indels
/// are the pair-HMM indel caller's; here they are only context.
#[derive(Debug, Clone, Default)]
pub struct ColumnEvidence {
    snps: [HaplotypeEvidence; 4],
    pub(super) snp_genotype: Option<super::genotype::DiploidPosterior<Base>>,
}

impl ColumnEvidence {
    pub(crate) fn snp(current: &MetricsForAlt) -> HaplotypeEvidence {
        current
            .metrics
            .phase
            .as_deref()
            .and_then(|phase| match phase {
                ColumnPhase::Likelihood(evidence) => current
                    .alt
                    .base
                    .known_index()
                    .and_then(|index| evidence.snps.get(index))
                    .copied(),
                _ => None,
            })
            .unwrap_or_default()
    }
}

pub(crate) fn parameters() -> Result<LikelihoodParams> {
    use seqair_types::Probability;
    Ok(LikelihoodParams {
        error_profile: super::ErrorProfile::default(),
        likelihood_engine: super::LikelihoodEngine::PairHmm,
        conversion_rate: Probability::new(0.95)?,
        false_conversion_rate: Probability::new(0.003)?,
        non_taps: false,
        methylation_prior_strength: 4.0,
        methylation_min_cpg_probability: Probability::new(0.9)?,
        quadrature_order: 16,
        snp_prior: Probability::new(0.001)?,
        indel_prior: Probability::new(0.0001)?,
        gap_open: Probability::new(10.0_f64.powf(-4.5))?,
        gap_extend: Probability::new(0.1)?,
        realignment_flank: 12,
        max_indel_length: 100,
    })
}

pub(crate) fn add_likelihoods(
    columns: &mut [PileupMetrics],
    segment: &Segment,
    fragments: &[Fragment],
    mask: &PhaseMask,
) -> Result<()> {
    add_likelihoods_with_profile(columns, segment, fragments, mask, &super::ErrorProfile::default())
}

pub(crate) fn add_likelihoods_with_profile(
    columns: &mut [PileupMetrics],
    segment: &Segment,
    fragments: &[Fragment],
    mask: &PhaseMask,
    profile: &super::ErrorProfile,
) -> Result<()> {
    profile.validate()?;
    let mut params = parameters()?;
    params.error_profile = profile.clone();
    let search = ExpansionParams { read_local: true, ..ExpansionParams::default() };
    let reference_end = segment
        .start
        .checked_add(u64::try_from(segment.sequence.len())?)
        .ok_or_eyre("Haplotype feature reference overflow")?;
    let mut total_regions = 0;
    let mut scored = 0;
    for domain in mask.unmasked(segment.contig.as_str(), segment.start..reference_end) {
        let reference = Reference {
            start: u32::try_from(domain.start)?,
            bases: segment
                .sequence
                .get(
                    usize::try_from(domain.start - segment.start)?
                        ..usize::try_from(domain.end - segment.start)?,
                )
                .ok_or_eyre("Invalid haplotype mask interval")?
                .iter()
                .copied()
                .map(Base::from)
                .collect(),
        };
        let observed_sites = columns
            .iter()
            .filter(|column| domain.contains(&u64::from(column.pos)))
            .filter_map(|column| {
                let alts: Vec<_> = column
                    .alts
                    .iter()
                    .filter(|alt| alt.metrics.depth >= 2)
                    .map(|alt| alt.base)
                    .collect();
                (!alts.is_empty() && column.reference_base != Base::Unknown)
                    .then(|| Site::new(column.pos, column.reference_base, alts))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut candidates = Vec::new();
        for (index, column) in columns.iter().enumerate().filter(|(_, column)| {
            domain.contains(&u64::from(column.pos)) && column.reference_base != Base::Unknown
        }) {
            // Only observed alleles enter the search. Proposing every base at
            // every candidate spent three of the beam's slots per site on
            // alleles nothing supports and cut regions at ~21 sites.
            let mut alleles = vec![column.reference_base];
            let mut snps = false;
            for alt in &column.alts {
                if alt.base.known_index().is_none()
                    || alleles.contains(&alt.base)
                    || !admits_alt(column, &alt.metrics)
                {
                    continue;
                }
                alleles.push(alt.base);
                snps |= column
                    .alt_metrics(alt.base)
                    .is_some_and(|candidate| crate::call::process::pre_ml_filter(&candidate));
            }
            let mut indels = Vec::new();
            if let Some(data) = &column.indel_data {
                for allele in data.counts.alleles.iter().filter(|allele| allele.total >= 2) {
                    let (deleted, inserted) = match &allele.allele {
                        IndelAllele::Insertion(bases) => (Vec::new(), bases.to_vec()),
                        IndelAllele::Deletion(bases) => (bases.to_vec(), Vec::new()),
                    };
                    let edit = Indel {
                        start: column.pos.checked_add(1).ok_or_eyre("Indel overflow")?,
                        deleted,
                        inserted,
                    };
                    if edit.len() > usize::from(params.max_indel_length)
                        || reference.base_at(edit.start).is_none()
                        || edit
                            .start
                            .checked_add(u32::try_from(edit.deleted.len())?)
                            .is_none_or(|end| reference.base_at(end).is_none())
                        || observed_sites.iter().any(|site| edit.deletes(site.pos))
                    {
                        continue;
                    }
                    indels.push(edit);
                }
            }
            if indels.len() >= 8 {
                indels.clear();
            }
            if snps || !indels.is_empty() || alleles.len() > 1 {
                candidates.push((
                    index,
                    super::regional::Candidate {
                        site: Site { pos: column.pos, reference: column.reference_base, alleles },
                        indels,
                    },
                    snps,
                ));
            }
        }
        // A site that is observed but not a candidate is context: it only
        // matters inside the window of a candidate's reads, and as a site of
        // its own it would only shorten regions.
        let scored_positions: Vec<u32> = candidates
            .iter()
            .filter(|(_, candidate, snps)| *snps || !candidate.indels.is_empty())
            .map(|(_, candidate, _)| candidate.site.pos)
            .collect();
        candidates.retain(|(_, candidate, snps)| {
            *snps || !candidate.indels.is_empty() || {
                let next = scored_positions.partition_point(|&pos| pos < candidate.site.pos);
                [next.checked_sub(1), Some(next)]
                    .into_iter()
                    .flatten()
                    .filter_map(|index| scored_positions.get(index))
                    .any(|pos| pos.abs_diff(candidate.site.pos) <= search.window_radius)
            }
        });
        let mut by_position = BTreeMap::<u32, Vec<usize>>::new();
        for (index, fragment) in fragments.iter().enumerate() {
            for read in &fragment.reads {
                if let Some((first, last)) = read
                    .sequence
                    .as_ref()
                    .and_then(|sequence| sequence.bases.first().zip(sequence.bases.last()))
                {
                    for bucket in first.position.pos / 1024..=last.position.pos / 1024 {
                        by_position.entry(bucket).or_default().push(index);
                    }
                }
            }
        }
        for indices in by_position.values_mut() {
            indices.sort_unstable();
            indices.dedup();
        }
        let mut regions = Vec::<std::ops::Range<usize>>::new();
        let mut begin = 0;
        for index in 1..candidates.len() {
            let previous =
                &candidates.get(index - 1).ok_or_eyre("Missing previous candidate")?.1.site;
            let site = &candidates.get(index).ok_or_eyre("Missing candidate")?.1.site;
            let first = &candidates.get(begin).ok_or_eyre("Missing region start")?.1.site;
            let linked = by_position
                .get(&(previous.pos / 1024))
                .into_iter()
                .flatten()
                .filter_map(|&i| fragments.get(i))
                .filter(|fragment| {
                    fragment
                        .reads
                        .iter()
                        .any(|read| read.covers(previous.pos, params.likelihood_engine))
                        && fragment
                            .reads
                            .iter()
                            .any(|read| read.covers(site.pos, params.likelihood_engine))
                })
                .take(usize::from(search.min_linking_fragments))
                .count()
                >= usize::from(search.min_linking_fragments);
            if index - begin >= usize::from(search.max_sites)
                || site.pos.saturating_sub(first.pos) > search.expansion_radius * 2
                || !linked
            {
                regions.push(begin..index);
                begin = index;
            }
        }
        if begin < candidates.len() {
            regions.push(begin..candidates.len());
        }
        let results = regions
            .par_iter()
            .map(|range| -> Result<_> {
                let offered =
                    candidates.get(range.clone()).ok_or_eyre("Invalid regional candidate range")?;
                let mut ids = BTreeSet::new();
                for (_, candidate, _) in offered {
                    let end = candidate
                        .indels
                        .iter()
                        .map(|edit| edit.equivalent_end(&reference))
                        .collect::<Result<Vec<_>>>()?
                        .into_iter()
                        .max()
                        .unwrap_or(candidate.site.pos);
                    for bucket in candidate.site.pos / 1024..=end / 1024 {
                        if let Some(indices) = by_position.get(&bucket) {
                            ids.extend(indices);
                        }
                    }
                }
                // Fragments are in position order, so the buckets' members
                // are one run of them; a fragment in between that reaches no
                // locus is left out of every matrix anyway.
                let local = match ids.first().zip(ids.last()) {
                    Some((&&first, &&last)) => {
                        fragments.get(first..=last).ok_or_eyre("Invalid regional fragment range")?
                    }
                    None => return Ok((0usize, Vec::new())),
                };
                let mut prepared = Vec::new();
                for (index, candidate, snps) in offered {
                    let mut candidate = candidate.clone();
                    let mut safe = Vec::new();
                    for edit in candidate.indels {
                        if !complex_repeat_context(&reference, &edit, &observed_sites, local)? {
                            safe.push(edit);
                        }
                    }
                    candidate.indels = safe;
                    if candidate.site.alleles.len() > 1 || !candidate.indels.is_empty() {
                        prepared.push((*index, candidate, *snps));
                    }
                }
                let mut scorer = super::Scorer::new(&params)?;
                let mut offset = 0;
                let mut result = Vec::new();
                let mut blocks = 0;
                let offered: Vec<_> =
                    prepared.iter().map(|(_, candidate, _)| candidate.clone()).collect();
                while offset < offered.len() {
                    let region = super::regional::infer(
                        &reference,
                        offered.get(offset..).ok_or_eyre("Invalid remaining candidates")?,
                        local,
                        &params,
                        &search,
                        &mut scorer,
                    )?;
                    blocks += 1;
                    for (index, candidate, snps) in prepared.iter().skip(offset).take(region.used) {
                        if !snps && candidate.indels.is_empty() {
                            continue;
                        }
                        let mut evidence = ColumnEvidence::default();
                        if region.realigned > 0 && *snps {
                            evidence.snp_genotype = Some(super::genotype::DiploidPosterior::pairs(
                                &region.haplotypes,
                                &region.pairs,
                                candidate.site.pos,
                                region.pruning_mass,
                            )?);
                            for &base in &candidate.site.alleles {
                                if base == candidate.site.reference {
                                    continue;
                                }
                                let slot = evidence
                                    .snps
                                    .get_mut(
                                        base.known_index().ok_or_eyre("Unknown regional allele")?,
                                    )
                                    .ok_or_eyre("Unknown feature allele")?;
                                *slot = regional_feature(
                                    &region,
                                    local,
                                    candidate.site.pos,
                                    |haplotype| haplotype.allele_at(candidate.site.pos) == Some(base),
                                )?;
                                slot.haplotype_genotype_entropy = evidence
                                    .snp_genotype
                                    .as_ref()
                                    .ok_or_eyre("Missing regional SNP posterior")?
                                    .entropy()
                                    .f();
                            }
                        }
                        result.push((*index, evidence));
                    }
                    offset += region.used;
                }
                Ok((blocks, result))
            })
            .collect::<Result<Vec<_>>>()?;
        for (blocks, results) in results {
            total_regions += blocks;
            for (index, evidence) in results {
                columns.get_mut(index).ok_or_eyre("Missing regional output column")?.phase =
                    Some(Box::new(ColumnPhase::Likelihood(evidence)));
                scored += 1;
            }
        }
    }
    tracing::debug!(
        scored,
        regions = total_regions,
        fragments = fragments.len(),
        "Calculated regional haplotype ML features"
    );
    Ok(())
}

/// GATK's informative-read threshold, log10 0.2, in nats.
const INFORMATIVE_NATS: f64 = 0.2 * std::f64::consts::LN_10;

/// Allele fraction the heterozygote test uses, under one half so that an allele
/// the informative strand under-observes still passes.
const HET_GATE_FRACTION: f64 = 0.4;
/// A column whose informative strand leaves a heterozygote this improbable is
/// not a regional site.
const HET_GATE_P: f64 = 0.01;
/// Fewer informative reads than this exclude nothing.
const HET_GATE_MIN_READS: u32 = 3;

/// Whether an observed allele earns a regional site. Admitting every column
/// with a mismatching read made the regional model the whole cost of a call:
/// on chr12 it scored 39,000 sites per Mb, 83 % of them with one alternate
/// read and nearly every methylated CpG among them, for 3 % of the records
/// that pass. Two reads are required, and the strand TAPS conversion leaves
/// informative for the allele must not rule a heterozygote out.
fn admits_alt(column: &PileupMetrics, alt: &AlleleMetrics) -> bool {
    if alt.depth < 2 {
        return false;
    }
    let (alt_reads, reads) = match informative_strand(column.reference_base, alt.base) {
        Some(strand) => {
            let on = |metrics: &AlleleMetrics| match strand {
                Strand::OT => metrics.strand_count.ot,
                Strand::OB | Strand::Unknown => metrics.strand_count.ob,
            };
            (on(alt), column.ref_alts_metrics().map(on).sum())
        }
        None => (alt.depth, column.ref_alts_metrics().map(|metrics| metrics.depth).sum()),
    };
    reads < HET_GATE_MIN_READS || binomial_cdf(alt_reads, reads, HET_GATE_FRACTION) >= HET_GATE_P
}

/// The strand that observes `alt` over `reference` without TAPS conversion in
/// the way: a `T` over a `C` may be a converted methylated `C` on OT and an
/// `A` over a `G` on OB, and a true `C` or `G` allele is itself converted
/// there. `None` when both strands count.
fn informative_strand(reference: Base, alt: Base) -> Option<Strand> {
    match (reference, alt) {
        (Base::C, Base::T) | (_, Base::C) => Some(Strand::OB),
        (Base::G, Base::A) | (_, Base::G) => Some(Strand::OT),
        _ => None,
    }
}

/// P(X <= k) for X ~ Binomial(n, p), summed upwards from the zero term.
fn binomial_cdf(k: u32, n: u32, p: f64) -> f64 {
    let mut term = (1.0 - p).powi(i32::try_from(n).unwrap_or(i32::MAX));
    let mut total = term;
    for i in 1..=k.min(n) {
        term *= f64::from(n - i + 1) / f64::from(i) * p / (1.0 - p);
        total += term;
    }
    total.min(1.0)
}

fn regional_feature(
    region: &super::regional::Inferred,
    fragments: &[Fragment],
    pos: u32,
    carries: impl Fn(&super::Haplotype) -> bool,
) -> Result<HaplotypeEvidence> {
    let mut probabilities = [0.0; 3];
    let (mut present, mut absent, mut present_prior, mut absent_prior) =
        (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for pair in &region.pairs {
        let a = region.haplotypes.get(pair.first).ok_or_eyre("Missing feature chromosome")?;
        let b = region.haplotypes.get(pair.second).ok_or_eyre("Missing feature chromosome")?;
        let copies = usize::from(carries(a)) + usize::from(carries(b));
        *probabilities.get_mut(copies).ok_or_eyre("Invalid allele copy count")? += *pair.posterior;
        if copies == 0 {
            absent = log_add(absent, pair.log_likelihood + pair.log_prior);
            absent_prior = log_add(absent_prior, pair.log_prior);
        } else {
            present = log_add(present, pair.log_likelihood + pair.log_prior);
            present_prior = log_add(present_prior, pair.log_prior);
        }
    }
    let [_, het, hom] = probabilities;
    let bayes = if present_prior.is_finite() && absent_prior.is_finite() {
        (present - present_prior) - (absent - absent_prior)
    } else {
        0.0
    };
    // The forest is trained at one depth and applied at others, so nothing
    // here may scale with coverage: the Bayes factor, a sum over molecules,
    // becomes a mean over the fragments covering the site, and the counts
    // below become fractions of the molecules that could have produced them.
    let covers = |fragment: &Fragment, at: u32| {
        fragment.reads.iter().any(|read| read.covers(at, super::LikelihoodEngine::PairHmm))
    };
    let (mut covering, mut linking) = (0usize, 0usize);
    for fragment in fragments.iter().filter(|fragment| covers(fragment, pos)) {
        covering += 1;
        if region.sites.iter().any(|site| site.pos != pos && covers(fragment, site.pos)) {
            linking += 1;
        }
    }
    let bayes_per_molecule = (bayes.f() / covering.max(1).f()).clamp(-50.0, 50.0);
    // The prior-free view the forest needs beside the posterior: which dosage
    // the molecules favour given the allele, and how many of them its best
    // pair assigns to the carrying chromosome (GATK's allele depth), counting
    // a fragment only where the two chromosomes' likelihoods differ by more
    // than a fifth of a decade.
    let conditional_homozygous = if het + hom > 0.0 { hom / (het + hom) } else { 0.0 };
    let (mut molecules, mut informative) = (0usize, 0usize);
    if let Some(best) = region
        .pairs
        .iter()
        .filter(|pair| {
            region.haplotypes.get(pair.first).is_some_and(&carries)
                || region.haplotypes.get(pair.second).is_some_and(&carries)
        })
        .max_by(|a, b| a.posterior.total_cmp(&b.posterior))
    {
        let first = region.likelihoods.get(best.first).ok_or_eyre("Missing chromosome row")?;
        let second = region.likelihoods.get(best.second).ok_or_eyre("Missing chromosome row")?;
        let carriers = [best.first, best.second]
            .map(|index| region.haplotypes.get(index).is_some_and(&carries));
        for (&a, &b) in first.iter().zip(second) {
            if (a - b).abs() <= INFORMATIVE_NATS {
                continue;
            }
            informative += 1;
            molecules += usize::from(if a > b { carriers[0] } else { carriers[1] });
        }
    }
    Ok(HaplotypeEvidence {
        haplotype_available: 1.0,
        haplotype_alt_probability: (het + hom).f(),
        haplotype_heterozygous_probability: het.f(),
        haplotype_homozygous_probability: hom.f(),
        haplotype_expected_fraction: (het / 2.0 + hom).f(),
        haplotype_log_bayes_factor_per_molecule: bayes_per_molecule,
        haplotype_genotype_entropy: (-probabilities
            .iter()
            .filter(|&&p| p > 0.0)
            .map(|p| p * p.ln())
            .sum::<f64>())
        .f(),
        haplotype_linking_fraction: if covering > 0 { linking.f() / covering.f() } else { 0.0 },
        haplotype_informative_fraction: informative.f() / region.fragments.max(1).f(),
        haplotype_realigned_fraction: region.realigned.f() / region.fragments.max(1).f(),
        haplotype_pruned_mass: region.pruning_mass.f(),
        haplotype_conditional_homozygous: conditional_homozygous.f(),
        haplotype_molecular_alt_fraction: if informative > 0 {
            molecules.f() / informative.f()
        } else {
            0.0
        },
    })
}

fn complex_repeat_context(
    reference: &Reference,
    indel: &Indel,
    sites: &[Site],
    fragments: &[Fragment],
) -> Result<bool> {
    // A normalized gap stands for every placement up to its reference-
    // equivalent end; a SNP inside that span is on an unknown chromosome.
    let end = indel.equivalent_end(reference)?;
    let can_be = |pos, base| {
        reference.base_at(pos) == Some(base)
            || sites.iter().any(|site| site.pos == pos && site.alleles.contains(&base))
    };
    for site in sites.iter().filter(|site| (indel.start..end).contains(&site.pos)) {
        let mut support = 0;
        for fragment in fragments {
            let observed = super::likelihood::consensus_at(fragment, site.pos)?;
            let Some(read) = observed.source.and_then(|source| fragment.reads.get(source)) else {
                continue;
            };
            if observed.base == site.reference
                || !observed.quality.is_some_and(|q| q >= 20)
                || !read.mapq.is_some_and(|q| q >= 20)
            {
                continue;
            }
            let converted = (site.alleles.contains(&Base::C)
                && observed.base == Base::T
                && read.strand != Strand::OB
                && site.pos.checked_add(1).is_some_and(|pos| can_be(pos, Base::G)))
                || (site.alleles.contains(&Base::G)
                    && observed.base == Base::A
                    && read.strand != Strand::OT
                    && site.pos.checked_sub(1).is_some_and(|pos| can_be(pos, Base::C)));
            if converted {
                continue;
            }
            support += 1;
            // Reference-equivalent gap shifts are not necessarily equivalent
            // on a SNP-bearing chromosome. This caller has no complex-allele
            // assembler, so it must not invent that chromosome's likelihood.
            if support >= 2 {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::super::reads::fixture_read;
    use super::*;
    use seqair::bam::cigar::CigarOpType::Match;

    fn reference(sequence: &[u8]) -> Reference {
        Reference { start: 1000, bases: sequence.iter().copied().map(Base::from).collect() }
    }

    /// Converted reads look like a `C>T` allele, so the regional model used to
    /// run at every methylated CpG; the strand TAPS cannot convert decides.
    #[test]
    fn a_methylated_cpg_is_not_a_regional_site_but_a_c_to_t_het_is() -> Result<()> {
        use crate::call::tests::utils::{ReadLine, create_pileups};
        let reference = reference(b"GACGAATG");
        let read = |bases: &[u8], strand: Strand, copies: usize| {
            let bases: Vec<Base> = bases.iter().copied().map(Base::from).collect();
            (0..copies).map(move |_| ReadLine { bases: bases.clone(), strand })
        };
        let available = |reads: Vec<ReadLine>| -> Result<f32> {
            let (segment, pileups) =
                create_pileups(b"GACGAATG".iter().copied().map(Base::from).collect(), reads);
            let mut columns =
                pileups.into_iter().map(PileupMetrics::new).collect::<Result<Vec<_>>>()?;
            let alt = fixture_read(&reference, b"GATGAATG", &[(Match, 8)])?;
            let wt = fixture_read(&reference, b"GACGAATG", &[(Match, 8)])?;
            let fragments = [alt.clone(), alt, wt.clone(), wt]
                .into_iter()
                .map(|read| Fragment { reads: vec![read] })
                .collect::<Vec<_>>();
            add_likelihoods(&mut columns, &segment, &fragments, &PhaseMask::default())?;
            let focal = columns.iter().find(|c| c.pos == 1002).ok_or_eyre("focal")?;
            let candidate = focal.alt_metrics(Base::T).ok_or_eyre("candidate")?;
            Ok(ColumnEvidence::snp(&candidate).haplotype_available)
        };
        // Every OT read converted, every OB read a C: methylation, not a variant.
        let methylated =
            read(b"GATGAATG", Strand::OT, 10).chain(read(b"GACGAATG", Strand::OB, 10)).collect();
        assert_eq!(available(methylated)?, 0.0);
        // Half the reads on both strands carry T: a heterozygote.
        let het = read(b"GATGAATG", Strand::OT, 5)
            .chain(read(b"GACGAATG", Strand::OT, 5))
            .chain(read(b"GATGAATG", Strand::OB, 5))
            .chain(read(b"GACGAATG", Strand::OB, 5))
            .collect();
        assert_eq!(available(het)?, 1.0);
        // Too few OB reads to rule a heterozygote out: still scored.
        let thin = read(b"GATGAATG", Strand::OT, 10).chain(read(b"GACGAATG", Strand::OB, 2)).collect();
        assert_eq!(available(thin)?, 1.0);
        // One mismatching read is no site at all.
        let single = read(b"GATGAATG", Strand::OB, 1).chain(read(b"GACGAATG", Strand::OB, 9)).collect();
        assert_eq!(available(single)?, 0.0);
        Ok(())
    }

    #[test]
    fn pre_filter_features_need_no_called_genotype_and_masks_only_remove_evidence() -> Result<()> {
        let reference = reference(b"GAAACATG");
        let (segment, pileups) = crate::pileups! {
            [G A A A C A T G] Ref,
            [G A T A C A T G] OT,
            [G A T A C A T G] OB,
            [G A A A C A T G] OT,
            [G A A A C A T G] OB,
        };
        let mut columns =
            pileups.into_iter().map(PileupMetrics::new).collect::<Result<Vec<_>>>()?;
        let alt = fixture_read(&reference, b"GATACATG", &[(Match, 8)])?;
        let wt = fixture_read(&reference, b"GAAACATG", &[(Match, 8)])?;
        let fragments = [alt.clone(), alt, wt.clone(), wt]
            .into_iter()
            .map(|read| Fragment { reads: vec![read] })
            .collect::<Vec<_>>();
        add_likelihoods(&mut columns, &segment, &fragments, &PhaseMask::default())?;
        let focal = columns.iter().find(|c| c.pos == 1002).ok_or_eyre("focal")?;
        let candidate = focal.alt_metrics(Base::T).ok_or_eyre("candidate")?;
        let row = ColumnEvidence::snp(&candidate);
        assert_eq!(row.haplotype_available, 1.0);
        assert_eq!(row.haplotype_realigned_fraction, 1.0);
        assert_eq!(row.haplotype_informative_fraction, 1.0, "{row:?}");
        assert_eq!(row.haplotype_molecular_alt_fraction, 0.5, "{row:?}");
        assert_eq!(row.haplotype_linking_fraction, 0.0, "{row:?}");
        assert!(row.haplotype_log_bayes_factor_per_molecule > 0.0, "{row:?}");
        assert_eq!(candidate.metrics.pos_metrics.depth, 4);
        for column in &mut columns {
            column.phase = None;
        }
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("mask.bed");
        std::fs::write(&path, "chr_test\t1002\t1003\n")?;
        let mask = PhaseMask::load(&[path])?;
        add_likelihoods(&mut columns, &segment, &fragments, &mask)?;
        let focal = columns.iter().find(|c| c.pos == 1002).ok_or_eyre("focal")?;
        let candidate = focal.alt_metrics(Base::T).ok_or_eyre("candidate")?;
        assert_eq!(ColumnEvidence::snp(&candidate).haplotype_available, 0.0);
        assert_eq!(candidate.alt.depth, 2);
        Ok(())
    }

    #[test]
    fn conversion_cannot_mark_a_repeat_as_genetically_broken() -> Result<()> {
        let reference = reference(b"ACGCGT");
        let indel = Indel { start: 1001, deleted: vec![Base::C, Base::G], inserted: Vec::new() };
        let read = fixture_read(&reference, b"ATGCGT", &[(Match, 6)])?;
        let sites = [Site::new(1001, Base::C, [Base::T])?];
        let mut fragments = vec![Fragment { reads: vec![read.clone()] }; 2];
        assert!(!complex_repeat_context(&reference, &indel, &sites, &fragments)?);
        for fragment in &mut fragments {
            for read in &mut fragment.reads {
                read.strand = Strand::OB;
            }
        }
        assert!(complex_repeat_context(&reference, &indel, &sites, &fragments)?);
        let mut only_one = fragments.first().ok_or_eyre("fragment")?.clone();
        only_one.reads.push(read);
        assert!(!complex_repeat_context(&reference, &indel, &sites, &[only_one])?);
        Ok(())
    }
}
