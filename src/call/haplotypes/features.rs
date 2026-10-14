use super::{
    ExpansionParams, Fragment, Indel, LikelihoodParams, Reference, Site, likelihood::log_add,
    reads::NormalizedObservation,
};
use crate::utils::IntoF32 as _;
use crate::{
    call::{
        phasing::{ColumnPhase, PhaseMask},
        pileup::{
            indels::{
                IndelAllele, IndelAlleleCounts, IndelCounts, IndelData, IndelPositionMetrics,
            },
            ref_features::{indel_ref_window_at, indel_tract_runs_at},
        },
    },
    metrics::{MetricsForAlt, MetricsForIndel, PileupMetrics, ml::features::HaplotypeEvidence},
    sequence::Segment,
};
use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, RmsAccumulator, Strand};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Default)]
pub struct ColumnEvidence {
    snps: [HaplotypeEvidence; 4],
    indels: Vec<(IndelAllele, HaplotypeEvidence)>,
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
    pub(crate) fn indel(current: &MetricsForIndel) -> HaplotypeEvidence {
        current
            .metrics
            .phase
            .as_deref()
            .and_then(|phase| match phase {
                ColumnPhase::Likelihood(evidence) => evidence
                    .indels
                    .iter()
                    .find(|(allele, _)| allele == &current.indel.allele)
                    .map(|(_, evidence)| *evidence),
                _ => None,
            })
            .unwrap_or_default()
    }
}

pub(crate) fn parameters() -> Result<LikelihoodParams> {
    use seqair_types::Probability;
    Ok(LikelihoodParams {
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

pub(crate) fn normalize_indels(
    columns: &mut [PileupMetrics],
    segment: &Segment,
    fragments: &[Fragment],
) -> Result<()> {
    if fragments.is_empty() {
        return Ok(());
    }
    let mut candidates = BTreeMap::<u32, Vec<(usize, NormalizedObservation)>>::new();
    for (index, fragment) in fragments.iter().enumerate() {
        let mut molecule = BTreeMap::<u32, NormalizedObservation>::new();
        for observation in fragment
            .reads
            .iter()
            .filter_map(|read| read.annotations.as_ref())
            .flat_map(|a| &a.indels)
        {
            let Some(anchor) = observation.indel.start.checked_sub(1) else { continue };
            if u64::from(anchor) < segment.start {
                continue;
            }
            let mut observation = observation.clone();
            if let Some(previous) = molecule.get_mut(&anchor) {
                if previous.indel == observation.indel {
                    if (observation.observation.mapq, observation.observation.base_qual)
                        > (previous.observation.mapq, previous.observation.base_qual)
                    {
                        std::mem::swap(previous, &mut observation);
                    }
                    previous.observation.noisy |= observation.observation.noisy;
                } else {
                    // Different alleles from the two mates are one uncertain molecule.
                    if observation.observation.base_qual > previous.observation.base_qual {
                        std::mem::swap(previous, &mut observation);
                    }
                    previous.observation.base_qual = previous
                        .observation
                        .base_qual
                        .saturating_sub(observation.observation.base_qual);
                    previous.observation.noisy = true;
                }
            } else {
                molecule.insert(anchor, observation);
            }
        }
        for (anchor, observation) in molecule {
            candidates.entry(anchor).or_default().push((index, observation));
        }
    }
    let positions: BTreeSet<_> = columns.iter().map(|p| p.pos).collect();
    let reference = Reference {
        start: u32::try_from(segment.start)?,
        bases: segment.sequence.iter().copied().map(Base::from).collect(),
    };
    // Keep evidence that cannot be moved into a represented column, rather than
    // silently deleting long or boundary-truncated events.
    for column in columns.iter_mut() {
        if let Some(data) = column.indel_data.as_mut() {
            let mut retained = Vec::new();
            for observation in &data.observations {
                let moved = match canonical_allele(&reference, column.pos, &observation.allele) {
                    Ok((anchor, allele)) => {
                        positions.contains(&anchor)
                            && candidates.get(&anchor).is_some_and(|offered| {
                                offered.iter().any(|(_, obs)| obs.observation.allele == allele)
                            })
                    }
                    Err(error) => {
                        tracing::debug!(
                            ?error,
                            pos = column.pos,
                            "Keeping an indel whose canonical position is outside the captured reference"
                        );
                        false
                    }
                };
                if !moved {
                    retained.push(observation.clone());
                }
            }
            data.observations = retained.into_iter().collect();
            if data.observations.is_empty() {
                column.indel_data = None;
            } else {
                data.counts = crate::metrics::aggregate_indels(
                    &data.observations,
                    column.pos_metrics.depth as usize,
                    data.counts.noisy_ref_count,
                    column.pos,
                );
            }
        }
    }
    for column in columns.iter_mut() {
        let Some(observations) = candidates.remove(&column.pos) else { continue };
        let mut alleles = Vec::<IndelAlleleCounts>::new();
        let mut source = Vec::new();
        let mut alternate_fragments = BTreeSet::new();
        let mut aligned_boundaries = 0usize;
        for (fragment, observation) in observations {
            // Collection filtered the original allele observation. Normalization
            // can move its boundary to a lower-quality or unsequenced base;
            // retain the molecule and let those boundary features express uncertainty.
            alternate_fragments.insert(fragment);
            aligned_boundaries += usize::from(observation.boundary_aligned);
            let allele = &observation.observation.allele;
            let index = if let Some(index) = alleles.iter().position(|a| &a.allele == allele) {
                index
            } else {
                alleles.push(IndelAlleleCounts {
                    allele: allele.clone(),
                    ot: 0,
                    ob: 0,
                    unknown_strand: 0,
                    noisy: 0,
                });
                alleles.len() - 1
            };
            let entry = alleles.get_mut(index).ok_or_eyre("Missing canonical indel allele")?;
            match observation.observation.strand {
                Strand::OT => entry.ot += 1,
                Strand::OB => entry.ob += 1,
                Strand::Unknown => entry.unknown_strand += 1,
            }
            entry.noisy += u32::from(observation.observation.noisy);
            source.push(observation.observation);
        }
        if source.is_empty() {
            continue;
        }
        let mut depth = 0u32;
        let mut noisy_ref = 0;
        let mut soft_clipped = 0;
        let mut mapq = RmsAccumulator::new();
        let mut baseq = RmsAccumulator::new();
        for (index, fragment) in fragments.iter().enumerate() {
            let base = super::likelihood::consensus_at(fragment, column.pos)?;
            if base.source.is_none() && !alternate_fragments.contains(&index) {
                continue;
            }
            depth += 1;
            if let Some(read) = base.source.and_then(|index| fragment.reads.get(index)) {
                if let Some(q) = read.mapq {
                    mapq.add(f64::from(q));
                }
            } else if let Some(q) = fragment.reads.iter().filter_map(|read| read.mapq).max() {
                mapq.add(f64::from(q));
            }
            if let Some(q) = base.quality {
                baseq.add(f64::from(q));
            }
            let annotations = fragment.reads.iter().filter_map(|read| read.annotations.as_ref());
            let noisy = annotations.clone().any(|a| a.noisy);
            soft_clipped += u32::from(annotations.clone().any(|a| a.soft_clipped));
            noisy_ref += u32::from(noisy && !alternate_fragments.contains(&index));
        }
        let idx = usize::try_from(u64::from(column.pos) - segment.start)?;
        let (ref_window, ref_anchor) = indel_ref_window_at(idx, segment);
        let runs = indel_tract_runs_at(u64::from(column.pos), segment);
        let boundary_aligned_fraction = aligned_boundaries.f() / source.len().max(1).f();
        if let Some(previous) = column.indel_data.as_ref() {
            source.extend(previous.observations.iter().cloned());
            for allele in &previous.counts.alleles {
                if !alleles.iter().any(|a| a.allele == allele.allele) {
                    alleles.push(allele.clone());
                }
            }
        }
        column.indel_data = Some(Box::new(IndelData {
            observations: source.into_iter().collect(),
            ref_window,
            ref_anchor,
            homopolymer_run: runs.homopolymer,
            dinucleotide_run: runs.dinucleotide,
            soft_clip_count: soft_clipped,
            counts: IndelCounts {
                alleles: alleles.into_iter().collect(),
                ref_count: depth.saturating_sub(u32::try_from(alternate_fragments.len())?),
                noisy_ref_count: noisy_ref,
            },
            calls: Vec::new(),
            canonical_position: Some(IndelPositionMetrics {
                depth,
                mapq: mapq.finish().f(),
                baseq: baseq.finish().f(),
                boundary_aligned_fraction,
            }),
        }));
    }
    Ok(())
}

pub(crate) fn canonical_allele(
    reference: &Reference,
    anchor: u32,
    allele: &IndelAllele,
) -> Result<(u32, IndelAllele)> {
    let start = anchor.checked_add(1).ok_or_eyre("Indel coordinate overflow")?;
    let (deleted, inserted) = match allele {
        IndelAllele::Insertion(bases) => (Vec::new(), bases.to_vec()),
        IndelAllele::Deletion(bases) => (bases.to_vec(), Vec::new()),
    };
    let indel = Indel { start, deleted, inserted }.normalized(reference)?;
    let anchor = indel.start.checked_sub(1).ok_or_eyre("Indel has no left anchor")?;
    ensure!(reference.base_at(anchor).is_some(), "Indel normalization reached reference boundary");
    let allele = if indel.inserted.is_empty() {
        IndelAllele::Deletion(indel.deleted.into_iter().collect())
    } else {
        IndelAllele::Insertion(indel.inserted.into_iter().collect())
    };
    Ok((anchor, allele))
}

pub(crate) fn add_likelihoods(
    columns: &mut [PileupMetrics],
    segment: &Segment,
    fragments: &[Fragment],
    mask: &PhaseMask,
) -> Result<()> {
    let reference_end = segment
        .start
        .checked_add(u64::try_from(segment.sequence.len())?)
        .ok_or_eyre("Haplotype feature reference overflow")?;
    let params = parameters()?;
    let search = ExpansionParams { read_local: true, ..ExpansionParams::default() };
    let mut scored = 0usize;
    for region in mask.unmasked(segment.contig.as_str(), segment.start..reference_end) {
        let reference = Reference {
            start: u32::try_from(region.start)?,
            bases: segment
                .sequence
                .get(
                    usize::try_from(region.start - segment.start)?
                        ..usize::try_from(region.end - segment.start)?,
                )
                .ok_or_eyre("Invalid haplotype mask interval")?
                .iter()
                .copied()
                .map(Base::from)
                .collect(),
        };
        let sites: Vec<_> = columns
            .iter()
            .filter(|column| region.contains(&u64::from(column.pos)))
            .filter_map(|column| {
                let alternatives = column
                    .alts
                    .iter()
                    .filter(|alt| alt.metrics.depth >= 2)
                    .map(|alt| alt.base)
                    .collect::<Vec<_>>();
                (!alternatives.is_empty() && column.reference_base != Base::Unknown)
                    .then(|| Site::new(column.pos, column.reference_base, alternatives))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut by_position = BTreeMap::<u32, Vec<usize>>::new();
        for (index, fragment) in fragments.iter().enumerate() {
            for read in &fragment.reads {
                if let Some(sequence) = &read.sequence
                    && let Some((first, last)) = sequence.bases.first().zip(sequence.bases.last())
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
        for column in columns.iter_mut().filter(|column| {
            region.contains(&u64::from(column.pos)) && column.reference_base != Base::Unknown
        }) {
            if column.alts.is_empty() && column.indel_data.is_none() {
                continue;
            }
            let mut evidence = ColumnEvidence::default();
            let mut indices = BTreeSet::new();
            let start = column.pos.saturating_sub(search.expansion_radius);
            let end = column.pos.saturating_add(search.expansion_radius);
            for bucket in start / 1024..=end / 1024 {
                if let Some(offered) = by_position.get(&bucket) {
                    indices.extend(offered.iter().copied());
                }
            }
            let local = match indices.first().zip(indices.last()) {
                Some((&first, &last)) => {
                    fragments.get(first..=last).ok_or_eyre("Invalid molecular feature interval")?
                }
                None => continue,
            };
            if local.is_empty() {
                continue;
            }
            let neighbors: Vec<_> = sites
                .iter()
                .filter(|site| site.pos.abs_diff(column.pos) <= search.expansion_radius)
                .cloned()
                .collect();
            let indels: Vec<_> = local
                .iter()
                .flat_map(|fragment| &fragment.reads)
                .filter_map(|read| read.sequence.as_ref())
                .flat_map(|sequence| &sequence.indels)
                .filter(|indel| {
                    indel.start > reference.start
                        && reference.base_at(indel.start).is_some()
                        && u32::try_from(indel.deleted.len())
                            .ok()
                            .and_then(|len| indel.start.checked_add(len))
                            .is_some_and(|end| reference.base_at(end).is_some())
                })
                .cloned()
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            let mut unsafe_indels = BTreeSet::new();
            for indel in &indels {
                if complex_repeat_context(&reference, indel, &neighbors, local)? {
                    unsafe_indels.insert(indel.clone());
                }
            }
            let score_snps = column.alts.iter().any(|alt| {
                column
                    .alt_metrics(alt.base)
                    .is_some_and(|candidate| crate::call::process::pre_ml_filter(&candidate))
            });
            if score_snps {
                if unsafe_indels.iter().any(|indel| indel.start.abs_diff(column.pos) <= 150) {
                    continue;
                }
                let result = super::infer_snp_expanded(
                    &reference, column.pos, &neighbors, &indels, local, &params, &search,
                )?;
                if result
                    .haplotypes
                    .iter()
                    .filter_map(|haplotype| haplotype.indel.as_ref())
                    .any(|indel| unsafe_indels.contains(indel))
                {
                    continue;
                }
                for allele in result.alleles {
                    let Some(index) = allele.base.known_index() else { continue };
                    let row = evidence
                        .snps
                        .get_mut(index)
                        .ok_or_eyre("Unknown haplotype feature allele")?;
                    *row = HaplotypeEvidence {
                        haplotype_available: 1.0,
                        haplotype_alt_probability: allele.posterior.f(),
                        haplotype_heterozygous_probability: allele.heterozygous.f(),
                        haplotype_homozygous_probability: allele.homozygous.f(),
                        haplotype_expected_fraction: allele.expected_fraction.f(),
                        haplotype_log_bayes_factor: allele
                            .log_bayes_factor
                            .clamp(-100.0, 100.0)
                            .f(),
                        haplotype_genotype_entropy: result.genotype_entropy.f(),
                        haplotype_context_sites: result.sites.len().f(),
                        haplotype_linking_fragments: result.retained_linking_fragments.f(),
                        haplotype_fragments: result.fragments.f(),
                        haplotype_realigned_fraction: result.realigned_fragments.f()
                            / result.fragments.max(1).f(),
                        haplotype_pruned_mass: result.maximum_pruned_mass.f(),
                        haplotype_boundary_resolved: 1.0,
                    };
                }
            }
            if let Some(data) = column.indel_data.as_ref() {
                for allele in &data.counts.alleles {
                    if allele.total() < 2 || allele.allele.len() > 100 {
                        continue;
                    }
                    let (deleted, inserted) = match &allele.allele {
                        IndelAllele::Insertion(bases) => (Vec::new(), bases.to_vec()),
                        IndelAllele::Deletion(bases) => (bases.to_vec(), Vec::new()),
                    };
                    let indel = Indel {
                        start: column.pos.checked_add(1).ok_or_eyre("Indel position overflow")?,
                        deleted,
                        inserted,
                    };
                    if reference.base_at(indel.start).is_none()
                        || u32::try_from(indel.deleted.len())
                            .ok()
                            .and_then(|len| indel.start.checked_add(len))
                            .is_none_or(|end| reference.base_at(end).is_none())
                    {
                        continue;
                    }
                    if complex_repeat_context(&reference, &indel, &neighbors, local)? {
                        continue;
                    }
                    let mut row = indel_likelihood(
                        &reference, column.pos, &indel, &neighbors, local, &params, &search,
                    )?;
                    row.haplotype_boundary_resolved =
                        data.canonical_position.map_or(0.0, |p| p.boundary_aligned_fraction);
                    evidence.indels.push((allele.allele.clone(), row));
                }
            }
            column.phase = Some(Box::new(ColumnPhase::Likelihood(evidence)));
            scored += 1;
        }
    }
    tracing::debug!(
        scored,
        fragments = fragments.len(),
        "Calculated pre-filter haplotype ML features"
    );
    Ok(())
}

fn complex_repeat_context(
    reference: &Reference,
    indel: &Indel,
    sites: &[Site],
    fragments: &[Fragment],
) -> Result<bool> {
    let end = fragments
        .iter()
        .flat_map(|fragment| &fragment.reads)
        .filter_map(|read| read.annotations.as_ref())
        .flat_map(|a| &a.indels)
        .filter(|observation| &observation.indel == indel)
        .map(|observation| observation.raw_start)
        .max()
        .unwrap_or(indel.start)
        .checked_add(u32::try_from(indel.deleted.len())?)
        .ok_or_eyre("Ambiguous indel span overflow")?;
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

fn indel_likelihood(
    reference: &Reference,
    anchor: u32,
    indel: &Indel,
    neighbors: &[Site],
    fragments: &[Fragment],
    params: &LikelihoodParams,
    search: &ExpansionParams,
) -> Result<HaplotypeEvidence> {
    let base = reference.base_at(anchor).ok_or_eyre("Indel anchor outside reference")?;
    let focal = Site { pos: anchor, reference: base, alleles: vec![base] };
    let neighbors: Vec<_> = neighbors
        .iter()
        .filter(|site| site.pos != anchor && !indel.deletes(site.pos))
        .cloned()
        .collect();
    let expanded = super::expansion::expand(
        reference,
        focal,
        &neighbors,
        Some(indel),
        fragments,
        params,
        search,
        true,
    )?;
    let matrix = super::likelihood_matrix_geometry(
        reference,
        &expanded.sites,
        &expanded.haplotypes,
        Some(indel),
        fragments,
        params,
        true,
    )?;
    if matrix.realigned_fragments == 0 {
        return Ok(HaplotypeEvidence::default());
    }
    let pairs = super::posterior_pairs(
        &expanded.haplotypes,
        &expanded.sites,
        &matrix.values,
        params,
        true,
    )?;
    let mut probabilities = [0.0; 3];
    let mut present = f64::NEG_INFINITY;
    let mut absent = f64::NEG_INFINITY;
    let mut present_prior = f64::NEG_INFINITY;
    let mut absent_prior = f64::NEG_INFINITY;
    for pair in pairs {
        let a = expanded.haplotypes.get(pair.first).ok_or_eyre("Missing first indel haplotype")?;
        let b =
            expanded.haplotypes.get(pair.second).ok_or_eyre("Missing second indel haplotype")?;
        let copies = usize::from(a.indel.is_some()) + usize::from(b.indel.is_some());
        *probabilities.get_mut(copies).ok_or_eyre("Invalid indel genotype")? += *pair.posterior;
        if copies == 0 {
            absent = log_add(absent, pair.log_likelihood + pair.log_prior);
            absent_prior = log_add(absent_prior, pair.log_prior);
        } else {
            present = log_add(present, pair.log_likelihood + pair.log_prior);
            present_prior = log_add(present_prior, pair.log_prior);
        }
    }
    let heterozygous =
        *probabilities.get(1).ok_or_eyre("Missing heterozygous indel probability")?;
    let homozygous = *probabilities.get(2).ok_or_eyre("Missing homozygous indel probability")?;
    Ok(HaplotypeEvidence {
        haplotype_available: 1.0,
        haplotype_alt_probability: (heterozygous + homozygous).f(),
        haplotype_heterozygous_probability: heterozygous.f(),
        haplotype_homozygous_probability: homozygous.f(),
        haplotype_expected_fraction: (heterozygous / 2.0 + homozygous).f(),
        haplotype_log_bayes_factor: ((present - present_prior) - (absent - absent_prior))
            .clamp(-100.0, 100.0)
            .f(),
        haplotype_genotype_entropy: -probabilities
            .iter()
            .filter(|&&p| p > 0.0)
            .map(|p| p * p.ln())
            .sum::<f64>()
            .f(),
        haplotype_context_sites: expanded.sites.len().f(),
        haplotype_linking_fragments: fragments
            .iter()
            .filter(|fragment| {
                fragment.reads.iter().any(|read| {
                    read.covers(anchor, params.likelihood_engine)
                        || read
                            .sequence
                            .as_ref()
                            .is_some_and(|sequence| sequence.indels.contains(indel))
                }) && fragment.reads.iter().any(|read| {
                    expanded.sites.iter().any(|site| {
                        site.pos != anchor && read.covers(site.pos, params.likelihood_engine)
                    })
                })
            })
            .count()
            .f(),
        haplotype_fragments: matrix.fragments.f(),
        haplotype_realigned_fraction: matrix.realigned_fragments.f() / matrix.fragments.max(1).f(),
        haplotype_pruned_mass: expanded.maximum_pruned_mass.f(),
        haplotype_boundary_resolved: 1.0,
    })
}

#[cfg(test)]
mod tests {
    use super::super::reads::fixture_read;
    use super::*;
    use seqair::bam::cigar::CigarOpType::{Insertion, Match};

    fn reference(sequence: &[u8]) -> Reference {
        Reference { start: 1000, bases: sequence.iter().copied().map(Base::from).collect() }
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
        assert_eq!(row.haplotype_fragments, 4.0);
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
    fn repeat_placements_merge_at_the_canonical_anchor_once_per_molecule() -> Result<()> {
        let reference = reference(b"GCAAAACTG");
        let (segment, pileups) = crate::pileups! {
            [G C A A A A C T G] Ref,
            [G C A A A A C T G] OT,
            [G C A A A A C T G] OB,
        };
        let mut columns =
            pileups.into_iter().map(PileupMetrics::new).collect::<Result<Vec<_>>>()?;
        let first =
            fixture_read(&reference, b"GCAAAAACTG", &[(Match, 5), (Insertion, 1), (Match, 4)])?;
        let second =
            fixture_read(&reference, b"GCAAAAACTG", &[(Match, 3), (Insertion, 1), (Match, 6)])?;
        let fragments =
            [Fragment { reads: vec![first.clone(), first] }, Fragment { reads: vec![second] }];
        normalize_indels(&mut columns, &segment, &fragments)?;
        let canonical = columns.iter().find(|c| c.pos == 1001).ok_or_eyre("canonical column")?;
        let data = canonical.indel_data.as_ref().ok_or_eyre("canonical indel")?;
        assert_eq!(data.observations.len(), 2);
        assert_eq!(data.counts.alleles.len(), 1);
        assert_eq!(data.counts.alleles.first().ok_or_eyre("allele")?.total(), 2);
        assert_eq!(data.canonical_position.ok_or_eyre("position")?.depth, 2);
        assert_eq!(data.homopolymer_run, 4);
        assert!(columns.iter().filter(|c| c.pos != 1001).all(|c| c.indel_data.is_none()));
        Ok(())
    }

    #[test]
    fn indel_probabilities_use_pair_hmm_and_fallback_has_no_indel_evidence() -> Result<()> {
        let reference = reference(b"GCAAAACTG");
        let alt =
            fixture_read(&reference, b"GCAAAAACTG", &[(Match, 5), (Insertion, 1), (Match, 4)])?;
        let indel = alt
            .annotations
            .as_ref()
            .and_then(|a| a.indels.first())
            .ok_or_eyre("indel")?
            .indel
            .clone();
        let mut fragments: Vec<_> =
            std::iter::repeat_n(Fragment { reads: vec![alt] }, 24).collect();
        let wt = fixture_read(&reference, b"GCAAAACTG", &[(Match, 9)])?;
        fragments.extend(std::iter::repeat_n(Fragment { reads: vec![wt] }, 24));
        let search = ExpansionParams { read_local: true, ..ExpansionParams::default() };
        let row =
            indel_likelihood(&reference, 1001, &indel, &[], &fragments, &parameters()?, &search)?;
        assert_eq!(row.haplotype_available, 1.0);
        assert_eq!(row.haplotype_realigned_fraction, 1.0);
        assert!(row.haplotype_heterozygous_probability > 0.95, "{row:?}");
        assert!(row.haplotype_homozygous_probability < 0.01, "{row:?}");
        for fragment in &mut fragments {
            for read in &mut fragment.reads {
                read.sequence.as_mut().ok_or_eyre("sequence")?.reference_skip = true;
            }
        }
        let row =
            indel_likelihood(&reference, 1001, &indel, &[], &fragments, &parameters()?, &search)?;
        assert_eq!(row.haplotype_available, 0.0);
        Ok(())
    }

    #[test]
    fn conversion_cannot_mark_a_repeat_as_genetically_broken() -> Result<()> {
        let reference = reference(b"ACGCGT");
        let indel = Indel { start: 1001, deleted: vec![Base::C, Base::G], inserted: Vec::new() };
        let mut read = fixture_read(&reference, b"ATGCGT", &[(Match, 6)])?;
        read.annotations.as_mut().ok_or_eyre("annotation")?.indels.push(NormalizedObservation {
            raw_start: 1003,
            indel: indel.clone(),
            boundary_aligned: true,
            observation: super::super::reads::fixture_read(
                &Reference {
                    start: 1000,
                    bases: b"GCAAAACTG".iter().copied().map(Base::from).collect(),
                },
                b"GCAAAAACTG",
                &[(Match, 5), (Insertion, 1), (Match, 4)],
            )?
            .annotations
            .as_ref()
            .and_then(|a| a.indels.first())
            .ok_or_eyre("indel")?
            .observation
            .clone(),
        });
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
