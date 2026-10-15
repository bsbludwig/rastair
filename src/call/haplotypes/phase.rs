use super::{
    Fragment, Haplotype, LikelihoodEngine, LikelihoodParams, Reference, Site,
    likelihood::{log_add, log_sum},
    methylation, pair_hmm,
};
use color_eyre::eyre::{OptionExt as _, Result, ensure};
use rayon::prelude::*;
use seqair_types::{Base, Probability};
use std::collections::BTreeMap;

pub(crate) fn parameters() -> Result<LikelihoodParams> {
    Ok(LikelihoodParams {
        likelihood_engine: LikelihoodEngine::PairHmm,
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
        max_indel_length: 20,
    })
}

#[derive(Clone, Debug)]
pub(crate) struct State {
    pub first: Vec<Base>,
}

pub(crate) struct Posterior {
    pub probabilities: Vec<f64>,
}

/// Genotypes are fixed; only complementary chromosome assignments are eligible.
/// The first site's lower genotype allele anchors the chromosome labels.
pub(crate) fn score(
    reference: &Reference,
    sites: &[Site],
    states: &[State],
    fragments: &[&Fragment],
    params: &LikelihoodParams,
    scorer: &mut pair_hmm::Scorer,
) -> Result<Posterior> {
    ensure!(!sites.is_empty() && !states.is_empty(), "Empty phase hypothesis set");
    let mut counts = BTreeMap::<super::Indel, usize>::new();
    for fragment in fragments {
        let unique: std::collections::BTreeSet<_> = fragment
            .reads
            .iter()
            .filter_map(|read| read.sequence.as_ref())
            .flat_map(|sequence| &sequence.indels)
            .collect();
        for indel in unique {
            *counts.entry(indel.clone()).or_default() += 1;
        }
    }
    let indel = counts
        .into_iter()
        .filter(|(indel, count)| {
            *count >= 2
                && indel.len() <= usize::from(params.max_indel_length)
                && reference.base_at(indel.start.saturating_sub(1)).is_some()
                && u32::try_from(indel.deleted.len())
                    .ok()
                    .and_then(|length| indel.start.checked_add(length))
                    .is_some_and(|end| reference.base_at(end).is_some())
                && !sites.iter().any(|site| indel.deletes(site.pos))
                && sites.iter().any(|site| site.pos.abs_diff(indel.start) <= 150)
        })
        .min_by_key(|(indel, _)| {
            sites.iter().map(|site| site.pos.abs_diff(indel.start)).min().unwrap_or(u32::MAX)
        })
        .map(|(indel, _)| indel);
    let configurations = if indel.is_some() {
        vec![
            (false, false, (1.0 - *params.indel_prior).ln()),
            (false, true, (*params.indel_prior / 4.0).ln()),
            (true, false, (*params.indel_prior / 4.0).ln()),
            (true, true, (*params.indel_prior / 2.0).ln()),
        ]
    } else {
        vec![(false, false, 0.0)]
    };
    let mut haplotypes = Vec::new();
    for state in states {
        ensure!(state.first.len() == sites.len(), "Incomplete phase hypothesis");
        let mut first = Haplotype { alleles: Vec::new(), indel: None };
        let mut second = first.clone();
        for (site, &base) in sites.iter().zip(&state.first) {
            let [a, b] = site.alleles.as_slice() else {
                color_eyre::eyre::bail!("Phasing requires two called alleles")
            };
            ensure!(base == *a || base == *b, "Phase allele outside called genotype");
            first.alleles.push((site.pos, base));
            second.alleles.push((site.pos, if base == *a { *b } else { *a }));
        }
        for &(a, b, _) in &configurations {
            let mut first = first.clone();
            let mut second = second.clone();
            first.indel = a.then(|| indel.clone()).flatten();
            second.indel = b.then(|| indel.clone()).flatten();
            first.alleles.sort_by_key(|(pos, _)| *pos);
            second.alleles.sort_by_key(|(pos, _)| *pos);
            haplotypes.extend([first, second]);
        }
    }
    let loci = super::scored_loci(reference, sites);
    let start = sites
        .iter()
        .map(|site| site.pos)
        .chain(indel.iter().map(|indel| indel.start))
        .min()
        .ok_or_eyre("Missing phase anchor")?
        .saturating_sub(u32::from(params.realignment_flank))
        .max(reference.start);
    let reference_end = reference
        .start
        .checked_add(u32::try_from(reference.bases.len())?)
        .ok_or_eyre("Reference coordinate overflow")?;
    let indel_end = indel
        .as_ref()
        .map(|indel| -> Result<u32> {
            indel
                .start
                .checked_add(u32::try_from(indel.deleted.len().max(1))?)
                .ok_or_eyre("Indel coordinate overflow")
        })
        .transpose()?;
    let end = sites
        .iter()
        .map(|site| site.pos.saturating_add(1))
        .chain(indel_end)
        .max()
        .ok_or_eyre("Missing phase end")?
        .saturating_add(u32::from(params.realignment_flank))
        .min(reference_end);
    let context = methylation::Context::new(reference, &haplotypes, &loci, Some(start..end))?;
    let local: Vec<_> = fragments
        .iter()
        .map(|fragment| {
            let selected: Vec<_> = sites
                .iter()
                .filter(|site| {
                    fragment
                        .reads
                        .iter()
                        .any(|read| read.covers(site.pos, LikelihoodEngine::PairHmm))
                })
                .cloned()
                .collect();
            ensure!(!selected.is_empty(), "Fragment does not cover phase block");
            pair_hmm::Prepared::new(reference, &selected, indel.as_ref(), &[*fragment], params)
                .map(pair_hmm::Prepared::localize)
        })
        .collect::<Result<_>>()?;
    let matrix = |rates: Option<&methylation::Rates>,
                  neutral: bool,
                  scorer: &mut pair_hmm::Scorer|
     -> Result<Vec<Vec<f64>>> {
        enum Slot {
            Planned(usize),
            Value(f64),
        }
        // Fragments are planned independently, on the pool when there is one;
        // the plans run as one batch afterwards, so the order never matters.
        let planned = fragments
            .par_iter()
            .zip(local.par_iter())
            .enumerate()
            .map(|(fragment_index, (fragment, prepared))| {
                let model = if neutral {
                    methylation::Scoring::Neutral(&context)
                } else if let Some(rates) = rates {
                    methylation::Scoring::Shrinkage(rates, fragment_index)
                } else {
                    methylation::Scoring::Ordinary
                };
                let mut cache = BTreeMap::<(Vec<(u32, u8)>, bool), usize>::new();
                let mut distinct: Vec<Slot> = Vec::new();
                let mut plans = Vec::new();
                let mut by_haplotype = Vec::with_capacity(haplotypes.len());
                for (index, haplotype) in haplotypes.iter().enumerate() {
                    let key = (prepared.signature(haplotype, params), haplotype.indel.is_some());
                    let slot = match cache.get(&key) {
                        Some(&slot) => slot,
                        None => {
                            let template =
                                prepared.template(reference, haplotype, &loci, params)?;
                            let slot =
                                match prepared.plan(0, fragment, &template, params, &model)? {
                                    Some(plan) => {
                                        plans.push(plan);
                                        Slot::Planned(plans.len() - 1)
                                    }
                                    None => Slot::Value(
                                        super::likelihood::fragment_likelihood_with_model(
                                            fragment,
                                            haplotype,
                                            reference,
                                            &loci,
                                            params,
                                            &model,
                                            context
                                                .haplotypes
                                                .get(index)
                                                .ok_or_eyre("Missing phase CpG context")?,
                                        )?,
                                    ),
                                };
                            distinct.push(slot);
                            cache.insert(key, distinct.len() - 1);
                            distinct.len() - 1
                        }
                    };
                    by_haplotype.push(slot);
                }
                Ok((distinct, by_haplotype, plans))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut plans = Vec::new();
        let mut slots = Vec::with_capacity(planned.len());
        for (distinct, by_haplotype, fragment_plans) in planned {
            let offset = plans.len();
            plans.extend(fragment_plans);
            slots.push((offset, distinct, by_haplotype));
        }
        // Shrinkage rates are refitted at every beam step, so only the other
        // models' forward results are worth remembering across calls.
        let scored = scorer.run(&plans, rates.is_none())?;
        let mut values = vec![vec![0.0; fragments.len()]; haplotypes.len()];
        for (fragment_index, (offset, distinct, by_haplotype)) in slots.iter().enumerate() {
            for (index, &slot) in by_haplotype.iter().enumerate() {
                let value = match distinct.get(slot).ok_or_eyre("Missing phase likelihood slot")? {
                    Slot::Planned(plan) => {
                        *scored.get(offset + plan).ok_or_eyre("Missing planned phase likelihood")?
                    }
                    Slot::Value(value) => *value,
                };
                *values
                    .get_mut(index)
                    .and_then(|row| row.get_mut(fragment_index))
                    .ok_or_eyre("Invalid phase likelihood index")? = value;
            }
        }
        Ok(values)
    };
    let joint = |values: &[Vec<f64>]| -> Result<Vec<(usize, usize, f64)>> {
        (0..states.len() * configurations.len())
            .map(|index| {
                let first = index * 2;
                let second = first + 1;
                let a = values.get(first).ok_or_eyre("Missing first chromosome")?;
                let b = values.get(second).ok_or_eyre("Missing second chromosome")?;
                Ok((
                    first,
                    second,
                    a.iter()
                        .zip(b)
                        .map(|(&a, &b)| log_add(a, b) - std::f64::consts::LN_2)
                        .sum::<f64>()
                        + configurations
                            .get(index % configurations.len())
                            .ok_or_eyre("Missing indel phase prior")?
                            .2,
                ))
            })
            .collect()
    };
    let values = if params.non_taps {
        matrix(None, false, scorer)?
    } else {
        let pilot = matrix(None, true, scorer)?;
        let rates = methylation::estimate_with_pairs(
            &context,
            &haplotypes,
            fragments,
            &pilot,
            params,
            &joint(&pilot)?,
        )?;
        matrix(Some(&rates), false, scorer)?
    };
    let pairs = joint(&values)?;
    let normalizer = log_sum(pairs.iter().map(|(_, _, value)| *value));
    let probabilities: Vec<_> = pairs
        .chunks(configurations.len())
        .map(|pairs| pairs.iter().map(|(_, _, value)| (value - normalizer).exp()).sum::<f64>())
        .collect();
    ensure!(probabilities.iter().all(|value| value.is_finite()), "Invalid phase posterior");
    Ok(Posterior { probabilities })
}

pub(crate) struct Linkage {
    pub edges: BTreeMap<(usize, usize), usize>,
    pub counts: Vec<[usize; 2]>,
    /// Per site, ascending, the fragments with an observation there on any
    /// read: what a block's scoring selects from.
    #[cfg_attr(
        not(feature = "experimental-seqair"),
        expect(dead_code, reason = "Read by phasing, which only the seqair backend runs")
    )]
    pub covering: Vec<Vec<usize>>,
}

pub(crate) fn links(
    sites: &[(u32, [Base; 2], Option<crate::metrics::methylation::CpgSide>)],
    fragments: &[Fragment],
) -> Result<Linkage> {
    use super::likelihood::{BaseMeasurement, consensus};
    let mut result = BTreeMap::new();
    let mut counts = vec![[0usize; 2]; sites.len()];
    let mut covering = vec![Vec::new(); sites.len()];
    for (fragment_index, fragment) in fragments.iter().enumerate() {
        let mut reached = Vec::new();
        for (index, (pos, bases, methylable)) in sites.iter().enumerate() {
            let mut observes = false;
            let observed =
                consensus(fragment.reads.iter().enumerate().filter_map(|(read_index, read)| {
                    let obs = read
                        .observations
                        .binary_search_by_key(pos, |obs| obs.pos)
                        .ok()
                        .and_then(|index| read.observations.get(index))?;
                    observes = true;
                    Some(BaseMeasurement {
                        base: obs.base,
                        quality: obs.quality,
                        mapq: read.mapq,
                        read_index,
                    })
                }))?;
            if observes && let Some(covering) = covering.get_mut(index) {
                covering.push(fragment_index);
            }
            let Some(source) = observed.source.and_then(|index| fragment.reads.get(index)) else {
                continue;
            };
            if !observed.quality.is_some_and(|quality| quality >= 20)
                || !source.mapq.is_some_and(|quality| quality >= 20)
                || !bases.contains(&observed.base)
            {
                continue;
            }
            if methylable.is_some_and(|side| {
                observed.base == side.mod_base()
                    && source.strand.ok().is_none_or(|strand| strand == side.strand())
            }) {
                continue;
            }
            let allele = usize::from(
                observed.base != bases.first().copied().ok_or_eyre("Missing linking genotype")?,
            );
            *counts
                .get_mut(index)
                .and_then(|counts| counts.get_mut(allele))
                .ok_or_eyre("Invalid allele support index")? += 1;
            reached.push(index);
        }
        for (offset, &a) in reached.iter().enumerate() {
            for &b in reached.iter().skip(offset + 1) {
                *result.entry((a, b)).or_insert(0) += 1;
            }
        }
    }
    Ok(Linkage { edges: result, counts, covering })
}

#[cfg(feature = "experimental-seqair")]
impl Linkage {
    /// The fragments observing any of these sites, ascending: the order the
    /// fragment slice has them in.
    pub(crate) fn fragments_covering(&self, sites: impl IntoIterator<Item = usize>) -> Vec<usize> {
        let mut selected: Vec<usize> = sites
            .into_iter()
            .filter_map(|site| self.covering.get(site))
            .flatten()
            .copied()
            .collect();
        selected.sort_unstable();
        selected.dedup();
        selected
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::haplotypes::{
        Observation, ReadEvidence, SequenceEvidence, SequenceObservation, SequencePosition,
    };
    use seqair_types::Strand;

    fn score(
        reference: &Reference,
        sites: &[Site],
        states: &[State],
        fragments: &[&Fragment],
        params: &LikelihoodParams,
    ) -> Result<Posterior> {
        let mut scorer = pair_hmm::Scorer::new(params)?;
        super::score(reference, sites, states, fragments, params, &mut scorer)
    }

    fn fixture() -> (Reference, Vec<Site>, Vec<State>) {
        let mut random = 0x7134_3afe_u32;
        let mut bases: Vec<_> = (0..60)
            .map(|_| {
                random = random.wrapping_mul(1664525).wrapping_add(1013904223);
                Base::KNOWN.get((random >> 30) as usize).copied().expect("DNA")
            })
            .collect();
        *bases.get_mut(20).expect("site") = Base::A;
        *bases.get_mut(35).expect("site") = Base::A;
        let reference = Reference { start: 1000, bases };
        let sites = [1020, 1035]
            .into_iter()
            .map(|pos| Site::new(pos, Base::A, [Base::G]).expect("site"))
            .collect();
        let states =
            vec![State { first: vec![Base::A, Base::A] }, State { first: vec![Base::A, Base::G] }];
        (reference, sites, states)
    }

    fn read(reference: &Reference, edits: &[(u32, Base)], start: u32, end: u32) -> ReadEvidence {
        let observations: Vec<_> = (start..end)
            .map(|pos| Observation {
                pos,
                base: edits
                    .iter()
                    .find(|(position, _)| *position == pos)
                    .map(|(_, base)| *base)
                    .unwrap_or_else(|| reference.base_at(pos).expect("reference")),
                quality: Some(40),
            })
            .collect();
        let sequence = SequenceEvidence {
            bases: observations
                .iter()
                .map(|obs| SequenceObservation {
                    position: SequencePosition::reference(obs.pos),
                    base: obs.base,
                    quality: obs.quality,
                })
                .collect(),
            reference_skip: false,
            indels: Vec::new(),
        };
        ReadEvidence {
            annotations: None,
            observations,
            sequence: Some(sequence),
            mapq: Some(60),
            strand: Strand::OT,
        }
    }

    #[test]
    fn a_joint_fragment_posterior_resolves_cis_and_trans() -> Result<()> {
        let (reference, sites, states) = fixture();
        let reads: Vec<_> = (0..4)
            .flat_map(|_| {
                [
                    Fragment { reads: vec![read(&reference, &[], 1005, 1050)] },
                    Fragment {
                        reads: vec![read(
                            &reference,
                            &[(1020, Base::G), (1035, Base::G)],
                            1005,
                            1050,
                        )],
                    },
                ]
            })
            .collect();
        let mut params = parameters()?;
        params.non_taps = true;
        let evidence =
            score(&reference, &sites, &states, &reads.iter().collect::<Vec<_>>(), &params)?;
        assert!(evidence.probabilities.first().copied().expect("cis") > 0.999);
        assert!((evidence.probabilities.iter().sum::<f64>() - 1.0).abs() < 1e-10);
        Ok(())
    }

    #[test]
    fn discordant_overlapping_mates_remove_the_phase_preference() -> Result<()> {
        let (reference, sites, states) = fixture();
        let first = read(&reference, &[], 1005, 1050);
        let conflicting = read(&reference, &[(1020, Base::G), (1035, Base::G)], 1005, 1050);
        let concordant = Fragment { reads: vec![first.clone(), first.clone()] };
        let discordant = Fragment { reads: vec![first, conflicting] };
        let mut params = parameters()?;
        params.non_taps = true;
        let good = score(&reference, &sites, &states, &[&concordant], &params)?;
        let bad = score(&reference, &sites, &states, &[&discordant], &params)?;
        assert!(good.probabilities.first().copied().expect("cis") > 0.99);
        assert!((bad.probabilities.first().copied().expect("cis") - 0.5).abs() < 1e-8);
        let linkage = links(
            &[(1020, [Base::A, Base::G], None), (1035, [Base::A, Base::G], None)],
            &[discordant],
        )?;
        assert!(linkage.edges.is_empty());
        assert_eq!(linkage.counts, vec![[0, 0], [0, 0]]);
        Ok(())
    }

    #[test]
    fn separate_single_site_molecules_do_not_identify_a_phase() -> Result<()> {
        let reference = Reference { start: 1000, bases: vec![Base::A; 300] };
        let sites =
            vec![Site::new(1020, Base::A, [Base::G])?, Site::new(1200, Base::A, [Base::G])?];
        let states =
            vec![State { first: vec![Base::A, Base::A] }, State { first: vec![Base::A, Base::G] }];
        let fragments: Vec<_> = [1020, 1200]
            .into_iter()
            .flat_map(|pos| {
                [Base::A, Base::G].map(|base| Fragment {
                    reads: vec![read(&reference, &[(pos, base)], pos, pos + 1)],
                })
            })
            .collect();
        let mut params = parameters()?;
        params.non_taps = true;
        let evidence =
            score(&reference, &sites, &states, &fragments.iter().collect::<Vec<_>>(), &params)?;
        assert!((evidence.probabilities.first().copied().expect("phase") - 0.5).abs() < 1e-10);
        Ok(())
    }

    #[test]
    fn scoring_is_independent_of_evidence_outside_unmasked_reference() -> Result<()> {
        let (full, sites, states) = fixture();
        let reference = Reference {
            start: 1018,
            bases: full.bases.get(18..37).expect("unmasked interval").to_vec(),
        };
        let fragments: Vec<_> =
            (0..2).map(|_| Fragment { reads: vec![read(&full, &[], 1005, 1050)] }).collect();
        let mut changed = fragments.clone();
        for fragment in &mut changed {
            for read in &mut fragment.reads {
                for observation in &mut read.observations {
                    if !(1018..1037).contains(&observation.pos) {
                        observation.base = Base::T;
                    }
                }
                let sequence = read.sequence.as_mut().expect("full sequence");
                for observation in &mut sequence.bases {
                    if !(1018..1037).contains(&observation.position.pos) {
                        observation.base = Base::T;
                    }
                }
                sequence.indels.push(super::super::Indel {
                    start: 1017,
                    deleted: vec![Base::A; 3],
                    inserted: Vec::new(),
                });
            }
        }
        let params = parameters()?;
        let a = score(&reference, &sites, &states, &fragments.iter().collect::<Vec<_>>(), &params)?;
        let b = score(&reference, &sites, &states, &changed.iter().collect::<Vec<_>>(), &params)?;
        assert_eq!(a.probabilities, b.probabilities);
        // Reference skips force the aligned fallback, which must respect the same bounds.
        for fragment in &mut changed {
            for read in &mut fragment.reads {
                read.sequence.as_mut().expect("full sequence").reference_skip = true;
            }
        }
        let mut unchanged_fallback = fragments;
        for fragment in &mut unchanged_fallback {
            for read in &mut fragment.reads {
                read.sequence.as_mut().expect("full sequence").reference_skip = true;
            }
        }
        let a = score(
            &reference,
            &sites,
            &states,
            &unchanged_fallback.iter().collect::<Vec<_>>(),
            &params,
        )?;
        let b = score(&reference, &sites, &states, &changed.iter().collect::<Vec<_>>(), &params)?;
        assert_eq!(a.probabilities, b.probabilities);
        Ok(())
    }
}
