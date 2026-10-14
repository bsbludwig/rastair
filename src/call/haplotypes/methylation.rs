use super::likelihood::{
    BaseMeasurement, consensus, log_add, log_sum, orientations_from, quadrature,
    retained_pair_priors, sequencing_emission,
};
use super::{Fragment, Haplotype, LikelihoodParams, Reference, SequencePosition, Site};
use crate::metrics::methylation::CpgSide;
use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, Strand};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct Cpg {
    pub c: SequencePosition,
    pub g: SequencePosition,
}

impl Cpg {
    pub fn contains(self, position: SequencePosition) -> bool {
        self.c == position || self.g == position
    }
}

pub(super) struct Context {
    pub haplotypes: Vec<BTreeSet<Cpg>>,
    ot: BTreeSet<SequencePosition>,
    ob: BTreeSet<SequencePosition>,
    calibration: CalibrationRegion,
}

enum CalibrationRegion {
    Aligned(BTreeSet<SequencePosition>),
    Realigned(std::ops::Range<u32>),
}

impl Context {
    pub fn new(
        reference: &Reference,
        haplotypes: &[Haplotype],
        loci: &[u32],
        realignment_window: Option<std::ops::Range<u32>>,
    ) -> Result<Self> {
        let reference_end = reference
            .start
            .checked_add(u32::try_from(reference.bases.len())?)
            .ok_or_eyre("Reference coordinate overflow")?;
        let mut start = loci.first().copied().ok_or_eyre("No scoring loci")?;
        let mut end = loci
            .last()
            .copied()
            .ok_or_eyre("No scoring loci")?
            .checked_add(1)
            .ok_or_eyre("Scoring coordinate overflow")?;
        let calibration = if let Some(window) = realignment_window {
            start = window.start;
            end = window.end;
            CalibrationRegion::Realigned(window)
        } else {
            CalibrationRegion::Aligned(
                loci.iter().copied().map(SequencePosition::reference).collect(),
            )
        };
        let start = start.saturating_sub(1).max(reference.start);
        let end = end.saturating_add(1).min(reference_end);
        let cpgs = haplotypes
            .iter()
            .map(|haplotype| super::pair_hmm::cpgs(reference, haplotype, start, end))
            .collect::<Result<Vec<_>>>()?;
        let ot = cpgs.iter().flatten().map(|cpg| cpg.c).collect();
        let ob = cpgs.iter().flatten().map(|cpg| cpg.g).collect();
        Ok(Self { haplotypes: cpgs, ot, ob, calibration })
    }

    fn accepts_observation(&self, position: SequencePosition) -> bool {
        match &self.calibration {
            CalibrationRegion::Aligned(positions) => positions.contains(&position),
            CalibrationRegion::Realigned(range) => range.contains(&position.pos),
        }
    }

    fn side(&self, position: SequencePosition, strand: Strand) -> Option<CpgSide> {
        match strand {
            Strand::OT if self.ot.contains(&position) => Some(CpgSide::C),
            Strand::OB if self.ob.contains(&position) => Some(CpgSide::G),
            _ => None,
        }
    }
}

pub(super) enum Scoring<'a> {
    #[cfg(test)]
    Fixed(f64),
    Ordinary,
    Neutral(&'a Context),
    Shrinkage(&'a Rates, usize),
}

impl Scoring<'_> {
    pub fn neutral_bases(
        &self,
        position: SequencePosition,
        strand: Strand,
        template: Base,
    ) -> Option<(Base, Base)> {
        let Self::Neutral(context) = self else { return None };
        let side = context.side(position, strand)?;
        (template == side.unmod_base() || template == side.mod_base())
            .then_some((side.unmod_base(), side.mod_base()))
    }

    pub fn rate(&self, cpg: Option<&Cpg>) -> f64 {
        match self {
            #[cfg(test)]
            Self::Fixed(rate) => *rate,
            Self::Shrinkage(rates, index) => cpg
                .and_then(|cpg| rates.values.get(cpg))
                .and_then(|values| values.get(*index))
                .copied()
                .unwrap_or(0.5),
            _ => 0.5,
        }
    }
}

pub(super) struct Rates {
    values: BTreeMap<Cpg, Vec<f64>>,
}

// Conversion observations do not enter the pilot genotype likelihoods. Fractional
// chromosome assignments retain uncertainty instead of committing to a pilot call.
pub(super) fn estimate(
    context: &Context,
    haplotypes: &[Haplotype],
    sites: &[Site],
    fragments: &[&Fragment],
    pilot: &[Vec<f64>],
    params: &LikelihoodParams,
    has_indel: bool,
) -> Result<Rates> {
    ensure!(
        context.haplotypes.len() == haplotypes.len()
            && pilot.len() == haplotypes.len()
            && pilot.iter().all(|values| values.len() == fragments.len()),
        "Pilot likelihood dimensions disagree"
    );
    let mut pairs = Vec::new();
    let priors = retained_pair_priors(haplotypes, sites, params, has_indel)?;
    for (a, _) in haplotypes.iter().enumerate() {
        for (b, _) in haplotypes.iter().enumerate().skip(a) {
            let likelihoods_a = pilot.get(a).ok_or_eyre("Missing pilot haplotype")?;
            let likelihoods_b = pilot.get(b).ok_or_eyre("Missing pilot haplotype")?;
            let joint = likelihoods_a
                .iter()
                .zip(likelihoods_b)
                .map(|(&a, &b)| log_add(a, b) - std::f64::consts::LN_2)
                .sum::<f64>()
                + *priors.get(pairs.len()).ok_or_eyre("Missing pilot pair prior")?;
            pairs.push((a, b, joint));
        }
    }
    estimate_with_pairs(context, haplotypes, fragments, pilot, params, &pairs)
}

pub(super) fn estimate_with_pairs(
    context: &Context,
    haplotypes: &[Haplotype],
    fragments: &[&Fragment],
    pilot: &[Vec<f64>],
    params: &LikelihoodParams,
    pairs: &[(usize, usize, f64)],
) -> Result<Rates> {
    ensure!(!pairs.is_empty(), "No eligible diploid pairs for CpG estimation");
    let normalizer = log_sum(pairs.iter().map(|&(_, _, joint)| joint));
    let mut assignments = vec![vec![0.0; fragments.len()]; haplotypes.len()];
    for &(a, b, joint) in pairs {
        let posterior = (joint - normalizer).exp();
        for index in 0..fragments.len() {
            let first = *pilot
                .get(a)
                .and_then(|values| values.get(index))
                .ok_or_eyre("Missing pilot fragment")?;
            let second = *pilot
                .get(b)
                .and_then(|values| values.get(index))
                .ok_or_eyre("Missing pilot fragment")?;
            let responsibility = (first - log_add(first, second)).exp();
            *assignments
                .get_mut(a)
                .and_then(|values| values.get_mut(index))
                .ok_or_eyre("Missing assignment")? += posterior * responsibility;
            *assignments
                .get_mut(b)
                .and_then(|values| values.get_mut(index))
                .ok_or_eyre("Missing assignment")? += posterior * (1.0 - responsibility);
        }
    }
    let nodes = quadrature(params.quadrature_order)?;
    let prior: Vec<_> = nodes
        .iter()
        .map(|node| {
            node.log_weight
                + (params.methylation_prior_strength / 2.0 - 1.0)
                    * (node.methylation.ln() + (-node.methylation).ln_1p())
        })
        .collect();
    let mut values = BTreeMap::new();
    for cpg in context.haplotypes.iter().flatten().copied().collect::<BTreeSet<_>>() {
        let mut total = prior.clone();
        let mut contributions = Vec::with_capacity(fragments.len());
        for (index, fragment) in fragments.iter().enumerate() {
            let mass = context
                .haplotypes
                .iter()
                .zip(&assignments)
                .filter(|(cpgs, _)| cpgs.contains(&cpg))
                .map(|(_, assignments)| assignments.get(index).copied().unwrap_or(0.0))
                .sum::<f64>();
            // Ambiguous CpG/non-CpG chromosome assignments cannot distinguish
            // conversion from a genetic T/A, so they must not fit a CpG rate.
            let contribution = if mass >= *params.methylation_min_cpg_probability {
                observation(fragment, cpg, context)?
            } else {
                None
            }
            .map(|observation| {
                let observed = observation.base;
                nodes
                    .iter()
                    .map(|node| {
                        let rate = params.conversion(Some(node.methylation));
                        let emission = (1.0 - rate)
                            * sequencing_emission(
                                observed.base,
                                observation.side.unmod_base(),
                                observed.quality,
                            )
                            + rate
                                * sequencing_emission(
                                    observed.base,
                                    observation.side.mod_base(),
                                    observed.quality,
                                );
                        let error =
                            observation.mapq.map_or(1.0, |q| 10.0_f64.powf(-f64::from(q) / 10.0));
                        mass * ((1.0 - error) * emission + error * 0.25).ln()
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_else(|| vec![0.0; nodes.len()]);
            for (total, contribution) in total.iter_mut().zip(&contribution) {
                *total += contribution;
            }
            contributions.push(contribution);
        }
        let rates = contributions
            .iter()
            .map(|own| {
                // Leaving this molecule out prevents its converted base from fitting
                // the rate used to score that very same observation.
                let normalizer = log_sum(total.iter().zip(own).map(|(total, own)| total - own));
                nodes
                    .iter()
                    .zip(&total)
                    .zip(own)
                    .map(|((node, total), own)| node.methylation * (total - own - normalizer).exp())
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        ensure!(
            rates.iter().all(|rate| rate.is_finite() && (0.0..=1.0).contains(rate)),
            "Invalid estimated methylation rate"
        );
        values.insert(cpg, rates);
    }
    Ok(Rates { values })
}

struct ConversionObservation {
    base: super::likelihood::MolecularBase,
    side: CpgSide,
    mapq: Option<u8>,
}

fn observation(
    fragment: &Fragment,
    cpg: Cpg,
    context: &Context,
) -> Result<Option<ConversionObservation>> {
    let present = |read: &super::ReadEvidence, position: SequencePosition| {
        if !context.accepts_observation(position) {
            return None;
        }
        if let Some(sequence) = &read.sequence {
            sequence
                .bases
                .binary_search_by_key(&position, |base| base.position)
                .ok()
                .and_then(|index| sequence.bases.get(index))
                .map(|base| (base.base, base.quality))
        } else if position.is_reference() {
            read.observations
                .binary_search_by_key(&position.pos, |base| base.pos)
                .ok()
                .and_then(|index| read.observations.get(index))
                .map(|base| (base.base, base.quality))
        } else {
            None
        }
    };
    let orientations = orientations_from(
        fragment
            .reads
            .iter()
            .filter(|read| present(read, cpg.c).is_some() || present(read, cpg.g).is_some())
            .map(|read| read.strand),
    );
    let [strand] = orientations else { return Ok(None) };
    let (position, side) = match strand {
        Strand::OT => (cpg.c, CpgSide::C),
        Strand::OB => (cpg.g, CpgSide::G),
        Strand::Unknown => return Ok(None),
    };
    let observed =
        consensus(fragment.reads.iter().enumerate().filter_map(|(read_index, read)| {
            present(read, position).map(|(base, quality)| BaseMeasurement {
                base,
                quality,
                mapq: read.mapq,
                read_index,
            })
        }))?;
    let Some(source) = observed.source else { return Ok(None) };
    let mapq = fragment.reads.get(source).ok_or_eyre("Missing consensus source")?.mapq;
    Ok(Some(ConversionObservation { base: observed, side, mapq }))
}

#[cfg(test)]
mod tests {
    use super::super::{LikelihoodEngine, Observation, ReadEvidence};
    use super::*;

    #[test]
    fn candidate_union_includes_inserted_and_deletion_created_cpgs() -> Result<()> {
        for (sequence, deleted, inserted, expected) in [
            (
                b"ATGA".as_slice(),
                vec![],
                vec![Base::C],
                Cpg {
                    c: SequencePosition::insertion(1002, 0),
                    g: SequencePosition::reference(1002),
                },
            ),
            (
                b"ACTTGA".as_slice(),
                vec![Base::T, Base::T],
                vec![],
                Cpg { c: SequencePosition::reference(1001), g: SequencePosition::reference(1004) },
            ),
        ] {
            let reference = Reference {
                start: 1000,
                bases: sequence.iter().copied().map(Base::from).collect(),
            };
            let haplotypes = [
                Haplotype { alleles: vec![], indel: None },
                Haplotype {
                    alleles: vec![],
                    indel: Some(super::super::Indel { start: 1002, deleted, inserted }),
                },
            ];
            let context = Context::new(&reference, &haplotypes, &[1001, expected.g.pos], None)?;
            assert!(context.haplotypes.first().ok_or_eyre("reference CpGs")?.is_empty());
            assert_eq!(
                context.haplotypes.get(1).ok_or_eyre("edited CpGs")?,
                &BTreeSet::from([expected])
            );
            assert_eq!(
                Scoring::Neutral(&context).neutral_bases(expected.g, Strand::OB, Base::A),
                Some((Base::G, Base::A))
            );
            assert_eq!(
                Scoring::Neutral(&context).neutral_bases(expected.c, Strand::OT, Base::T),
                Some((Base::C, Base::T))
            );
        }
        Ok(())
    }

    #[test]
    fn distant_indel_cpgs_are_neutralized_without_using_an_inactive_mate() -> Result<()> {
        let mut bases = vec![Base::A; 200];
        *bases.get_mut(180).ok_or_eyre("fixture G")? = Base::G;
        let reference = Reference { start: 1000, bases };
        let haplotypes = [
            Haplotype { alleles: vec![], indel: None },
            Haplotype {
                alleles: vec![],
                indel: Some(super::super::Indel {
                    start: 1180,
                    deleted: vec![],
                    inserted: vec![Base::C],
                }),
            },
        ];
        let mut params = super::super::tests::params(false);
        params.likelihood_engine = LikelihoodEngine::PairHmm;
        let sites = [Site::new(1020, Base::A, [Base::C])?];
        let prepared = super::super::pair_hmm::Prepared::new(
            &reference,
            &sites,
            haplotypes.get(1).and_then(|haplotype| haplotype.indel.as_ref()),
            &[],
            &params,
        )?;
        let context = Context::new(&reference, &haplotypes, &[1020], Some(prepared.interval()))?;
        assert_eq!(
            Scoring::Neutral(&context).neutral_bases(
                SequencePosition::reference(1180),
                Strand::OB,
                Base::A
            ),
            Some((Base::G, Base::A))
        );
        let mut edge = reference.clone();
        *edge.bases.get_mut(32).ok_or_eyre("edge C")? = Base::C;
        *edge.bases.get_mut(33).ok_or_eyre("edge G")? = Base::G;
        let haplotypes = [Haplotype { alleles: vec![], indel: None }];
        let prepared = super::super::pair_hmm::Prepared::new(&edge, &sites, None, &[], &params)?;
        let context = Context::new(&edge, &haplotypes, &[1020, 1021], Some(prepared.interval()))?;
        let cpg = context.haplotypes.first().and_then(BTreeSet::first).ok_or_eyre("edge CpG")?;
        let mate = Fragment {
            reads: vec![ReadEvidence {
                annotations: None,
                mapq: Some(60),
                strand: Strand::OB,
                sequence: None,
                observations: vec![Observation { pos: 1033, base: Base::A, quality: Some(40) }],
            }],
        };
        assert!(observation(&mate, *cpg, &context)?.is_none());
        Ok(())
    }

    fn molecule(first: Base, second: Base) -> Fragment {
        Fragment {
            reads: vec![ReadEvidence {
                annotations: None,
                strand: Strand::OT,
                mapq: Some(60),
                sequence: None,
                observations: vec![
                    Observation { pos: 1001, base: first, quality: Some(40) },
                    Observation { pos: 1005, base: second, quality: Some(40) },
                ],
            }],
        }
    }

    #[test]
    fn separate_cpg_rates_shrink_and_exclude_the_scored_molecule() -> Result<()> {
        let reference =
            Reference { start: 1000, bases: b"ACGAACGA".iter().copied().map(Base::from).collect() };
        let sites = [Site::new(1001, Base::C, [Base::T])?];
        let haplotypes = [Haplotype { alleles: vec![(1001, Base::C)], indel: None }];
        let params = super::super::tests::params(false);
        let context = Context::new(&reference, &haplotypes, &[1001, 1005], None)?;
        let cpgs = context.haplotypes.first().ok_or_eyre("cpgs")?;
        let first = cpgs.first().ok_or_eyre("first CpG")?;
        let second = cpgs.last().ok_or_eyre("second CpG")?;
        let singleton = molecule(Base::T, Base::C);
        let rates =
            estimate(&context, &haplotypes, &sites, &[&singleton], &[vec![0.0]], &params, false)?;
        assert!((Scoring::Shrinkage(&rates, 0).rate(Some(first)) - 0.5).abs() < 1e-12);
        let mut paired = singleton.clone();
        paired.reads.extend(singleton.reads.clone());
        let doubled =
            estimate(&context, &haplotypes, &sites, &[&paired], &[vec![0.0]], &params, false)?;
        assert_eq!(rates.values, doubled.values);
        let molecules = vec![singleton; 24];
        let fragments = molecules.iter().collect::<Vec<_>>();
        let rates =
            estimate(&context, &haplotypes, &sites, &fragments, &[vec![0.0; 24]], &params, false)?;
        assert!(Scoring::Shrinkage(&rates, 0).rate(Some(first)) > 0.8);
        assert!(Scoring::Shrinkage(&rates, 0).rate(Some(second)) < 0.2);
        let paired = molecules
            .iter()
            .cloned()
            .map(|mut fragment| {
                fragment.reads.extend(fragment.reads.clone());
                fragment
            })
            .collect::<Vec<_>>();
        let doubled = estimate(
            &context,
            &haplotypes,
            &sites,
            &paired.iter().collect::<Vec<_>>(),
            &[vec![0.0; 24]],
            &params,
            false,
        )?;
        assert_eq!(rates.values, doubled.values);
        let mut unknown = molecules;
        for read in unknown.iter_mut().flat_map(|fragment| &mut fragment.reads) {
            for base in &mut read.observations {
                base.quality = None;
            }
        }
        let rates = estimate(
            &context,
            &haplotypes,
            &sites,
            &unknown.iter().collect::<Vec<_>>(),
            &[vec![0.0; 24]],
            &params,
            false,
        )?;
        assert!((Scoring::Shrinkage(&rates, 0).rate(Some(first)) - 0.5).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn ambiguous_cpg_carriers_cannot_fit_conversion_rates() -> Result<()> {
        let reference =
            Reference { start: 1000, bases: b"ACGAACGA".iter().copied().map(Base::from).collect() };
        let sites = [Site::new(1001, Base::C, [Base::T])?];
        let haplotypes =
            [Base::C, Base::T].map(|base| Haplotype { alleles: vec![(1001, base)], indel: None });
        let mut params = super::super::tests::params(false);
        params.likelihood_engine = LikelihoodEngine::Aligned;
        let context = Context::new(&reference, &haplotypes, &[1001], None)?;
        let mut molecules = vec![molecule(Base::T, Base::Unknown); 16];
        for base in [Base::C, Base::T] {
            for _ in 0..8 {
                let mut fragment = molecule(base, Base::Unknown);
                for read in &mut fragment.reads {
                    read.strand = Strand::OB;
                }
                molecules.push(fragment);
            }
        }
        let pilot = haplotypes
            .iter()
            .map(|haplotype| {
                molecules
                    .iter()
                    .map(|fragment| {
                        let read = fragment.reads.first().ok_or_eyre("read")?;
                        let observation = read.observations.first().ok_or_eyre("base")?;
                        Ok(if read.strand == Strand::OT {
                            0.0
                        } else {
                            sequencing_emission(
                                observation.base,
                                haplotype.allele_at(1001).ok_or_eyre("allele")?,
                                Some(40),
                            )
                            .ln()
                        })
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let rates = estimate(
            &context,
            &haplotypes,
            &sites,
            &molecules.iter().collect::<Vec<_>>(),
            &pilot,
            &params,
            false,
        )?;
        let cpg = context.haplotypes.first().and_then(BTreeSet::first).ok_or_eyre("CpG")?;
        assert!((Scoring::Shrinkage(&rates, 0).rate(Some(cpg)) - 0.5).abs() < 1e-12);
        Ok(())
    }
}
