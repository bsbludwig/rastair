use super::{Haplotype, PairEvidence};
use crate::call::{
    phasing::ColumnPhase,
    pileup::indels::IndelAllele,
    variant_calling::{ErrorModel, EstimatedGenotype, GenotypeTag},
};
use crate::metrics::{AltCall, PileupMetrics};
use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, Phred, Probability};
use std::num::NonZeroU8;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum GenotypeModel {
    #[default]
    Counts,
    PairHmm,
}

impl GenotypeModel {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Counts => "counts",
            Self::PairHmm => "pairhmm",
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct DiploidPosterior<A> {
    alleles: Vec<A>,
    probabilities: Vec<Probability>,
    pruning_mass: f64,
}

impl<A> DiploidPosterior<A> {
    fn new(alleles: Vec<A>, probabilities: Vec<Probability>, pruning_mass: f64) -> Result<Self> {
        let count = alleles.len();
        ensure!(count > 0 && count <= 256, "Invalid molecular genotype allele count");
        ensure!(
            probabilities.len() == count * (count + 1) / 2,
            "Invalid diploid posterior dimensions"
        );
        let total: f64 = probabilities.iter().map(|value| **value).sum();
        ensure!((total - 1.0).abs() < 1e-8, "Molecular genotype posteriors do not normalize");
        ensure!((0.0..=1.0).contains(&pruning_mass), "Invalid genotype pruning mass");
        Ok(Self { alleles, probabilities, pruning_mass })
    }

    pub(super) fn entropy(&self) -> f64 {
        -self
            .probabilities
            .iter()
            .map(|p| **p)
            .filter(|p| *p > 0.0)
            .map(|p| p * p.ln())
            .sum::<f64>()
    }

    /// The forest decides whether an allele is present; the likelihood decides
    /// its dosage. A genotype carrying an admissible non-reference allele is
    /// chosen over homozygous reference whenever there is one, however the
    /// prior-weighted posterior ranks them: with the per-locus prior and a
    /// library whose homopolymer stutter is this common, real heterozygotes
    /// would otherwise be emitted as `0/0` under a passing ML score. Its
    /// confidence is conditional on that choice, the dosage error given
    /// presence. Disallowed ML alleles retain their probability mass as
    /// uncertainty either way: the forest is a gate on emitted alleles, not
    /// independent likelihood evidence.
    pub(super) fn select(
        &self,
        map: impl Fn(&A) -> Option<u8>,
    ) -> Result<Option<EstimatedGenotype>> {
        let mut reference: Option<f64> = None;
        let mut best: Option<(f64, u8, u8)> = None;
        let mut admissible = 0.0;
        for (higher, b) in self.alleles.iter().enumerate() {
            let Some(b) = map(b) else { continue };
            for (lower, a) in self.alleles.iter().take(higher + 1).enumerate() {
                let Some(a) = map(a) else { continue };
                let index = higher * (higher + 1) / 2 + lower;
                let probability =
                    **self.probabilities.get(index).ok_or_eyre("Missing diploid posterior")?;
                if a == 0 && b == 0 {
                    reference = Some(reference.map_or(probability, |p| p.max(probability)));
                    continue;
                }
                admissible += probability;
                if best.is_none_or(|(previous, _, _)| probability > previous) {
                    best = Some((probability, a.min(b), a.max(b)));
                }
            }
        }
        let (posterior, a, b) = match (best, reference) {
            (Some((probability, a, b)), _) => ((probability / admissible).clamp(0.0, 1.0), a, b),
            (None, Some(probability)) => (probability, 0, 0),
            (None, None) => return Ok(None),
        };
        let genotype = match (NonZeroU8::new(a), NonZeroU8::new(b)) {
            (None, None) => GenotypeTag::HomRef,
            (None, Some(b)) => GenotypeTag::ref_het(b),
            (Some(a), Some(b)) if a == b => GenotypeTag::hom_alt(a),
            (Some(a), Some(b)) => GenotypeTag::alt_het(a, b),
            (Some(_), None) => color_eyre::eyre::bail!("Unordered molecular genotype"),
        };
        // Search mass is an observed approximation diagnostic, not a calibrated
        // bound on all omitted hypotheses. It can lower, never raise, the GQ.
        let error = Probability::new((1.0 - posterior).max(self.pruning_mass).clamp(1e-10, 1.0))?;
        Ok(Some(EstimatedGenotype { genotype, likelihood: error, confidence: error }))
    }
}

impl DiploidPosterior<Base> {
    pub(super) fn pairs(
        haplotypes: &[Haplotype],
        pairs: &[PairEvidence],
        focal: u32,
        pruning_mass: f64,
    ) -> Result<Self> {
        let mut probabilities = [0.0; 10];
        for pair in pairs {
            let (a, b) = super::pair_bases(pair, haplotypes, focal)?;
            let index = super::genotype_index(a, b).ok_or_eyre("Unknown molecular SNP genotype")?;
            *probabilities.get_mut(index).ok_or_eyre("Invalid molecular SNP genotype")? +=
                *pair.posterior;
        }
        Self::new(
            Base::KNOWN.to_vec(),
            probabilities.into_iter().map(super::probability).collect::<Result<_>>()?,
            pruning_mass,
        )
    }
}

impl DiploidPosterior<Option<IndelAllele>> {
    #[cfg(test)]
    pub(super) fn indels(
        haplotypes: &[Haplotype],
        pairs: &[PairEvidence],
        pruning_mass: f64,
    ) -> Result<Self> {
        let starts = haplotypes
            .iter()
            .flat_map(|h| &h.indels)
            .map(|e| e.start)
            .collect::<std::collections::BTreeSet<_>>();
        ensure!(starts.len() <= 1, "Use position-specific marginalization for regional indels");
        Self::indels_at(haplotypes, pairs, starts.first().copied().unwrap_or(0), pruning_mass)
    }

    pub(super) fn indels_at(
        haplotypes: &[Haplotype],
        pairs: &[PairEvidence],
        start: u32,
        pruning_mass: f64,
    ) -> Result<Self> {
        let edits = haplotypes
            .iter()
            .flat_map(|haplotype| &haplotype.indels)
            .filter(|edit| edit.start == start)
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let mut alleles = vec![None];
        for edit in &edits {
            alleles.push(Some(if edit.inserted.is_empty() {
                IndelAllele::Deletion(edit.deleted.iter().copied().collect())
            } else {
                IndelAllele::Insertion(edit.inserted.iter().copied().collect())
            }));
        }
        let mut probabilities = vec![0.0; alleles.len() * (alleles.len() + 1) / 2];
        let index = |haplotype: usize| -> Result<usize> {
            match haplotypes
                .get(haplotype)
                .ok_or_eyre("Missing genotype haplotype")?
                .indels
                .iter()
                .find(|edit| edit.start == start)
            {
                None => Ok(0),
                Some(edit) => edits
                    .binary_search(&edit)
                    .map(|index| index + 1)
                    .map_err(|_| color_eyre::eyre::eyre!("Missing genotype indel")),
            }
        };
        for pair in pairs {
            let a = index(pair.first)?;
            let b = index(pair.second)?;
            let slot = a.max(b) * (a.max(b) + 1) / 2 + a.min(b);
            *probabilities.get_mut(slot).ok_or_eyre("Invalid joint indel genotype")? +=
                *pair.posterior;
        }
        Self::new(
            alleles,
            probabilities.into_iter().map(super::probability).collect::<Result<_>>()?,
            pruning_mass,
        )
    }
}

pub(crate) fn model(pileup: &PileupMetrics) -> GenotypeModel {
    if matches!(pileup.phase.as_deref(), Some(ColumnPhase::Genotyped(_))) {
        GenotypeModel::PairHmm
    } else {
        GenotypeModel::Counts
    }
}

pub(crate) fn assign_snp(
    pileup: &mut PileupMetrics,
    threshold: Option<Probability>,
    error: ErrorModel,
) -> Result<()> {
    let evidence = match pileup.phase.as_deref() {
        Some(ColumnPhase::Likelihood(evidence)) => evidence.snp_genotype.as_ref(),
        _ => None,
    };
    let estimate = evidence
        .map(|posterior| {
            posterior.select(|base| {
                if *base == pileup.reference_base {
                    return Some(0);
                }
                let (index, _) = pileup.alts.iter().enumerate().find(|(_, alt)| {
                    alt.base == *base
                        && matches!(alt.call, AltCall::RealVariant)
                        && threshold.is_none_or(|threshold| {
                            alt.filters.ml.is_some_and(|score| score >= threshold)
                        })
                })?;
                u8::try_from(index + 1).ok()
            })
        })
        .transpose()?
        .flatten();
    let molecular = estimate.is_some();
    pileup.pos_metrics.extended.genotype =
        estimate.or_else(|| pileup.estimate_genotype(threshold, error));
    if molecular {
        pileup.phase = Some(Box::new(ColumnPhase::Genotyped(None)));
    } else if matches!(pileup.phase.as_deref(), Some(ColumnPhase::Likelihood(_))) {
        pileup.phase = None;
    }
    Ok(())
}

pub(crate) fn assign_indels(
    pileup: &mut PileupMetrics,
    threshold: Option<Probability>,
) -> Result<()> {
    let Some(ColumnPhase::Likelihood(evidence)) = pileup.phase.as_deref() else { return Ok(()) };
    let Some(posterior) = &evidence.indel_genotype else { return Ok(()) };
    let Some(data) = pileup.indel_data.as_mut() else { return Ok(()) };
    let estimate = posterior.select(|allele| match allele {
        None => Some(0),
        Some(allele) => {
            let (index, _) = data.calls.iter().enumerate().find(|(_, call)| {
                call.allele == *allele
                    && !call.one_sided
                    && !call.ml.zip(threshold).is_some_and(|(score, threshold)| score < threshold)
            })?;
            u8::try_from(index + 1).ok()
        }
    })?;
    let Some(estimate) = estimate else { return Ok(()) };
    let copies = |index: u8| -> u8 {
        match estimate.genotype {
            GenotypeTag::HomRef => 0,
            GenotypeTag::RefHet(a) => u8::from(a.get() == index),
            GenotypeTag::HomAlt(a) => 2 * u8::from(a.get() == index),
            GenotypeTag::AltHet(a, b) => u8::from(a.get() == index || b.get() == index),
        }
    };
    let one = NonZeroU8::new(1).ok_or_eyre("Invalid genotype allele")?;
    for (index, call) in data.calls.iter_mut().enumerate() {
        let count = copies(u8::try_from(index + 1)?);
        call.genotype = match count {
            0 => GenotypeTag::HomRef,
            1 => GenotypeTag::ref_het(one),
            2 => GenotypeTag::hom_alt(one),
            _ => color_eyre::eyre::bail!("Invalid diploid allele dosage"),
        };
        call.quality = Phred::from(estimate.confidence);
        call.genotype_model = GenotypeModel::PairHmm;
    }
    if let GenotypeTag::AltHet(a, b) = estimate.genotype {
        let selected = [a, b].map(|index| data.calls.get(usize::from(index.get()) - 1).cloned());
        let [Some(mut first), Some(mut second)] = selected else {
            color_eyre::eyre::bail!("Missing called indel genotype")
        };
        let two = NonZeroU8::new(2).ok_or_eyre("Invalid second genotype allele")?;
        first.genotype = GenotypeTag::alt_het(one, two);
        second.genotype = first.genotype;
        data.calls = vec![first, second];
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn distribution(
        alleles: Vec<Base>,
        values: &[f64],
        pruning: f64,
    ) -> Result<DiploidPosterior<Base>> {
        DiploidPosterior::new(
            alleles,
            values.iter().copied().map(Probability::new).collect::<Result<_, _>>()?,
            pruning,
        )
    }

    #[test]
    fn joint_genotypes_keep_compound_hets_and_genotype_indices() -> Result<()> {
        let posterior = distribution(
            vec![Base::A, Base::C, Base::G],
            &[0.001, 0.001, 0.001, 0.001, 0.995, 0.001],
            0.0,
        )?;
        let selected = posterior
            .select(|base| {
                Some(match base {
                    Base::A => 0,
                    Base::C => 2,
                    _ => 1,
                })
            })?
            .ok_or_eyre("genotype")?;
        assert_eq!(
            selected.genotype,
            GenotypeTag::alt_het(
                NonZeroU8::new(1).ok_or_eyre("allele")?,
                NonZeroU8::new(2).ok_or_eyre("allele")?
            )
        );
        // Conditional on an admissible allele being present: the homozygous
        // reference mass is not part of the dosage decision.
        assert!((*selected.confidence - (1.0 - 0.995 / 0.999)).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn an_accepted_allele_is_never_vetoed_by_the_reference_genotype() -> Result<()> {
        let posterior = distribution(vec![Base::A, Base::C], &[0.99, 0.009, 0.001], 0.0)?;
        let selected =
            posterior.select(|base| Some(u8::from(*base == Base::C)))?.ok_or_eyre("genotype")?;
        assert_eq!(
            selected.genotype,
            GenotypeTag::ref_het(NonZeroU8::new(1).ok_or_eyre("allele")?)
        );
        assert!((*selected.confidence - 0.1).abs() < 1e-12);
        let rejected =
            posterior.select(|base| (*base == Base::A).then_some(0))?.ok_or_eyre("genotype")?;
        assert_eq!(rejected.genotype, GenotypeTag::HomRef);
        assert!((*rejected.confidence - 0.01).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn filtering_an_allele_never_renormalizes_its_evidence_into_confidence() -> Result<()> {
        let posterior = distribution(vec![Base::A, Base::C], &[0.005, 0.005, 0.99], 0.0)?;
        let selected =
            posterior.select(|base| (*base == Base::A).then_some(0))?.ok_or_eyre("genotype")?;
        assert_eq!(selected.genotype, GenotypeTag::HomRef);
        assert!((*selected.confidence - 0.995).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn genotype_quality_is_limited_by_search_uncertainty() -> Result<()> {
        let posterior = distribution(vec![Base::A, Base::C], &[0.00001, 0.99998, 0.00001], 0.1)?;
        let selected =
            posterior.select(|base| Some(u8::from(*base == Base::C)))?.ok_or_eyre("genotype")?;
        assert_eq!(Phred::from(selected.confidence).as_int(), 10);
        assert!(distribution(vec![Base::A, Base::C], &[0.1, 0.2, 0.2], 0.0).is_err());
        Ok(())
    }
}
