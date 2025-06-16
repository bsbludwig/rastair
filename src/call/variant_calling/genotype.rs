// Adapted from rastair1, AGPL-3.0-only, (C) Benjamin Schuster-Boeckler
// source: https://bitbucket.org/bsblabludwig/rastair/src/306bf046f14c64992c06b9000c60113e35a0f766/src/operations/count_variants/mod.rs#lines-394

use crate::{
    call::{variant_calling::ErrorModel, variants::VariantCandidatePileup},
    utils::{Base, Strand},
};
use color_eyre::eyre::{Result, ensure};
use probability::prelude::{Binomial, Discrete as _, Distribution as _};
use rastair2_vcf::standard_fields::GenotypeAllele;
use tracing::{debug, instrument, trace};

impl VariantCandidatePileup {
    pub fn estimate_genotype(&self, error_model: ErrorModel) -> Option<EstimatedGenotype> {
        let (nosnp, snp) = if self.reference_base == Base::C {
            (
                self.bases.iter().filter(|b| b.base == Base::C && b.strand == Strand::OB).count(),
                self.bases.iter().filter(|b| b.base == Base::T && b.strand == Strand::OB).count(),
            )
        } else if self.reference_base == Base::G {
            (
                self.bases.iter().filter(|b| b.base == Base::G && b.strand == Strand::OT).count(),
                self.bases.iter().filter(|b| b.base == Base::A && b.strand == Strand::OT).count(),
            )
        } else {
            return None;
        };
        match EstimatedGenotype::calculate(nosnp, snp, error_model) {
            Ok(gt) => Some(gt),
            Err(error) => {
                debug!(%error, "Failed to calculate genotype");
                None
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GenotypeTag {
    /// Homozygous reference (CC)
    #[default]
    CC,
    /// Heterozygous (CT)
    CT,
    /// Homozygous alternative (TT)
    TT,
}

impl From<GenotypeTag> for [GenotypeAllele; 2] {
    fn from(value: GenotypeTag) -> Self {
        match value {
            GenotypeTag::CC => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)],
            GenotypeTag::CT => [GenotypeAllele::Phased(0), GenotypeAllele::Phased(1)],
            GenotypeTag::TT => [GenotypeAllele::Phased(1), GenotypeAllele::Phased(1)],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EstimatedGenotype {
    pub genotype: GenotypeTag,
    pub likelihood: f64,
    pub confidence: f64,
}

impl EstimatedGenotype {
    #[allow(clippy::cast_possible_truncation)]
    #[instrument(level = "trace", name = "calculate_genotype")]
    pub fn calculate(ref_count: usize, alt_count: usize, error_model: ErrorModel) -> Result<Self> {
        let error_rate = error_model.error_rate();
        ensure!(
            ref_count > 0 || alt_count > 0,
            "No ref or alt read counts, cannot compute likelihood"
        );
        ensure!(error_rate > f64::MIN, "Error rate too small, cannot calculate likelihood");

        // This is a simple estimate of genotype, based on the following consideration:
        // A site is either het or hom, where hom could be CC or TT.
        // If alt_count > ref_count, the latter is more likely, otherwise the former.

        // First, I calculate the likelihood to observe this many alt_reads
        // under the assumption that ref and alt are equally likely, ie this is a het position.
        // TODO This assumes a simple diploid sample with no purity issues. For
        // cancer samples, we could make this a setting to allow for different cancer fraction?

        let mut binom = Binomial::new(ref_count + alt_count, 0.5); // 0.5 because a het position
        let p_het = binom.mass(alt_count);
        let p_het_max = binom.mass(((alt_count + ref_count) as f64 / 2.0).round() as usize);

        // Then, I calculate the probability that this many or more alt_count/ref_count reads
        // are observed by error, assuming independence of reads and errors.
        binom = Binomial::new(ref_count + alt_count, error_rate);

        if ref_count >= alt_count {
            let p_hom = binom.mass(alt_count) + (1.0 - binom.distribution(alt_count as f64));

            if p_het < p_hom {
                trace!("Assuming CC: ({ref_count} vs {alt_count}) -> ({p_het:.5} < {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::CC,
                    likelihood: p_hom,
                    confidence: (p_hom - p_het) / p_hom,
                })
            } else {
                trace!("Assuming CT: ({ref_count} vs {alt_count}) -> ({p_het:.5} >= {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::CT,
                    likelihood: p_het / p_het_max,
                    confidence: (p_het - p_hom) / p_het,
                })
            }
        } else {
            let p_hom = binom.mass(ref_count) + (1.0 - binom.distribution(ref_count as f64));
            if p_het < p_hom {
                trace!("Assuming TT: ({ref_count} vs {alt_count}) -> ({p_het:.5} < {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::TT,
                    likelihood: p_hom,
                    confidence: (p_hom - p_het) / p_hom,
                })
            } else {
                trace!("Assuming TC: ({ref_count} vs {alt_count}) -> ({p_het:.5} >= {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::CT,
                    likelihood: p_het / p_het_max,
                    confidence: (p_het - p_hom) / p_het,
                })
            }
        }
    }
}
