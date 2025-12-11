// Adapted from rastair1, (c) Benjamin Schuster-Boeckler
use crate::{
    call::variant_calling::ErrorModel,
    metrics::{AltCall, PileupMetrics},
    utils::{Base, Base::*, IntoF64 as _, Strand::*},
};
use color_eyre::eyre::{ContextCompat, Result, bail, ensure};
use probability::prelude::{Binomial, Discrete as _, Distribution as _};
use rastair_types::Probability;
use rastair_vcf::standard_fields::{Genotype, GenotypeAllele};
use std::num::{NonZeroI32, NonZeroU8};
use tracing::{instrument, trace};

impl PileupMetrics {
    /// Estimate genotype using ML scores and variant calls.
    ///
    /// This method implements the genotyping logic:
    /// - If ML score < `ml_threshold`: genotype is `0/0` (`HomRef`)
    /// - If ML score >= `ml_threshold`: use binomial model to determine `0/1`, `1/1`, or compound het (`1/2`)
    ///
    /// Only considers alts marked as [`AltCall::RealVariant`].
    pub fn estimate_genotype(
        &self,
        ml_threshold: Option<Probability>,
        error_model: ErrorModel,
    ) -> Option<EstimatedGenotype> {
        // Filter to only real variants
        let real_variant_alts: Vec<_> = self
            .alts
            .iter()
            .enumerate()
            .filter(|(_, alt)| matches!(alt.call, AltCall::RealVariant))
            .collect();

        if real_variant_alts.is_empty() {
            return None;
        }

        // If ML threshold is set, filter by ML score
        let passing_alts: Vec<_> = if let Some(threshold) = ml_threshold {
            real_variant_alts
                .into_iter()
                .filter(|(_, alt)| alt.filters.ml.is_some_and(|ml| ml >= threshold))
                .collect()
        } else {
            real_variant_alts
        };

        // If no alts pass ML threshold, return 0/0 with confidence based on distance from threshold
        if passing_alts.is_empty() {
            if let Some(threshold) = ml_threshold {
                let max_ml = self
                    .alts
                    .iter()
                    .filter_map(|alt| alt.filters.ml)
                    .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or(Probability::ZERO);

                let confidence = (*threshold - *max_ml) / *threshold;
                return Some(EstimatedGenotype {
                    genotype: GenotypeTag::HomRef,
                    likelihood: Probability::ONE,
                    confidence: Probability::new(confidence).ok()?,
                });
            } else {
                return None;
            }
        }

        // Calculate genotype for each passing alt using binomial model
        // Note: We compare each alt independently against ref (Option A).
        // Alternative approaches: compare alt1 vs alt2 for compound het (Option B),
        // or compare each alt against total depth minus its count (Option C).
        let mut alt_genotypes = Vec::new();
        for (alt_idx, alt) in &passing_alts {
            let alt_base = alt.base;
            let (ref_count, alt_count) = self.get_counts_for_alt(alt_base);

            #[allow(
                clippy::cast_possible_truncation,
                reason = "alt_idx is from enumerate() on smallvec, fits in u8"
            )]
            let Some(alt_index) = NonZeroU8::new((alt_idx + 1) as u8) else {
                continue;
            };

            match EstimatedGenotype::calculate_for_alt(ref_count, alt_count, alt_index, error_model)
            {
                Ok(gt) => alt_genotypes.push(gt),
                Err(error) => {
                    trace!(%error, alt_base=%alt_base, "Failed to calculate genotype for alt");
                }
            }
        }

        if alt_genotypes.is_empty() {
            return None;
        }

        // Determine final genotype based on calculated genotypes
        if alt_genotypes.len() == 1 {
            // Single alt: return its genotype (0/1 or 1/1)
            alt_genotypes.into_iter().next()
        } else {
            // Multiple alts: determine best genotype
            // Sort by likelihood descending
            alt_genotypes.sort_by(|a, b| {
                b.likelihood.partial_cmp(&a.likelihood).unwrap_or(std::cmp::Ordering::Equal)
            });

            let top1 = alt_genotypes.first()?;
            let Some(top2) = alt_genotypes.get(1) else {
                return alt_genotypes.into_iter().next();
            };

            // Extract alt indices from genotypes
            let (GenotypeTag::RefHet(alt1_idx) | GenotypeTag::HomAlt(alt1_idx)) = top1.genotype
            else {
                return alt_genotypes.into_iter().next();
            };
            let (GenotypeTag::RefHet(alt2_idx) | GenotypeTag::HomAlt(alt2_idx)) = top2.genotype
            else {
                return alt_genotypes.into_iter().next();
            };

            // Check if we should call compound heterozygous (1/2)
            // Special case: if both are HomAlt, always call compound het
            let should_call_compound_het = match (top1.genotype, top2.genotype) {
                (GenotypeTag::HomAlt(_), GenotypeTag::HomAlt(_)) => {
                    // Both alts look like hom alt (no ref reads) -> must be compound het
                    true
                }
                (GenotypeTag::RefHet(_), GenotypeTag::RefHet(_)) => {
                    // Both are het with ref reads -> only call compound het if very similar
                    // Use a higher threshold (0.8) to be more conservative
                    let likelihood_ratio = *top2.likelihood / *top1.likelihood;
                    likelihood_ratio > 0.8
                }
                _ => {
                    // Mixed case (one HomAlt, one RefHet) -> use moderate threshold
                    let likelihood_ratio = *top2.likelihood / *top1.likelihood;
                    // Handle NaN (both likelihoods are 0) as similar
                    likelihood_ratio.is_nan() || likelihood_ratio > 0.5
                }
            };

            if should_call_compound_het {
                // Both alts have reasonable support, call compound het
                // Use combined probability for likelihood and confidence
                let combined_likelihood =
                    Probability::new((*top1.likelihood + *top2.likelihood) / 2.0)
                        .unwrap_or(top1.likelihood);
                let combined_confidence =
                    Probability::new((*top1.confidence + *top2.confidence) / 2.0)
                        .unwrap_or(top1.confidence);

                Some(EstimatedGenotype {
                    genotype: GenotypeTag::alt_het(alt1_idx, alt2_idx),
                    likelihood: combined_likelihood,
                    confidence: combined_confidence,
                })
            } else {
                // Top alt has much better support, use its genotype
                alt_genotypes.into_iter().next()
            }
        }
    }

    /// Get reference and alt read counts for a specific alt base.
    ///
    /// For C→T and G→A variants, uses strand-specific counting to avoid
    /// confounding with methylation. For other variants, uses both strands.
    fn get_counts_for_alt(&self, alt_base: Base) -> (usize, usize) {
        let ref_base = self.ref_base();
        let reads = || self.pileup.reads.iter();

        match (ref_base, alt_base) {
            (C, T) => (
                reads().filter(|r| r.base == C && r.strand == OB).count(),
                reads().filter(|r| r.base == T && r.strand == OB).count(),
            ),
            (G, A) => (
                reads().filter(|r| r.base == G && r.strand == OT).count(),
                reads().filter(|r| r.base == A && r.strand == OT).count(),
            ),
            _ => (
                reads().filter(|r| r.base == ref_base).count(),
                reads().filter(|r| r.base == alt_base).count(),
            ),
        }
    }
}

/// Represents a diploid genotype with support for multiple alternative alleles.
///
/// This enum provides a more type-safe and readable representation than raw VCF
/// genotype indices, while still supporting the full range of diploid genotypes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[must_use]
pub enum GenotypeTag {
    /// Homozygous reference (0/0)
    HomRef,
    /// Heterozygous with one reference and one alt allele (0/n where n > 0)
    RefHet(NonZeroU8),
    /// Compound heterozygous with two different alt alleles (m/n where m < n and m,n > 0)
    ///
    /// The alleles are always stored in sorted order (smaller first) to ensure uniqueness.
    AltHet(NonZeroU8, NonZeroU8),
    /// Homozygous alternative (n/n where n > 0)
    HomAlt(NonZeroU8),
}

impl GenotypeTag {
    /// Creates a homozygous reference genotype (0/0).
    ///
    /// This is equivalent to the legacy `CC` genotype.
    pub const fn hom_ref() -> Self {
        Self::HomRef
    }

    /// Creates a heterozygous genotype with one reference and one alt allele (0/n).
    ///
    /// When `alt_allele` is 1, this is equivalent to the legacy `CT` genotype.
    pub const fn ref_het(alt_allele: NonZeroU8) -> Self {
        Self::RefHet(alt_allele)
    }

    /// Creates a homozygous alternate genotype (n/n).
    pub const fn hom_alt(alt_allele: NonZeroU8) -> Self {
        Self::HomAlt(alt_allele)
    }

    /// Creates a compound heterozygous genotype with two different alt alleles (m/n).
    ///
    /// The alleles are automatically sorted to ensure canonical representation.
    /// If both alleles are the same, use [`GenotypeTag::hom_alt()`] instead.
    #[inline]
    pub fn alt_het(allele1: NonZeroU8, allele2: NonZeroU8) -> Self {
        assert_ne!(
            allele1, allele2,
            "Genotype with two different alleles expect, but got `{allele1}` twice."
        );
        // Store in sorted order for canonical representation
        if allele1 < allele2 {
            Self::AltHet(allele1, allele2)
        } else {
            Self::AltHet(allele2, allele1)
        }
    }

    /// Is this genotype heterozygous (any combination of different alleles)?
    pub const fn is_heterozygous(&self) -> bool {
        matches!(self, Self::RefHet(_) | Self::AltHet(_, _))
    }

    /// Is this genotype homozygous (both alleles the same)?
    pub const fn is_homozygous(&self) -> bool {
        matches!(self, Self::HomRef | Self::HomAlt(_))
    }

    pub const CC: Self = Self::HomRef;
    pub const CT: Self = Self::RefHet(NonZeroU8::new(1).unwrap());
    pub const TT: Self = Self::HomAlt(NonZeroU8::new(1).unwrap());
}

impl From<GenotypeTag> for [GenotypeAllele; 2] {
    fn from(value: GenotypeTag) -> Self {
        match value {
            GenotypeTag::HomRef => [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)],
            GenotypeTag::RefHet(n) => {
                [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(NonZeroI32::from(n).get())]
            }
            GenotypeTag::AltHet(m, n) => [
                GenotypeAllele::Unphased(NonZeroI32::from(m).get()),
                GenotypeAllele::Unphased(NonZeroI32::from(n).get()),
            ],
            GenotypeTag::HomAlt(n) => [
                GenotypeAllele::Unphased(NonZeroI32::from(n).get()),
                GenotypeAllele::Unphased(NonZeroI32::from(n).get()),
            ],
        }
    }
}

impl From<GenotypeTag> for Genotype {
    fn from(value: GenotypeTag) -> Self {
        let genotype: [GenotypeAllele; 2] = value.into();
        Genotype(genotype.into())
    }
}

impl TryFrom<&[GenotypeAllele]> for GenotypeTag {
    type Error = color_eyre::Report;

    fn try_from(alleles: &[GenotypeAllele]) -> Result<Self> {
        use color_eyre::eyre::bail;

        match alleles {
            // Homozygous reference: 0/0
            [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]
            | [GenotypeAllele::Phased(0), GenotypeAllele::Phased(0)] => Ok(Self::HomRef),

            // Heterozygous with reference: 0/n or n/0
            [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(n)]
            | [GenotypeAllele::Unphased(n), GenotypeAllele::Unphased(0)]
            | [GenotypeAllele::Phased(0), GenotypeAllele::Phased(n)]
            | [GenotypeAllele::Phased(n), GenotypeAllele::Phased(0)]
                if *n > 0 && *n <= 255 =>
            {
                #[allow(clippy::cast_possible_truncation, reason = "safe: n is in 1..=255")]
                let n = *n as u8;
                Ok(Self::RefHet(NonZeroU8::new(n).wrap_err("expected non-null value")?))
            }

            // Both alleles are alt (non-zero)
            [GenotypeAllele::Unphased(m), GenotypeAllele::Unphased(n)]
            | [GenotypeAllele::Phased(m), GenotypeAllele::Phased(n)]
                if *m > 0 && *m <= 255 && *n > 0 && *n <= 255 =>
            {
                #[allow(clippy::cast_possible_truncation, reason = "n,m in 1..=255")]
                let (m, n) = (*m as u8, *n as u8);
                let m = NonZeroU8::new(m).wrap_err("expected non-null value")?;
                let n = NonZeroU8::new(n).wrap_err("expected non-null value")?;

                if m == n {
                    // Homozygous alt: n/n
                    Ok(Self::HomAlt(m))
                } else {
                    // Compound heterozygous: m/n (sorted)
                    Ok(Self::alt_het(m, n))
                }
            }

            _ => bail!("Unsupported or invalid genotype format: {alleles:?}"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
#[must_use]
pub struct EstimatedGenotype {
    pub genotype: GenotypeTag,
    pub likelihood: Probability,
    pub confidence: Probability,
}

impl EstimatedGenotype {
    /// Calculate genotype for a specific alt allele using binomial model.
    ///
    /// Compares the probability of observing the read counts under different
    /// genotype hypotheses (hom ref, het, hom alt) and returns the most likely.
    #[instrument(level = "trace", name = "calculate_genotype")]
    pub fn calculate_for_alt(
        ref_count: usize,
        alt_count: usize,
        alt_index: NonZeroU8,
        error_model: ErrorModel,
    ) -> Result<Self> {
        let error_rate = error_model.error_rate();
        ensure!(
            ref_count > 0 || alt_count > 0,
            "No ref or alt read counts, cannot compute likelihood"
        );
        ensure!(*error_rate > f64::MIN, "Error rate too small, cannot calculate likelihood");

        // This is a simple estimate of genotype, based on the following consideration:
        // A site is either het or hom, where hom could be CC or TT.
        // If alt_count > ref_count, the latter is more likely, otherwise the former.

        // First, I calculate the likelihood to observe this many alt_reads
        // under the assumption that ref and alt are equally likely, ie this is a het position.

        // TODO: This assumes a simple diploid sample with no purity issues. For
        // cancer samples, we could make this a setting to allow for different cancer fraction?

        let mut binom = Binomial::new(ref_count + alt_count, 0.5); // 0.5 because a het position
        let p_het = binom.mass(alt_count);
        #[allow(
            clippy::cast_possible_truncation,
            reason = "safe: counts are usize, division result fits in usize"
        )]
        let p_het_max = binom.mass(((alt_count + ref_count).f() / 2.0).round() as usize);

        // Then, I calculate the probability that this many or more alt_count/ref_count reads
        // are observed by error, assuming independence of reads and errors.
        binom = Binomial::new(ref_count + alt_count, *error_rate);

        if ref_count >= alt_count {
            let p_hom = binom.mass(alt_count) + (1.0 - binom.distribution(alt_count.f()));

            if p_het < p_hom {
                trace!("Assuming 0/0: ({ref_count} vs {alt_count}) -> ({p_het:.5} < {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::hom_ref(),
                    likelihood: Probability::new(p_hom)?.inverted(),
                    confidence: Probability::new(((p_hom - p_het) / p_hom).clamp(0.0, 1.0))?
                        .inverted(),
                })
            } else {
                trace!(
                    "Assuming 0/{alt_index}: ({ref_count} vs {alt_count}) -> ({p_het:.5} >= {p_hom:.5})"
                );
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::ref_het(alt_index),
                    likelihood: Probability::new(p_het / p_het_max)?.inverted(),
                    confidence: Probability::new(((p_het - p_hom) / p_het).clamp(0.0, 1.0))?
                        .inverted(),
                })
            }
        } else {
            let p_hom = binom.mass(ref_count) + (1.0 - binom.distribution(ref_count.f()));
            if p_het < p_hom {
                trace!(
                    "Assuming {alt_index}/{alt_index}: ({ref_count} vs {alt_count}) -> ({p_het:.5} < {p_hom:.5})"
                );
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::hom_alt(alt_index),
                    likelihood: Probability::new(p_hom)?.inverted(),
                    confidence: Probability::new(((p_hom - p_het) / p_hom).clamp(0.0, 1.0))?
                        .inverted(),
                })
            } else {
                trace!(
                    "Assuming 0/{alt_index}: ({ref_count} vs {alt_count}) -> ({p_het:.5} >= {p_hom:.5})"
                );
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::ref_het(alt_index),
                    likelihood: Probability::new(p_het / p_het_max)?.inverted(),
                    confidence: Probability::new(((p_het - p_hom) / p_het).clamp(0.0, 1.0))?
                        .inverted(),
                })
            }
        }
    }
}
