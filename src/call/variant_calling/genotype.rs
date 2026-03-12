// Adapted from rastair1, (c) Benjamin Schuster-Boeckler
use crate::{
    call::variant_calling::ErrorModel,
    metrics::{AltCall, FormsDenovo, PileupMetrics},
    utils::{Base, Base::*, IntoF64 as _, Strand::*},
    vcf::InCpG,
};
use color_eyre::eyre::{ContextCompat, Result, ensure};
use probability::prelude::{Binomial, Discrete as _, Distribution as _};
use rastair_types::{Probability, SmallVec};
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
        let real_variant_alts: SmallVec<_, 2> = self
            .alts
            .iter()
            .enumerate()
            .filter(|(_, alt)| matches!(alt.call, AltCall::RealVariant))
            .collect();

        if real_variant_alts.is_empty() {
            // No real variants - calculate HomRef genotype using no_snp and snp counts
            // like rastair1 did.
            let counts = self.pos_metrics.extended.methylation_strand_info;
            let ref_count = counts.no_snp as usize;
            let alt_count = counts.snp as usize;

            if ref_count == 0 && alt_count == 0 {
                return None;
            }

            const ONE: NonZeroU8 = NonZeroU8::new(1).expect("1 > 0");

            return EstimatedGenotype::calculate_for_alt(ref_count, alt_count, ONE, error_model)
                .ok()
                .map(|mut gt| {
                    // Force genotype to HomRef since we have no real variants
                    gt.genotype = GenotypeTag::hom_ref();
                    gt
                });
        }

        // If ML threshold is set, filter by ML score
        let passing_alts: SmallVec<_, 2> = if let Some(threshold) = ml_threshold {
            real_variant_alts
                .into_iter()
                .filter(|(_, alt)| alt.filters.ml.is_some_and(|ml| ml >= threshold))
                .collect()
        } else {
            real_variant_alts
        };

        // Calculate genotype for each passing alt using binomial model
        // Note: We compare each alt independently against ref.
        // Alternative approaches:
        // - compare alt1 vs alt2 for compound het
        // - compare each alt against total depth minus its count
        let mut alt_genotypes: SmallVec<_, 2> = SmallVec::new();

        for (alt_idx, alt) in &passing_alts {
            let alt_base = alt.base;
            let forms_denovo = alt.metrics.denovo;
            let (ref_count, alt_count) = self.get_counts_for_alt(alt_base, forms_denovo);

            // Add methylation evidence to ref count for genotype calculation.
            // This is the same adjustment for all alts at this position,
            // because we're comparing each alt independently against the same
            // reference.
            let ref_count_with_methylation = ref_count;

            #[allow(
                clippy::cast_possible_truncation,
                reason = "alt_idx is from enumerate() on smallvec, fits in u8"
            )]
            let Some(alt_index) = NonZeroU8::new((alt_idx + 1) as u8) else {
                continue;
            };

            match EstimatedGenotype::calculate_for_alt(
                ref_count_with_methylation,
                alt_count,
                alt_index,
                error_model,
            ) {
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
    /// For CpG positions (reference C or G), uses strand-specific counting to avoid
    /// confounding with methylation. For non-CpG positions, uses both strands.
    ///
    /// De-novo CpGs (X→C or X→G that create new CpG contexts) are treated
    /// the same as regular CpG variants:
    /// - X→C (ThisBecomesC): Uses OB strand only (like C→T)
    /// - X→G (ThisBecomesG): Uses OT strand only (like G→A)
    ///
    /// This prevents methylation evidence from confounding genotype calls.
    fn get_counts_for_alt(&self, alt_base: Base, forms_denovo: FormsDenovo) -> (usize, usize) {
        let ref_base = self.ref_base();
        let reads = || self.pileup.reads.iter();
        let obs = || reads().filter(|r| r.strand == OB);
        let ots = || reads().filter(|r| r.strand == OT);
        let count_all = || {
            (
                reads().filter(|r| r.base == ref_base).count(),
                reads().filter(|r| r.base == alt_base).count(),
            )
        };
        let count_ob = || {
            (
                obs().filter(|r| r.base == ref_base).count(),
                obs().filter(|r| r.base == alt_base).count(),
            )
        };
        let count_ot = || {
            (
                ots().filter(|r| r.base == ref_base).count(),
                ots().filter(|r| r.base == alt_base).count(),
            )
        };

        // For de-novo CpGs, treat them like CpG SNPs to avoid methylation confounding:
        // - ThisBecomesC: new C can be methylated (C→T on OT strand), so use OB strand only
        // - ThisBecomesG: partner C can be methylated (shows as A on OB strand), so use OT strand only
        match forms_denovo {
            FormsDenovo::ThisBecomesC if ref_base == T => count_ob(),
            FormsDenovo::ThisBecomesC => {
                let ref_count = reads().filter(|r| r.base == ref_base).count();
                let alt_count = obs().filter(|r| r.base == C).count()
                    + ots().filter(|r| matches!(r.base, C | T)).count();
                (ref_count, alt_count)
            }
            FormsDenovo::ThisBecomesG if ref_base == A => count_ot(),
            FormsDenovo::ThisBecomesG => {
                let ref_count = reads().filter(|r| r.base == ref_base).count();
                let alt_count = ots().filter(|r| r.base == G).count()
                    + obs().filter(|r| matches!(r.base, G | A)).count();
                (ref_count, alt_count)
            }
            // For regular variants (not de-novo CpGs), use CpG-aware counting when needed.
            FormsDenovo::No => match self.pos_metrics.cpg {
                InCpG::C if alt_base == T => count_ob(),
                InCpG::C => {
                    let ref_count = obs().filter(|r| r.base == C).count()
                        + ots().filter(|r| matches!(r.base, C | T)).count();
                    let alt_count = reads().filter(|r| r.base == alt_base).count();
                    (ref_count, alt_count)
                }
                InCpG::G if alt_base == A => count_ot(),
                InCpG::G => {
                    let ref_count = ots().filter(|r| r.base == G).count()
                        + obs().filter(|r| matches!(r.base, G | A)).count();
                    let alt_count = reads().filter(|r| r.base == alt_base).count();
                    (ref_count, alt_count)
                }
                InCpG::No => count_all(),
            },
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

    /// Is this genotype homozygous reference (0/0)?
    pub const fn is_hom_ref(&self) -> bool {
        matches!(self, Self::HomRef)
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

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const ALT_INDEX_1: NonZeroU8 = NonZeroU8::new(1).unwrap();

    // Strategy: realistic read counts (1..=500 each, not both zero)
    fn read_counts() -> impl Strategy<Value = (usize, usize)> {
        (1usize..=500, 0usize..=500)
            .prop_filter("at least one read required", |&(ref_c, alt_c)| ref_c > 0 || alt_c > 0)
    }

    // Strategy: a valid error rate well above f64::MIN
    fn error_model() -> impl Strategy<Value = ErrorModel> {
        (1u32..=100u32).prop_map(|n| {
            ErrorModel::Custom(Probability::new(f64::from(n) / 10_000.0).expect("valid rate"))
        })
    }

    proptest! {
        /// Likelihood and confidence are always valid probabilities (in [0,1], not NaN).
        #[test]
        fn outputs_are_valid_probabilities(
            (ref_count, alt_count) in read_counts(),
            error in error_model(),
        ) {
            let result = EstimatedGenotype::calculate_for_alt(ref_count, alt_count, ALT_INDEX_1, error);
            let Ok(gt) = result else { return Ok(()); };

            prop_assert!(*gt.likelihood >= 0.0 && *gt.likelihood <= 1.0,
                "likelihood out of range: {}", *gt.likelihood);
            prop_assert!(*gt.confidence >= 0.0 && *gt.confidence <= 1.0,
                "confidence out of range: {}", *gt.confidence);
            prop_assert!(gt.likelihood.is_finite(), "likelihood is NaN/inf");
            prop_assert!(gt.confidence.is_finite(), "confidence is NaN/inf");
        }

        /// When all reads support the reference (alt_count == 0), genotype must be HomRef.
        #[test]
        fn all_ref_reads_gives_hom_ref(
            ref_count in 10usize..=500,
            error in error_model(),
        ) {
            let gt = EstimatedGenotype::calculate_for_alt(ref_count, 0, ALT_INDEX_1, error)
                .expect("should succeed with nonzero counts");
            prop_assert_eq!(gt.genotype, GenotypeTag::HomRef,
                "expected HomRef with ref={} alt=0", ref_count);
        }

        /// When all reads support the alt (ref_count == 0), genotype must be HomAlt.
        #[test]
        fn all_alt_reads_gives_hom_alt(
            alt_count in 10usize..=500,
            error in error_model(),
        ) {
            let gt = EstimatedGenotype::calculate_for_alt(0, alt_count, ALT_INDEX_1, error)
                .expect("should succeed with nonzero counts");
            prop_assert_eq!(gt.genotype, GenotypeTag::HomAlt(ALT_INDEX_1),
                "expected HomAlt with ref=0 alt={}", alt_count);
        }

        /// ref > alt never produces HomAlt; alt > ref never produces HomRef.
        ///
        /// The majority base should always "win" — the model must not call the
        /// minority allele as homozygous.
        #[test]
        fn majority_base_determines_hom_direction(
            ref_count in 1usize..=500,
            extra in 10usize..=200,
            error in error_model(),
        ) {
            let alt_count = ref_count + extra; // alt clearly dominates
            let gt_alt_wins = EstimatedGenotype::calculate_for_alt(ref_count, alt_count, ALT_INDEX_1, error)
                .expect("nonzero counts");
            prop_assert_ne!(gt_alt_wins.genotype, GenotypeTag::HomRef,
                "got HomRef even though alt ({}) > ref ({})", alt_count, ref_count);

            let gt_ref_wins = EstimatedGenotype::calculate_for_alt(alt_count, ref_count, ALT_INDEX_1, error)
                .expect("nonzero counts");
            prop_assert_ne!(gt_ref_wins.genotype, GenotypeTag::HomAlt(ALT_INDEX_1),
                "got HomAlt even though ref ({}) > alt ({})", alt_count, ref_count);
        }

        /// Confidence increases as the read counts become more extreme.
        ///
        /// A site with 490 ref / 10 alt should have higher confidence in HomRef
        /// than a site with 300 ref / 200 alt.
        #[test]
        fn more_extreme_counts_give_higher_confidence(
            base in 100usize..=400,
            small_imbalance in 10usize..=40,
            large_imbalance in 80usize..=200,
            error in error_model(),
        ) {
            prop_assume!(base > large_imbalance);
            let ref_large = base + large_imbalance;
            let alt_large = base - large_imbalance;
            let ref_small = base + small_imbalance;
            let alt_small = base - small_imbalance;

            let gt_extreme = EstimatedGenotype::calculate_for_alt(ref_large, alt_large, ALT_INDEX_1, error)
                .expect("nonzero counts");
            let gt_moderate = EstimatedGenotype::calculate_for_alt(ref_small, alt_small, ALT_INDEX_1, error)
                .expect("nonzero counts");

            // Only compare when both agree on the same genotype class to avoid
            // comparing confidence values across different decision boundaries.
            if gt_extreme.genotype == gt_moderate.genotype {
                prop_assert!(
                    *gt_extreme.confidence >= *gt_moderate.confidence,
                    "more extreme counts ({ref_large}/{alt_large}) should have confidence >= moderate ({ref_small}/{alt_small}): \
                     extreme={:.4} moderate={:.4}",
                    *gt_extreme.confidence,
                    *gt_moderate.confidence,
                );
            }
        }

        /// Zero counts for both ref and alt must return an error.
        #[test]
        fn zero_counts_is_error(error in error_model()) {
            let result = EstimatedGenotype::calculate_for_alt(0, 0, ALT_INDEX_1, error);
            prop_assert!(result.is_err(), "expected error for (0, 0) counts");
        }
    }

    // --- GenotypeTag property tests ---

    fn nonzero_u8_strategy() -> impl Strategy<Value = NonZeroU8> {
        (1u8..=255).prop_map(|n| NonZeroU8::new(n).expect("nonzero"))
    }

    proptest! {
        /// GenotypeTag → [GenotypeAllele; 2] → GenotypeTag must roundtrip.
        #[test]
        fn genotype_tag_roundtrips_through_alleles(
            n in nonzero_u8_strategy(),
        ) {
            for tag in [
                GenotypeTag::HomRef,
                GenotypeTag::RefHet(n),
                GenotypeTag::HomAlt(n),
            ] {
                let alleles: [GenotypeAllele; 2] = tag.into();
                let roundtripped = GenotypeTag::try_from(alleles.as_slice())
                    .expect("valid alleles should roundtrip");
                prop_assert_eq!(roundtripped, tag, "roundtrip failed for {:?}", tag);
            }
        }

        /// AltHet alleles are always stored in sorted (canonical) order.
        #[test]
        fn alt_het_is_canonically_sorted(
            a in nonzero_u8_strategy(),
            b in nonzero_u8_strategy(),
        ) {
            prop_assume!(a != b);
            let tag = GenotypeTag::alt_het(a, b);
            let GenotypeTag::AltHet(lo, hi) = tag else {
                prop_assert!(false, "expected AltHet");
                return Ok(());
            };
            prop_assert!(lo <= hi, "alleles not sorted: {lo} > {hi}");
        }

        /// is_heterozygous and is_homozygous are always opposite.
        #[test]
        fn het_and_hom_are_exclusive(n in nonzero_u8_strategy()) {
            for tag in [
                GenotypeTag::HomRef,
                GenotypeTag::RefHet(n),
                GenotypeTag::HomAlt(n),
            ] {
                prop_assert_ne!(
                    tag.is_heterozygous(), tag.is_homozygous(),
                    "het and hom should be mutually exclusive for {:?}", tag
                );
            }
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
