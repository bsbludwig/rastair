// Adapted from rastair1, (c) Benjamin Schuster-Boeckler
use crate::{
    call::{pileup::Pileup, variant_calling::ErrorModel},
    utils::{Base::*, IntoF64 as _, Strand},
};
use color_eyre::eyre::{ContextCompat, Result, ensure};
use probability::prelude::{Binomial, Discrete as _, Distribution as _};
use rastair_types::Probability;
use rastair_vcf::standard_fields::{Genotype, GenotypeAllele};
use std::num::{NonZeroI32, NonZeroU8};
use tracing::{instrument, trace};

impl Pileup {
    pub fn estimate_genotype(&self, error_model: ErrorModel) -> Option<EstimatedGenotype> {
        let (nosnp, snp) = if self.reference_base == C {
            (
                self.reads.iter().filter(|b| b.base == C && b.strand == Strand::OB).count(),
                self.reads.iter().filter(|b| b.base == T && b.strand == Strand::OB).count(),
            )
        } else if self.reference_base == G {
            (
                self.reads.iter().filter(|b| b.base == G && b.strand == Strand::OT).count(),
                self.reads.iter().filter(|b| b.base == A && b.strand == Strand::OT).count(),
            )
        } else {
            return None;
        };
        match EstimatedGenotype::calculate(nosnp, snp, error_model) {
            Ok(gt) => Some(gt),
            Err(error) => {
                trace!(%error, "Failed to calculate genotype");
                None
            }
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
        binom = Binomial::new(ref_count + alt_count, error_rate);

        if ref_count >= alt_count {
            let p_hom = binom.mass(alt_count) + (1.0 - binom.distribution(alt_count.f()));

            if p_het < p_hom {
                trace!("Assuming CC: ({ref_count} vs {alt_count}) -> ({p_het:.5} < {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::hom_ref(),
                    likelihood: Probability::new(1.0 - p_hom)?,
                    confidence: Probability::new(1.0 - (p_hom - p_het) / p_hom)?,
                })
            } else {
                trace!("Assuming CT: ({ref_count} vs {alt_count}) -> ({p_het:.5} >= {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::CT,
                    likelihood: Probability::new(1.0 - p_het / p_het_max)?,
                    confidence: Probability::new(1.0 - (p_het - p_hom) / p_het)?,
                })
            }
        } else {
            let p_hom = binom.mass(ref_count) + (1.0 - binom.distribution(ref_count.f()));
            if p_het < p_hom {
                trace!("Assuming TT: ({ref_count} vs {alt_count}) -> ({p_het:.5} < {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::TT,
                    likelihood: Probability::new(1.0 - p_hom)?,
                    confidence: Probability::new(1.0 - (p_hom - p_het) / p_hom)?,
                })
            } else {
                trace!("Assuming TC: ({ref_count} vs {alt_count}) -> ({p_het:.5} >= {p_hom:.5})");
                Ok(EstimatedGenotype {
                    genotype: GenotypeTag::CT,
                    likelihood: Probability::new(1.0 - p_het / p_het_max)?,
                    confidence: Probability::new((1.0 - p_het - p_hom) / p_het)?,
                })
            }
        }
    }
}
