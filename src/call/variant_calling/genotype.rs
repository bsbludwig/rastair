//! SNV genotypes: the SNV engine's dosage for a called alt, the binomial of a
//! CpG's reads for a column with none (`estimate_genotype`).
//!
//! - **A called alt takes the engine's genotype.** The forests decide which
//!   alts are present, so the pairs a genotype is chosen from carry at least
//!   one of the column's called alts and none of its others; the cluster's
//!   other events stay at their margins ([`Verdict::dosage_given`]). The
//!   likeliest such pair gives each called alt its copies and GQ as `GC`; one
//!   carrying both called alts is the compound `1/2`, at the lower GQ. The
//!   cascade scores every alt the forests can pass
//!   (`from_seqair::cascade_events`), so only a called alt at an `N`
//!   reference has no dosage: `./.` with no `GC`.
//! - **Without forests** (`--no-ml`) the engine's margins decide presence too
//!   ([`Verdict::dosage`]), and a called alt it finds absent is `0/0` at the
//!   verdict's `quality`.
//! - **A column with no called alt is `0/0`**, its `GC` the binomial of the
//!   CpG's no-SNP against SNP reads (`no_snp`/`snp`) under the error model: no
//!   event was scored there.
//! - **Two variants can form a CpG together** (`TA>CG`). `FormsDenovo` reads
//!   the reference context, so before ML neither alt forms one; once both are
//!   called, `process::mark_variant_pair_cpgs` marks them
//!   `ThisBecomesC`/`ThisBecomesG` and `cpg_origin` accepts that as the
//!   partner, for methylation calling. It assumes the two sit on one
//!   chromosome.
//!
//! Methylation calling reads the genotype (`GenotypeAdjustment::HomAlt`), so a
//! genotype flip moves that site's beta too. The `pileups!` fixtures in
//! `vcf_tests/genotyping.rs` are scored by the pair-HMM on their read lines
//! (`tests::utils::score_snv_lines`).
//!
//! [`Verdict::dosage`]: crate::call::pair_hmm::Verdict::dosage
//! [`Verdict::dosage_given`]: crate::call::pair_hmm::Verdict::dosage_given
use crate::{
    call::{
        pair_hmm::{Allele, Event, EventIdx, HapPair, Margin},
        variant_calling::ErrorModel,
    },
    metrics::{AltCall, PileupMetrics},
    utils::IntoF64 as _,
};
use color_eyre::eyre::{ContextCompat, Result};
use probability::prelude::{Binomial, Discrete as _, Distribution as _};
use rust_htslib::bcf::record::GenotypeAllele as HtslibGenotypeAllele;
use seqair_types::{Base, Probability, SmallVec};
use std::num::{NonZeroI32, NonZeroU8};

impl PileupMetrics {
    /// The genotype of the alts marked [`AltCall::RealVariant`], from the SNV
    /// engine's verdict; `0/0` from the reads' binomial when there are none.
    /// With `forests` they called the alts present; without, the engine does.
    /// `None` for a called alt the engine did not score, and for a column with
    /// no called alt and no reads to weigh.
    pub fn estimate_genotype(
        &self,
        error_model: ErrorModel,
        forests: bool,
    ) -> Option<EstimatedGenotype> {
        let called: SmallVec<(NonZeroU8, Base), 2> = self
            .alts
            .iter()
            .zip(1u8..)
            .filter(|(alt, _)| alt.call == AltCall::RealVariant)
            .filter_map(|(alt, index)| Some((NonZeroU8::new(index)?, alt.base)))
            .collect();
        if called.is_empty() {
            return self.hom_ref_genotype(error_model);
        }
        let scores = self.snv_scores.as_deref()?;
        let verdict = &scores.verdict;
        // The forests decide which of this column's alts are real, so a pair
        // carries at least one of them and no other alt here; which it carries
        // is the reads' call, as is every other event of the cluster.
        let here = |event: EventIdx| match scores.event(event) {
            Some(Event { anchor, allele: Allele::Base(base) })
                if forests && *anchor == self.pos =>
            {
                Some(called.iter().any(|(_, called)| called == base))
            }
            _ => None,
        };
        let admits = |pair: HapPair| {
            let mut any_called = false;
            for event in (0..scores.cluster().events().len()).filter_map(EventIdx::new) {
                let carried = pair.carries(event);
                match here(event) {
                    Some(true) => any_called |= carried,
                    Some(false) if carried => return false,
                    Some(false) => {}
                    None => {
                        if verdict.margin(event).is_some_and(Margin::is_present) != carried {
                            return false;
                        }
                    }
                }
            }
            any_called || !forests
        };
        let dosages: SmallVec<_, 2> = called
            .iter()
            .filter_map(|&(index, base)| {
                Some((index, verdict.dosage_given(scores.base_event_at(self.pos, base)?, admits)?))
            })
            .collect();
        let (genotype, gq) = match dosages.as_ref() {
            [] if forests => return None,
            [] => (GenotypeTag::hom_ref(), verdict.quality),
            [(index, dosage)] if dosage.copies >= 2 => (GenotypeTag::hom_alt(*index), dosage.gq),
            [(index, dosage)] => (GenotypeTag::ref_het(*index), dosage.gq),
            [(first, a), (second, b), ..] => {
                (GenotypeTag::alt_het(*first, *second), if *a.gq < *b.gq { a.gq } else { b.gq })
            }
        };
        let confidence = Probability::new(10f64.powf(-*gq / 10.0)).ok()?;
        Some(EstimatedGenotype { genotype, confidence })
    }

    /// `0/0`, with `GC` from the binomial of the CpG's no-SNP
    /// against SNP reads: het against the error rate's hom.
    fn hom_ref_genotype(&self, error_model: ErrorModel) -> Option<EstimatedGenotype> {
        let counts = self.pos_metrics.extended.methylation_strand_info;
        let (ref_count, alt_count) = (counts.no_snp as usize, counts.snp as usize);
        let reads = ref_count + alt_count;
        if reads == 0 {
            return None;
        }
        let het = Binomial::new(reads, 0.5);
        let p_het = het.mass(alt_count);
        let minority = ref_count.min(alt_count);
        let errors = Binomial::new(reads, *error_model.error_rate());
        let p_hom = errors.mass(minority) + (1.0 - errors.distribution(minority.f()));
        let confidence =
            if p_het < p_hom { (p_hom - p_het) / p_hom } else { (p_het - p_hom) / p_het };
        Some(EstimatedGenotype {
            genotype: GenotypeTag::hom_ref(),
            confidence: Probability::new(confidence.clamp(0.0, 1.0)).ok()?.inverted(),
        })
    }
}

/// Represents a diploid genotype with support for multiple alternative alleles.
///
/// This enum provides a more type-safe and readable representation than raw VCF
/// genotype indices, while still supporting the full range of diploid genotypes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

/// A single genotype allele, mirroring [`rust_htslib::bcf::record::GenotypeAllele`]
/// with serde support. Used when reading genotypes back from a VCF (e.g. BED
/// conversion); writing goes through seqair's [`Genotype`](seqair::vcf::Genotype).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenotypeAllele {
    /// Unphased allele with index `i`.
    Unphased(i32),
    /// Phased allele with index `i`.
    Phased(i32),
    /// Unphased missing allele.
    UnphasedMissing,
    /// Phased missing allele.
    PhasedMissing,
}

impl From<HtslibGenotypeAllele> for GenotypeAllele {
    fn from(allele: HtslibGenotypeAllele) -> Self {
        match allele {
            HtslibGenotypeAllele::Unphased(i) => GenotypeAllele::Unphased(i),
            HtslibGenotypeAllele::Phased(i) => GenotypeAllele::Phased(i),
            HtslibGenotypeAllele::UnphasedMissing => GenotypeAllele::UnphasedMissing,
            HtslibGenotypeAllele::PhasedMissing => GenotypeAllele::PhasedMissing,
        }
    }
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

    /// Three reads of each allele, and an engine whose six reads all fit the
    /// alt haplotype.
    fn hom_by_engine(alt: Base) -> PileupMetrics {
        use crate::call::{
            pair_hmm::{
                Allele, Cluster, Event, HaplotypeScores,
                fixtures::{NO, rows},
            },
            phasing::fixtures::column,
        };
        use std::sync::Arc;

        let mut site = column(1500, Base::A, &[(alt, AltCall::RealVariant)], None, &[]);
        let snv = Event { anchor: site.pos, allele: Allele::Base(alt) };
        let cluster = Cluster::new([snv].into_iter().collect()).expect("a cluster");
        let scores =
            HaplotypeScores::new(cluster, rows(&[([-4.0, -1.0, NO], 6)]).into()).expect("scores");
        site.snv_scores = Some(Arc::new(scores));
        site
    }

    /// A called alt takes the engine's dosage, and is `./.` where the engine
    /// did not score it.
    #[test]
    fn the_engine_decides_the_dosage() {
        let mut site = hom_by_engine(Base::G);
        let genotype = site.estimate_genotype(ErrorModel::default(), true).expect("a genotype");
        assert_eq!(genotype.genotype, GenotypeTag::hom_alt(ALT_INDEX_1));
        site.snv_scores = None;
        assert_eq!(site.estimate_genotype(ErrorModel::default(), true), None);
    }

    /// The SNV shares its cluster with a het deletion, and every read covering
    /// it carries the deletion as well: no read says whether the other
    /// chromosome carries the SNV, so `(ref, del+SNV)` and `(SNV, del+SNV)`
    /// tie: the earlier in VCF genotype order, the het, is written, at a `GC` that
    /// admits the tie.
    #[test]
    fn a_linked_snv_the_reads_cannot_place_is_a_low_quality_het() {
        use crate::call::{
            pair_hmm::{Allele, Cluster, Event, HaplotypeScores, ReadScores, fixtures::NO},
            pileup::indels::IndelAllele,
        };
        use std::sync::Arc;

        let mut site = hom_by_engine(Base::G);
        let deletion = IndelAllele::Deletion(std::iter::once(Base::T).collect());
        let events = [
            Event { anchor: site.pos + 14, allele: deletion.into() },
            Event { anchor: site.pos, allele: Allele::Base(Base::G) },
        ];
        let cluster = Cluster::new(events.into_iter().collect()).expect("a cluster");
        // By haplotype: reference, deletion, SNV, both.
        let both = ReadScores::new(None, 60, [-9.0, -4.5, -9.0, 0.0, NO, NO, NO, NO]);
        let other = ReadScores::new(None, 60, [0.0, -4.5, 0.0, -4.5, NO, NO, NO, NO]);
        let reads: Vec<_> = std::iter::repeat_n(both, 6).chain(std::iter::once(other)).collect();
        site.snv_scores =
            Some(Arc::new(HaplotypeScores::new(cluster, reads.into()).expect("scores")));
        let engine = site.estimate_genotype(ErrorModel::default(), true).expect("a genotype");
        assert_eq!(engine.genotype, GenotypeTag::ref_het(ALT_INDEX_1));
        let gc = -10.0 * engine.confidence.log10();
        assert!(gc <= 5.0, "{gc}");
    }

    /// A column with no called alt is `0/0`, surer the fewer SNP reads it has.
    #[test]
    fn a_column_without_a_called_alt_is_hom_ref() {
        use crate::call::phasing::fixtures::column;

        let gc = |no_snp, snp| {
            let mut site = column(1500, Base::A, &[], None, &[]);
            site.pos_metrics.extended.methylation_strand_info.no_snp = no_snp;
            site.pos_metrics.extended.methylation_strand_info.snp = snp;
            let genotype = site.estimate_genotype(ErrorModel::default(), true).expect("a genotype");
            assert_eq!(genotype.genotype, GenotypeTag::HomRef);
            -10.0 * genotype.confidence.log10()
        };
        assert!(gc(20, 0) > gc(10, 0), "{} vs {}", gc(20, 0), gc(10, 0));
        assert_eq!(
            column(1500, Base::A, &[], None, &[]).estimate_genotype(ErrorModel::default(), true),
            None
        );
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

#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct EstimatedGenotype {
    pub genotype: GenotypeTag,
    pub confidence: Probability,
}
