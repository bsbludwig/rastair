//! SNV genotypes: the SNV engine's genotype for a called alt
//! (`call_engine_genotype`, in `set_alt_calls`), the binomial of the column's
//! reads for a column with none (`hom_ref_genotype`).
//!
//! - **A called alt takes the engine's genotype.** The record's likelihoods
//!   are over its called alts, the column's other alts absent and the
//!   cluster's other events held at their calls
//!   ([`Verdict::genotype_likelihoods`]); `GT`, `GQ` and `PL` are read off
//!   them as the **Genotypes** section of `call::pair_hmm` describes, the
//!   forests deciding which alts are present.
//! - **An alt is called when the forests pass it and the genotype holds a
//!   copy** (`call_engine_genotype`, in `set_alt_calls`): an
//!   alt the engine did not score, or one the likeliest genotype leaves out, is
//!   not called, as an indel allele the pair-HMM does not call is not.
//! - **A column with no called alt is `0/0`**, its `GQ` from the binomials of
//!   its reads under the error model (at a CpG the no-SNP against SNP reads,
//!   `no_snp`/`snp`): no event was scored there. No `PL`.
//! - **Two variants can form a CpG together** (`TA>CG`). The forests see
//!   only the reference context, where neither alt forms one; once both are
//!   called, `process::mark_called_cpgs` records each beside the other
//!   (`PileupMetrics::cpg_side_of`), for methylation calling, when the
//!   pair-HMM puts the two on one chromosome (`Verdict::has_haplotype`). A
//!   reference partner whose own called alt sits on every chromosome carrying
//!   the de-novo alt makes no CpG either.
//!
//! Methylation calling reads the genotype (`GenotypeAdjustment`), so a
//! genotype flip moves that site's beta too. The `pileups!` fixtures in
//! `vcf_tests/genotyping.rs` are scored by the pair-HMM on their read lines
//! (`tests::utils::score_snv_lines`).
//!
//! [`Verdict::genotype_likelihoods`]: crate::call::pair_hmm::Verdict::genotype_likelihoods
use crate::{
    call::{
        pair_hmm::{Allele, Event, EventIdx, GenotypeLikelihoods, HaplotypeScores},
        variant_calling::ErrorModel,
    },
    metrics::{AltVerdict, MethylationEvidenceStrandInfo, PileupMetrics, Rejection},
    utils::IntoF64 as _,
};
use color_eyre::eyre::{ContextCompat, Result};
use rust_htslib::bcf::record::GenotypeAllele;
use seqair_types::{Base, Phred, Probability, SmallVec};
use std::num::NonZeroU8;

impl PileupMetrics {
    /// The SNV engine's genotype of the called alts, and the one place a
    /// passed alt is uncalled for it: a called alt the likeliest genotype
    /// carries no copy of (one the engine did not score, at an `N` or with no
    /// read it could score, or the second of two passed alts when the
    /// genotype carries only the other) is rejected for it, and the
    /// genotype is then of the alts left. `None` when none is left or the
    /// engine gives none.
    pub(crate) fn call_engine_genotype(&mut self) -> Option<EstimatedGenotype> {
        let called = self.called_alts();
        if called.is_empty() {
            return None;
        }
        let genotype = self.engine_genotype(&called);
        let mut uncalled = false;
        for &(index, base) in &called {
            let carried = genotype.as_ref().is_some_and(|genotype| {
                genotype.genotype.alleles().contains(&GtAllele::Alt(index))
            });
            if !carried
                && let Some(alt) = self.find_alt_mut(base)
                && let AltVerdict::Called(ml) = alt.verdict
            {
                alt.verdict = AltVerdict::Rejected(Rejection::NoGenotypeCopy(ml));
                uncalled = true;
            }
        }
        if !uncalled {
            return genotype;
        }
        let called = self.called_alts();
        if called.is_empty() {
            return None;
        }
        self.engine_genotype(&called)
    }

    fn engine_genotype(&self, called: &[(NonZeroU8, Base)]) -> Option<EstimatedGenotype> {
        let likelihoods = self.engine_likelihoods(called)?;
        let genotype = likelihoods.likeliest_with_an_alt()?;
        Some(EstimatedGenotype {
            genotype,
            quality: likelihoods.quality(genotype)?,
            source: GenotypeSource::Engine(Box::new(likelihoods)),
        })
    }

    /// The record's genotype likelihoods over the called alts the engine
    /// scored. The cluster's other events are held where the margins call
    /// them, this column's uncalled alts absent: the forests decide what is
    /// present here, the reads only which genotype.
    fn engine_likelihoods(&self, called: &[(NonZeroU8, Base)]) -> Option<AltLikelihoods> {
        let scores = self.snv_scores.as_deref()?;
        let verdict = scores.verdict();
        let scored: SmallVec<(NonZeroU8, EventIdx), 2> = called
            .iter()
            .filter_map(|&(index, base)| Some((index, scores.base_event_at(self.pos, base)?)))
            .collect();
        let events: SmallVec<EventIdx, 2> = scored.iter().map(|&(_, event)| event).collect();
        let log10 = verdict.genotype_likelihoods(&events, |event| {
            self.alt_here(scores, called, event).is_none() && verdict.is_called(event)
        })?;
        Some(AltLikelihoods { alts: scored.iter().map(|&(index, _)| index).collect(), log10 })
    }

    /// The alts marked `called`, by their 1-based index.
    fn called_alts(&self) -> SmallVec<(NonZeroU8, Base), 2> {
        self.alts
            .iter()
            .zip(1u8..)
            .filter(|(alt, _)| alt.called())
            .filter_map(|(alt, index)| Some((NonZeroU8::new(index)?, alt.base)))
            .collect()
    }

    /// Whether `event` is one of this column's alts, and if so whether it is
    /// called; `None` for the cluster's other events.
    fn alt_here(
        &self,
        scores: &HaplotypeScores,
        called: &[(NonZeroU8, Base)],
        event: EventIdx,
    ) -> Option<bool> {
        match scores.event(event) {
            Some(Event { anchor, allele: Allele::Base(base) }) if *anchor == self.pos => {
                Some(called.iter().any(|(_, called)| called == base))
            }
            _ => None,
        }
    }

    /// `0/0`, its `GQ` the chance the column is het (a binomial at one half)
    /// or hom-alt (the reference reads errors) rather than hom-ref (the alt
    /// reads errors), under a flat prior over the three. At a CpG the reads
    /// are the CpG's no-SNP and SNP reads on the strand TAPS leaves
    /// informative; elsewhere nothing converts and every read weighs. A
    /// de-novo partner is in a CpG only once its partner's alt is called, as
    /// methylation calling has it.
    /// `None` for a column with no reads to weigh.
    pub fn hom_ref_genotype(&self, error_model: ErrorModel) -> Option<EstimatedGenotype> {
        let (ref_count, alt_count) = if self.in_cpg() {
            let counts = MethylationEvidenceStrandInfo::from_pileup(self);
            (counts.no_snp as usize, counts.snp as usize)
        } else {
            let reads = |base| {
                self.allele(base)
                    .map_or(0, |allele| allele.strand_count.ot + allele.strand_count.ob)
            };
            let alts: u32 = self.alts.iter().map(|alt| reads(alt.base)).sum();
            (reads(self.ref_base()) as usize, alts as usize)
        };
        let reads = ref_count + alt_count;
        if reads == 0 {
            return None;
        }
        // The three binomials share their coefficient, so the log10 of each
        // likelihood without it is enough, and stays finite at any depth.
        let (alts, refs) = (alt_count.f(), ref_count.f());
        let error = error_model.error_rate().clamp(f64::MIN_POSITIVE, 0.5);
        let (wrong_read, right_read) = (error.log10(), (1.0 - error).log10());
        let hom_ref = alts * wrong_read + refs * right_read;
        let het = (alts + refs) * 0.5f64.log10();
        let hom_alt = refs * wrong_read + alts * right_read;
        let others = 10f64.powf(het - hom_ref) + 10f64.powf(hom_alt - hom_ref);
        // Written so an `inf` of others is certainly wrong, not NaN.
        let wrong = 1.0 - 1.0 / (1.0 + others);
        Some(EstimatedGenotype {
            genotype: GenotypeTag::HomRef,
            quality: Phred::from(Probability::new(wrong).ok()?),
            source: GenotypeSource::Binomial,
        })
    }
}

/// An allele as a genotype numbers it: the reference, or the `n`-th alt (of
/// `PileupMetrics::alts` for an SNV record, of the record's alleles for an
/// indel one). The VCF writer renumbers an SNV record's onto its `ALT` order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum GtAllele {
    Ref,
    Alt(NonZeroU8),
}

impl GtAllele {
    /// The allele's number in `GT`.
    #[must_use]
    pub const fn number(self) -> u8 {
        match self {
            Self::Ref => 0,
            Self::Alt(n) => n.get(),
        }
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
    /// The genotype of two alleles, in either order.
    pub fn from_alleles(a: GtAllele, b: GtAllele) -> Self {
        match (a, b) {
            (GtAllele::Ref, GtAllele::Ref) => Self::HomRef,
            (GtAllele::Ref, GtAllele::Alt(n)) | (GtAllele::Alt(n), GtAllele::Ref) => {
                Self::RefHet(n)
            }
            (GtAllele::Alt(m), GtAllele::Alt(n)) if m == n => Self::HomAlt(m),
            (GtAllele::Alt(m), GtAllele::Alt(n)) => Self::alt_het(m, n),
        }
    }

    /// The genotype of two alt alleles (m/n), stored sorted; one allele twice
    /// is n/n, [`GenotypeTag::HomAlt`].
    #[inline]
    pub fn alt_het(allele1: NonZeroU8, allele2: NonZeroU8) -> Self {
        match allele1.cmp(&allele2) {
            std::cmp::Ordering::Less => Self::AltHet(allele1, allele2),
            std::cmp::Ordering::Greater => Self::AltHet(allele2, allele1),
            std::cmp::Ordering::Equal => Self::HomAlt(allele1),
        }
    }

    /// Is this genotype heterozygous (any combination of different alleles)?
    pub const fn is_heterozygous(&self) -> bool {
        matches!(self, Self::RefHet(_) | Self::AltHet(_, _))
    }

    /// Is this genotype homozygous (both alleles the same)?
    #[cfg(test)]
    pub const fn is_homozygous(&self) -> bool {
        matches!(self, Self::HomRef | Self::HomAlt(_))
    }

    /// Whether either chromosome carries an alt.
    #[must_use]
    pub const fn carries_an_alt(self) -> bool {
        !matches!(self, Self::HomRef)
    }

    /// The two alleles, lower first.
    #[must_use]
    pub const fn alleles(self) -> [GtAllele; 2] {
        match self {
            Self::HomRef => [GtAllele::Ref; 2],
            Self::RefHet(n) => [GtAllele::Ref, GtAllele::Alt(n)],
            Self::AltHet(m, n) => [GtAllele::Alt(m), GtAllele::Alt(n)],
            Self::HomAlt(n) => [GtAllele::Alt(n); 2],
        }
    }

    #[cfg(test)]
    pub const CT: Self = Self::RefHet(NonZeroU8::new(1).unwrap());
    #[cfg(test)]
    pub const TT: Self = Self::HomAlt(NonZeroU8::new(1).unwrap());
}

/// Reads a genotype back from a VCF, as `convert` does.
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

#[derive(Debug, Clone)]
#[must_use]
pub struct EstimatedGenotype {
    pub genotype: GenotypeTag,
    /// `GQ`: Phred of the chance the genotype is wrong.
    pub quality: Phred,
    pub source: GenotypeSource,
}

/// What a genotype was estimated from.
#[derive(Debug, Clone)]
pub enum GenotypeSource {
    /// The SNV engine, with the likelihoods (`PL`) behind the `GQ`.
    Engine(Box<AltLikelihoods>),
    /// The binomial of a `0/0` column's reads, which has no `PL`.
    Binomial,
}

impl GenotypeSource {
    #[must_use]
    pub fn likelihoods(&self) -> Option<&AltLikelihoods> {
        match self {
            Self::Engine(likelihoods) => Some(likelihoods),
            Self::Binomial => None,
        }
    }
}

/// A record's genotype likelihoods and the alts they are over, by their
/// index in the column's alts, as [`GenotypeTag`] numbers them.
#[derive(Debug, Clone, PartialEq)]
pub struct AltLikelihoods {
    pub alts: SmallVec<NonZeroU8, 2>,
    pub log10: GenotypeLikelihoods,
}

impl AltLikelihoods {
    /// The likeliest genotype carrying one of the alts
    /// ([`GenotypeLikelihoods::likeliest_with_an_alt`]), numbered as the
    /// column numbers them.
    fn likeliest_with_an_alt(&self) -> Option<GenotypeTag> {
        let column = |allele: GtAllele| match allele {
            GtAllele::Ref => Some(GtAllele::Ref),
            GtAllele::Alt(n) => {
                Some(GtAllele::Alt(*self.alts.get(usize::from(n.get().checked_sub(1)?))?))
            }
        };
        let [a, b] = self.log10.likeliest_with_an_alt()?.alleles();
        Some(GenotypeTag::from_alleles(column(a)?, column(b)?))
    }

    /// `allele`, numbered as the column numbers it, as the likelihoods do.
    fn numbered(&self, allele: GtAllele) -> Option<GtAllele> {
        let GtAllele::Alt(n) = allele else { return Some(GtAllele::Ref) };
        let position = self.alts.iter().position(|&alt| alt == n)? + 1;
        Some(GtAllele::Alt(NonZeroU8::new(u8::try_from(position).ok()?)?))
    }

    /// `GQ` of `genotype`; `None` when it names an alt the likelihoods are
    /// not over.
    fn quality(&self, genotype: GenotypeTag) -> Option<Phred> {
        let [a, b] = genotype.alleles();
        Some(self.log10.quality(self.numbered(a)?, self.numbered(b)?))
    }

    /// The likelihoods over `vcf_alts` in that order, each a column alt
    /// index; `None` unless they are exactly this record's alts.
    #[must_use]
    pub fn over(&self, vcf_alts: &[NonZeroU8]) -> Option<GenotypeLikelihoods> {
        if vcf_alts.len() != self.alts.len() {
            return None;
        }
        let order: Option<SmallVec<usize, 3>> = std::iter::once(Some(0))
            .chain(vcf_alts.iter().map(|&n| Some(self.alts.iter().position(|&alt| alt == n)? + 1)))
            .collect();
        self.log10.reordered(&order?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::test_helpers::pos0;
    use crate::utils::MappingQuality;
    use hegel::{TestCase, generators as gs};

    const ALT_INDEX_1: NonZeroU8 = NonZeroU8::new(1).unwrap();

    #[test]
    fn alt_het_of_one_allele_twice_is_hom_alt() {
        assert_eq!(
            GenotypeTag::alt_het(ALT_INDEX_1, ALT_INDEX_1),
            GenotypeTag::HomAlt(ALT_INDEX_1)
        );
    }

    /// Three reads of each allele, and an engine whose six reads all fit the
    /// alt haplotype.
    fn hom_by_engine(alt: Base) -> PileupMetrics {
        use crate::call::{
            pair_hmm::{
                Allele, Cluster, Event, HaplotypeScores,
                fixtures::{NO, rows},
            },
            tests::fixtures::column,
        };
        use std::sync::Arc;

        let mut site = column(1500, Base::A, &[(alt, true)], None, &[]);
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
        let genotype = site.call_engine_genotype().expect("a genotype");
        assert_eq!(genotype.genotype, GenotypeTag::HomAlt(ALT_INDEX_1));
        site.snv_scores = None;
        assert!(site.call_engine_genotype().is_none());
    }

    /// The SNV shares its cluster with a het deletion, and every read covering
    /// it carries the deletion as well: no read says whether the other
    /// chromosome carries the SNV, so `(ref, del+SNV)` and `(SNV, del+SNV)`
    /// tie, and the het and the homozygote are written at a `GQ` that admits
    /// the tie.
    #[test]
    fn a_linked_snv_the_reads_cannot_place_is_a_low_quality_call() {
        use crate::call::{
            pair_hmm::{Allele, Cluster, Event, HaplotypeScores, ReadScores, fixtures::NO},
            pileup::indels::IndelAllele,
        };
        use std::sync::Arc;

        let mut site = hom_by_engine(Base::G);
        let deletion = IndelAllele::Deletion(std::iter::once(Base::T).collect());
        let events = [
            Event { anchor: pos0(site.pos.as_u32() + 14), allele: deletion.into() },
            Event { anchor: site.pos, allele: Allele::Base(Base::G) },
        ];
        let cluster = Cluster::new(events.into_iter().collect()).expect("a cluster");
        // By haplotype: reference, deletion, SNV, both.
        let both = ReadScores::new(
            None,
            MappingQuality::from_byte(60),
            [-9.0, -4.5, -9.0, 0.0, NO, NO, NO, NO],
        );
        let other = ReadScores::new(
            None,
            MappingQuality::from_byte(60),
            [0.0, -4.5, 0.0, -4.5, NO, NO, NO, NO],
        );
        let reads: Vec<_> = std::iter::repeat_n(both, 6).chain(std::iter::once(other)).collect();
        site.snv_scores =
            Some(Arc::new(HaplotypeScores::new(cluster, reads.into()).expect("scores")));
        let engine = site.call_engine_genotype().expect("a genotype");
        assert_ne!(engine.genotype, GenotypeTag::HomRef);
        let gq = *engine.quality;
        assert!(gq <= 5.0, "{gq}");
        // The tie is in `PL` too: the het and the homozygote equally likely.
        let pl = engine.source.likelihoods().expect("likelihoods").log10.phred();
        assert_eq!(pl.get(1), pl.get(2), "{pl:?}");
    }

    /// A column with no called alt is `0/0`, surer the fewer SNP reads it has.
    #[test]
    fn a_column_without_a_called_alt_is_hom_ref() {
        use crate::call::tests::fixtures::column;

        // At a CpG `C`, the OB reads: `no_snp` C, `snp` T.
        let gq = |no_snp, snp| {
            let mut site = column(1500, Base::C, &[(Base::T, false)], None, &[]);
            site.context.after_1 = Some(Base::G);
            site.ref_metrics.strand_count.ob = no_snp;
            if let Some(t) = site.alts.first_mut() {
                t.metrics.strand_count.ob = snp;
            }
            let genotype = site.hom_ref_genotype(ErrorModel::default()).expect("a genotype");
            assert_eq!(genotype.genotype, GenotypeTag::HomRef);
            assert!(matches!(genotype.source, GenotypeSource::Binomial));
            *genotype.quality
        };
        assert!(gq(20, 0) > gq(10, 0), "{} vs {}", gq(20, 0), gq(10, 0));
        assert!(gq(10, 0) > 20.0, "{}", gq(10, 0));
        // Reads a het explains are no support for `0/0`, nor are SNP reads alone.
        assert!(gq(10, 10) < 1.0, "{}", gq(10, 10));
        assert!(gq(0, 10) < 1.0, "{}", gq(0, 10));
        // Depths where the binomials themselves leave `f64`.
        assert_eq!(gq(0, 2_000), 0.0);
        assert!(gq(2_000, 0) > 99.0, "{}", gq(2_000, 0));
        let mut uncovered = column(1500, Base::C, &[], None, &[]);
        uncovered.context.after_1 = Some(Base::G);
        assert!(uncovered.hom_ref_genotype(ErrorModel::default()).is_none());
    }

    /// Outside a CpG nothing converts, so every read weighs.
    #[test]
    fn a_non_cpg_column_weighs_every_read() {
        use crate::call::tests::fixtures::column;

        let gq = |reference: u32, alt: u32| {
            let mut site = column(1500, Base::A, &[(Base::T, false)], None, &[]);
            site.ref_metrics.strand_count.ot = reference;
            if let Some(t) = site.alts.first_mut() {
                t.metrics.strand_count.ob = alt;
            }
            *site.hom_ref_genotype(ErrorModel::default()).expect("a genotype").quality
        };
        assert!(gq(12, 0) > 20.0, "{}", gq(12, 0));
        assert!(gq(6, 6) < 1.0, "{}", gq(6, 6));
    }

    // --- GenotypeTag property tests ---

    #[hegel::composite]
    fn nonzero_u8s(tc: &TestCase) -> NonZeroU8 {
        NonZeroU8::new(tc.draw(gs::integers::<u8>().min_value(1))).expect("nonzero")
    }

    /// `GT` is the likeliest of the record's own genotypes that carry a
    /// called alt: its `PL` is the least of theirs, whatever a linked
    /// event's phase does to how many pairs make each genotype.
    #[hegel::test]
    fn the_genotype_is_the_likeliest_of_its_likelihoods(tc: TestCase) {
        use crate::call::{
            pair_hmm::{Allele, Cluster, Event, HaplotypeScores, ReadScores, fixtures::NO},
            pileup::indels::IndelAllele,
        };
        use std::sync::Arc;

        let reads = tc.draw(gs::integers::<usize>().min_value(1).max_value(11));
        let rows: Vec<[f32; 4]> = tc.draw(
            gs::vecs(gs::arrays(gs::floats::<f32>().min_value(-12.0).max_value_exclusive(0.0)))
                .min_size(reads)
                .max_size(reads),
        );

        let mut site = hom_by_engine(Base::G);
        let deletion = IndelAllele::Deletion(std::iter::once(Base::T).collect());
        let events = [
            Event { anchor: pos0(site.pos.as_u32() + 14), allele: deletion.into() },
            Event { anchor: site.pos, allele: Allele::Base(Base::G) },
        ];
        let cluster = Cluster::new(events.into_iter().collect()).expect("a cluster");
        let reads: Vec<_> = rows
            .into_iter()
            .map(|[a, b, c, d]| {
                ReadScores::new(None, MappingQuality::from_byte(60), [a, b, c, d, NO, NO, NO, NO])
            })
            .collect();
        site.snv_scores =
            Some(Arc::new(HaplotypeScores::new(cluster, reads.into()).expect("scores")));
        let Some(engine) = site.call_engine_genotype() else {
            return;
        };
        let pl = engine.source.likelihoods().expect("likelihoods").log10.phred();
        let [a, b] = engine.genotype.alleles();
        let written = pl.get(GenotypeLikelihoods::index(a, b)).copied();
        let least = pl.iter().skip(1).copied().min();
        assert_eq!(written, least, "{:?} {:?}", engine.genotype, pl);
    }

    /// A written `GenotypeTag` parses back to itself.
    #[hegel::test]
    fn genotype_tag_roundtrips_through_alleles(tc: TestCase) {
        let n = tc.draw(nonzero_u8s());
        for tag in [GenotypeTag::HomRef, GenotypeTag::RefHet(n), GenotypeTag::HomAlt(n)] {
            let alleles = tag.alleles().map(|a| GenotypeAllele::Unphased(a.number().into()));
            let roundtripped =
                GenotypeTag::try_from(alleles.as_slice()).expect("valid alleles should roundtrip");
            assert_eq!(roundtripped, tag, "roundtrip failed for {tag:?}");
        }
    }

    /// `AltHet` alleles are always stored in sorted (canonical) order.
    #[hegel::test]
    fn alt_het_is_canonically_sorted(tc: TestCase) {
        let a = tc.draw(gs::integers::<u8>().min_value(1));
        // Any other allele: skip over `a` instead of rejecting `b == a`.
        let b = tc.draw(gs::integers::<u8>().min_value(1).max_value(254));
        let b = if b >= a { b + 1 } else { b };
        let (a, b) = (NonZeroU8::new(a).expect("nonzero"), NonZeroU8::new(b).expect("nonzero"));
        let tag = GenotypeTag::alt_het(a, b);
        let GenotypeTag::AltHet(lo, hi) = tag else {
            panic!("expected AltHet, got {tag:?}");
        };
        assert!(lo <= hi, "alleles not sorted: {lo} > {hi}");
    }

    /// `is_heterozygous` and `is_homozygous` are always opposite.
    #[hegel::test]
    fn het_and_hom_are_exclusive(tc: TestCase) {
        let n = tc.draw(nonzero_u8s());
        for tag in [GenotypeTag::HomRef, GenotypeTag::RefHet(n), GenotypeTag::HomAlt(n)] {
            assert_ne!(
                tag.is_heterozygous(),
                tag.is_homozygous(),
                "het and hom should be mutually exclusive for {tag:?}"
            );
        }
    }
}
