use crate::{
    metrics::PileupMetrics,
    utils::{Base, Base::*, ByStrand},
};
use enumset::EnumSet;
use seqair_types::Strand;

/// One side of a CpG, reference or de-novo.
#[derive(Debug, enumset::EnumSetType)]
pub enum CpgSide {
    C,
    G,
}

/// A CpG side's reads, on the strand that converts it and on the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ByConversion {
    pub(crate) converting: u32,
    pub(crate) plain: u32,
}

impl CpgSide {
    /// The side a `C` or a `G` takes in a CpG.
    pub(crate) fn of_base(base: Base) -> Option<CpgSide> {
        match base {
            C => Some(CpgSide::C),
            G => Some(CpgSide::G),
            _ => None,
        }
    }

    /// The sides a base takes in a CpG with its genomic neighbours, given the
    /// sides of the bases there ([`Self::of_base`]): a `G` after a `C`, a `C`
    /// before a `G`.
    ///
    /// The one CpG rule. The neighbours' bases are the reference's for a
    /// reference CpG or an alt making one with the reference
    /// ([`TapsContext`]); a neighbour's candidate alts while the column loop
    /// decides which columns to keep; its called alts once the calls are in
    /// (`process::mark_called_cpgs`).
    pub(crate) fn beside(before: EnumSet<CpgSide>, after: EnumSet<CpgSide>) -> EnumSet<CpgSide> {
        let mut sides = EnumSet::empty();
        if before.contains(CpgSide::C) {
            sides.insert(CpgSide::G);
        }
        if after.contains(CpgSide::G) {
            sides.insert(CpgSide::C);
        }
        sides
    }

    /// The side `base` takes in a CpG, if `sides` ([`Self::beside`]) has it.
    pub(crate) fn taken_by(base: Base, sides: EnumSet<CpgSide>) -> Option<CpgSide> {
        Self::of_base(base).filter(|&side| sides.contains(side))
    }

    /// The side whose methylation TAPS reads as `observed` over `reference`:
    /// `T` over `C`, `A` over `G`. Says nothing about the context; a CpG may
    /// be the reference's or one a variant makes.
    pub(crate) fn converting(reference: Base, observed: Base) -> Option<CpgSide> {
        [CpgSide::C, CpgSide::G]
            .into_iter()
            .find(|side| side.unmod_base() == reference && side.mod_base() == observed)
    }

    /// The strand on which TAPS converts this side: OT reads the `C`, OB the
    /// `G` (as the `C` of the other strand).
    pub(crate) fn strand(self) -> Strand {
        match self {
            CpgSide::C => Strand::OT,
            CpgSide::G => Strand::OB,
        }
    }

    /// The side TAPS converts on `strand`, the inverse of [`Self::strand`].
    pub(crate) fn converted_on(strand: Strand) -> CpgSide {
        match strand {
            Strand::OT => CpgSide::C,
            Strand::OB => CpgSide::G,
        }
    }

    /// Whether one read on the strand converting this side shows a CpG here
    /// methylated: its modified base is, its unmodified base is not, and any
    /// other base (the other side's, a variant's) says nothing. Per read, from
    /// its base alone, never from the position's β.
    pub(crate) fn read_methylated(self, observed: Base) -> Option<bool> {
        if observed == self.mod_base() {
            Some(true)
        } else if observed == self.unmod_base() {
            Some(false)
        } else {
            None
        }
    }

    /// `reads` split into those on the strand that converts this side, and
    /// those on the other, which read the base as it is.
    pub(crate) fn converting_and_plain(self, reads: ByStrand<u32>) -> ByConversion {
        let (converting, plain) = match self.strand() {
            Strand::OT => (reads.ot, reads.ob),
            Strand::OB => (reads.ob, reads.ot),
        };
        ByConversion { converting, plain }
    }

    /// The base a methylated read shows (T for C-side, A for G-side).
    pub(crate) fn mod_base(self) -> Base {
        match self {
            CpgSide::C => T,
            CpgSide::G => A,
        }
    }

    /// The base an unmethylated read shows (C for C-side, G for G-side).
    pub(crate) fn unmod_base(self) -> Base {
        match self {
            CpgSide::C => C,
            CpgSide::G => G,
        }
    }

    /// The partner base the side needs (G after C, C before G).
    pub(crate) fn adjacent_base(self) -> Base {
        match self {
            CpgSide::C => G,
            CpgSide::G => C,
        }
    }
}

/// One reference position as TAPS sees it: the reference base and its
/// neighbours, which decide whether a substitution there can be methylation
/// rather than a variant.
///
/// TAPS reads a methylated CpG `C` as `T`. The `C` side converts on OT; the
/// `G` side converts on OB, where it is the other strand's `C` and shows as
/// `A`. The opposite strand of each reads the base as it is, so it is the
/// informative one. Outside a CpG, reference or de-novo, nothing converts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TapsContext {
    reference: Base,
    before: Option<Base>,
    after: Option<Base>,
}

impl TapsContext {
    pub fn new(reference: Base, before: Option<Base>, after: Option<Base>) -> Self {
        Self { reference, before, after }
    }

    pub fn of(pileup: &PileupMetrics) -> Self {
        Self::new(pileup.reference_base, pileup.context.before_1, pileup.context.after_1)
    }

    /// The CpG sides a base here takes with the reference neighbours.
    pub(crate) fn sides(self) -> EnumSet<CpgSide> {
        let side_of = |base: Option<Base>| base.and_then(CpgSide::of_base).into_iter().collect();
        CpgSide::beside(side_of(self.before), side_of(self.after))
    }

    /// The reference CpG side here, if any.
    pub fn in_cpg(self) -> Option<CpgSide> {
        CpgSide::taken_by(self.reference, self.sides())
    }

    /// The side of a CpG `alt` in place of the reference would form with a
    /// reference neighbour.
    pub fn forms_cpg(self, alt: Base) -> Option<CpgSide> {
        if alt == self.reference { None } else { CpgSide::taken_by(alt, self.sides()) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::test_helpers::variant_pileup;
    use color_eyre::Result;

    fn ctx(reference: Base, before: Option<Base>, after: Option<Base>) -> Option<CpgSide> {
        TapsContext::new(reference, before, after).in_cpg()
    }

    #[test]
    fn test_in_cpg() {
        assert_eq!(ctx(A, None, None), None);

        assert_eq!(ctx(C, None, Some(G)), Some(CpgSide::C));
        assert_eq!(ctx(C, Some(A), Some(G)), Some(CpgSide::C));

        assert_eq!(ctx(G, Some(C), None), Some(CpgSide::G));
        assert_eq!(ctx(G, Some(C), Some(T)), Some(CpgSide::G));

        assert_eq!(ctx(T, Some(C), Some(G)), None);
        assert_eq!(ctx(C, Some(T), Some(A)), None);
    }

    #[test]
    fn test_cpg_detection() -> Result<()> {
        let (_, pileup) = variant_pileup("chr19", 6105711)?;
        assert_eq!(TapsContext::of(&pileup).in_cpg(), Some(CpgSide::C));

        let (_, pileup) = variant_pileup("chr19", 6105712)?;
        assert_eq!(TapsContext::of(&pileup).in_cpg(), Some(CpgSide::G));

        Ok(())
    }

    #[test]
    fn an_alt_forms_a_cpg_with_its_reference_neighbour() {
        assert_eq!(TapsContext::new(T, Some(A), Some(G)).forms_cpg(C), Some(CpgSide::C));
        assert_eq!(TapsContext::new(A, Some(C), Some(T)).forms_cpg(G), Some(CpgSide::G));
        assert_eq!(TapsContext::new(A, Some(T), Some(T)).forms_cpg(C), None);
        assert_eq!(TapsContext::new(G, Some(C), Some(T)).forms_cpg(G), None);
    }

    #[test]
    fn a_g_follows_a_c_and_a_c_precedes_a_g() {
        let (c, g, none) = (EnumSet::only(CpgSide::C), EnumSet::only(CpgSide::G), EnumSet::empty());
        assert_eq!(CpgSide::beside(c, none), g);
        assert_eq!(CpgSide::beside(none, g), c);
        assert_eq!(CpgSide::beside(g, c), none);
        assert_eq!(CpgSide::beside(c, g), EnumSet::all());
        assert_eq!(CpgSide::taken_by(G, CpgSide::beside(c, none)), Some(CpgSide::G));
        assert_eq!(CpgSide::taken_by(C, CpgSide::beside(c, none)), None);
        assert_eq!(CpgSide::taken_by(A, EnumSet::all()), None);
    }

    #[test]
    fn only_t_over_c_and_a_over_g_are_conversions() {
        assert_eq!(CpgSide::converting(C, T), Some(CpgSide::C));
        assert_eq!(CpgSide::converting(G, A), Some(CpgSide::G));
        for (reference, observed) in [(C, A), (G, T), (T, C), (A, G), (C, C), (A, T)] {
            assert_eq!(
                CpgSide::converting(reference, observed),
                None,
                "{reference:?}>{observed:?}"
            );
        }
        let reads = ByStrand { ot: 3, ob: 5 };
        assert_eq!(
            CpgSide::C.converting_and_plain(reads),
            ByConversion { converting: 3, plain: 5 }
        );
        assert_eq!(
            CpgSide::G.converting_and_plain(reads),
            ByConversion { converting: 5, plain: 3 }
        );
    }
}
