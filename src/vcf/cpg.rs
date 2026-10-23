use crate::{
    metrics::{FormsDenovo, PileupMetrics},
    utils::{Base, Base::*},
};
use seqair_types::Strand;
use std::fmt;
use std::ops::Deref;

/// Is this a CpG site?
#[derive(Clone, Copy, Default, PartialEq, Eq)]
pub enum InCpG {
    /// No, this is not a CpG site.
    #[default]
    No,
    /// Yes, first position in CpG.
    C,
    /// Yes, second position in CpG.
    G,
}

impl InCpG {
    /// The strand on which TAPS converts this side: OT reads the `C`, OB the
    /// `G` (as the `C` of the other strand).
    pub(crate) fn strand(self) -> Option<Strand> {
        match self {
            InCpG::C => Some(Strand::OT),
            InCpG::G => Some(Strand::OB),
            InCpG::No => None,
        }
    }

    /// The base a methylated read shows (T for C-side, A for G-side).
    pub(crate) fn mod_base(self) -> Option<Base> {
        match self {
            InCpG::C => Some(T),
            InCpG::G => Some(A),
            InCpG::No => None,
        }
    }

    /// The base an unmethylated read shows (C for C-side, G for G-side).
    pub(crate) fn unmod_base(self) -> Option<Base> {
        match self {
            InCpG::C => Some(C),
            InCpG::G => Some(G),
            InCpG::No => None,
        }
    }

    /// The partner base the side needs (G after C, C before G).
    pub(crate) fn adjacent_base(self) -> Option<Base> {
        match self {
            InCpG::C => Some(G),
            InCpG::G => Some(C),
            InCpG::No => None,
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

    /// The reference CpG side here, if any.
    pub fn in_cpg(self) -> InCpG {
        if self.reference == C && self.after == Some(G) {
            InCpG::C
        } else if self.reference == G && self.before == Some(C) {
            InCpG::G
        } else {
            InCpG::No
        }
    }

    /// Would `alt` in place of the reference form a CpG with a reference
    /// neighbour?
    pub fn forms_cpg(self, alt: Base) -> FormsDenovo {
        if alt == self.reference {
            FormsDenovo::No
        } else if self.before == Some(C) && alt == G {
            FormsDenovo::ThisBecomesG
        } else if self.after == Some(G) && alt == C {
            FormsDenovo::ThisBecomesC
        } else {
            FormsDenovo::No
        }
    }

    /// Is `alt` what a methylated reference CpG base reads as?
    pub fn is_conversion(self, alt: Base) -> bool {
        self.in_cpg().mod_base() == Some(alt)
    }

    /// The strand on which reads of `alt` may be methylation rather than the
    /// allele: the converting strand of a reference CpG read as its
    /// methylated base, or of the CpG `alt` itself forms, whose own reads
    /// convert there.
    pub fn confounded_strand(self, alt: Base) -> Option<Strand> {
        if self.is_conversion(alt) {
            return self.in_cpg().strand();
        }
        match self.forms_cpg(alt) {
            FormsDenovo::ThisBecomesC => InCpG::C.strand(),
            FormsDenovo::ThisBecomesG => InCpG::G.strand(),
            FormsDenovo::No => None,
        }
    }

    /// The strand that reads `alt` as it is where the other may convert it;
    /// `None` where neither does.
    pub fn informative_strand(self, alt: Base) -> Option<Strand> {
        self.confounded_strand(alt).map(|strand| match strand {
            Strand::OT => Strand::OB,
            Strand::OB => Strand::OT,
            Strand::Unknown => Strand::Unknown,
        })
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for InCpG {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InCpG::No => write!(f, "NoCpg"),
            InCpG::C => write!(f, "CpG::C"),
            InCpG::G => write!(f, "CpG::G"),
        }
    }
}

impl Deref for InCpG {
    type Target = bool;

    fn deref(&self) -> &Self::Target {
        match self {
            InCpG::No => &false,
            _ => &true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::test_helpers::variant_pileup;
    use color_eyre::Result;

    fn ctx(reference: Base, before: Option<Base>, after: Option<Base>) -> InCpG {
        TapsContext::new(reference, before, after).in_cpg()
    }

    #[test]
    fn test_in_cpg() {
        assert_eq!(ctx(A, None, None), InCpG::No);

        assert_eq!(ctx(C, None, Some(G)), InCpG::C);
        assert_eq!(ctx(C, Some(A), Some(G)), InCpG::C);

        assert_eq!(ctx(G, Some(C), None), InCpG::G);
        assert_eq!(ctx(G, Some(C), Some(T)), InCpG::G);

        assert_eq!(ctx(T, Some(C), Some(G)), InCpG::No);
        assert_eq!(ctx(C, Some(T), Some(A)), InCpG::No);
    }

    #[test]
    fn test_cpg_detection() -> Result<()> {
        let (_, pileup) = variant_pileup("chr19", 6105711)?;
        assert_eq!(TapsContext::of(&pileup).in_cpg(), InCpG::C);

        let (_, pileup) = variant_pileup("chr19", 6105712)?;
        assert_eq!(TapsContext::of(&pileup).in_cpg(), InCpG::G);

        Ok(())
    }

    #[test]
    fn conversion_needs_a_cpg() {
        let cpg_c = TapsContext::new(C, Some(A), Some(G));
        assert!(cpg_c.is_conversion(T));
        assert_eq!(cpg_c.confounded_strand(T), Some(Strand::OT));
        assert_eq!(cpg_c.informative_strand(T), Some(Strand::OB));
        assert!(!cpg_c.is_conversion(A));

        let cpg_g = TapsContext::new(G, Some(C), Some(T));
        assert_eq!(cpg_g.informative_strand(A), Some(Strand::OT));

        let plain_c = TapsContext::new(C, Some(A), Some(T));
        assert!(!plain_c.is_conversion(T));
        assert_eq!(plain_c.informative_strand(T), None);
    }

    #[test]
    fn a_denovo_cpg_allele_converts_on_its_own_side() {
        let before_g = TapsContext::new(T, Some(A), Some(G));
        assert_eq!(before_g.forms_cpg(C), FormsDenovo::ThisBecomesC);
        assert_eq!(before_g.informative_strand(C), Some(Strand::OB));

        let after_c = TapsContext::new(A, Some(C), Some(T));
        assert_eq!(after_c.forms_cpg(G), FormsDenovo::ThisBecomesG);
        assert_eq!(after_c.informative_strand(G), Some(Strand::OT));

        assert_eq!(TapsContext::new(A, Some(T), Some(T)).informative_strand(C), None);
    }
}
