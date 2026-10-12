use seqair_types::Strand;

/// Helper struct to hold values for top and bottom strands
#[derive(Debug, Clone)]
pub struct ByStrand<T> {
    /// Value for the top strand
    pub ot: T,
    /// Value for the bottom strand
    pub ob: T,
}

impl<T: Copy> Copy for ByStrand<T> {}

impl<T: Default> Default for ByStrand<T> {
    fn default() -> Self {
        ByStrand { ot: T::default(), ob: T::default() }
    }
}

impl<T> ByStrand<T> {
    /// The value for a known strand; `None` for [`Strand::Unknown`].
    pub const fn get(&self, strand: Strand) -> Option<&T> {
        match strand {
            Strand::OT => Some(&self.ot),
            Strand::OB => Some(&self.ob),
            Strand::Unknown => None,
        }
    }

    pub const fn get_mut(&mut self, strand: Strand) -> Option<&mut T> {
        match strand {
            Strand::OT => Some(&mut self.ot),
            Strand::OB => Some(&mut self.ob),
            Strand::Unknown => None,
        }
    }
}

impl ByStrand<u32> {
    pub const fn total(&self) -> u32 {
        self.ot.saturating_add(self.ob)
    }
}

/// One of a phase set's two haplotypes.
///
/// Which is which is only meaningful inside one set: haplotype 1 is the
/// chromosome the set's anchor put its first allele on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Haplotype {
    First,
    Second,
}

/// A value per haplotype of one phase set, in the shape of [`ByStrand`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ByHaplotype<T> {
    pub first: T,
    pub second: T,
}

impl<T> ByHaplotype<T> {
    pub const fn new(first: T, second: T) -> Self {
        Self { first, second }
    }

    pub const fn get(&self, haplotype: Haplotype) -> &T {
        match haplotype {
            Haplotype::First => &self.first,
            Haplotype::Second => &self.second,
        }
    }

    pub const fn get_mut(&mut self, haplotype: Haplotype) -> &mut T {
        match haplotype {
            Haplotype::First => &mut self.first,
            Haplotype::Second => &mut self.second,
        }
    }
}
