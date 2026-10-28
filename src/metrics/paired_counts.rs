use seqair_types::{Base, Strand};

/// Key for looking up or incrementing a paired count.
///
/// `current` is the base at the position being queried; `adj` is the adjacent
/// base (before or after, depending on which `PairedCounts` field is accessed).
/// Unknown bases are silently ignored on both `get` and `increment`.
pub struct ReadKey {
    pub strand: Strand,
    pub current: Base,
    pub adj: Base,
}

/// Counts of (`current_base`, `adjacent_base`) pairs by strand.
///
/// Array-backed for O(1) access. Indexed as `[strand][current_base][adj_base]`.
///
/// Counters are `u32`, like every other per-position read count:
/// `--max-coverage 0` loads every read, so a deep amplicon can pass `u16::MAX`
/// reads at one position.
#[derive(Debug, Clone, Default)]
pub struct PairedCounts {
    inner: [[[u32; 4]; 4]; 2],
}

impl PairedCounts {
    pub fn get(&self, key: ReadKey) -> u32 {
        key.slot().and_then(|slot| self.flat().get(slot).copied()).unwrap_or(0)
    }

    pub fn increment(&mut self, key: ReadKey) {
        if let Some(slot) = key.slot()
            && let Some(count) = self.flat_mut().get_mut(slot)
        {
            *count = count.saturating_add(1);
        }
    }

    /// The table as one contiguous run of 32 counters. Every access is by a
    /// computed slot: three nested index expressions would be three bounds
    /// checks and three dependent address computations for one increment that
    /// runs per read per column.
    fn flat(&self) -> &[u32] {
        self.inner.as_flattened().as_flattened()
    }

    fn flat_mut(&mut self) -> &mut [u32] {
        self.inner.as_flattened_mut().as_flattened_mut()
    }
}

impl ReadKey {
    /// This key's slot in the flattened `[strand][current][adj]` table, or
    /// `None` when either base is unknown — those observations are silently
    /// not counted, on `get` as on `increment`.
    #[inline]
    fn slot(&self) -> Option<usize> {
        let strand = match self.strand {
            Strand::OT => 0,
            Strand::OB => 1,
        };
        // `known_index` is a table lookup, so this is three loads and a
        // multiply-add rather than three branches on data-dependent bases.
        Some(strand * 16 + self.current.known_index()? * 4 + self.adj.known_index()?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(strand: Strand, current: Base, adj: Base) -> ReadKey {
        ReadKey { strand, current, adj }
    }

    /// Every (strand, current, adj) triple must have its own slot, and no
    /// increment may land in another's: the flattened index arithmetic keeps
    /// the nested table's layout.
    #[test]
    fn every_known_triple_has_its_own_slot() {
        for strand in [Strand::OT, Strand::OB] {
            for current in Base::KNOWN {
                for adj in Base::KNOWN {
                    let mut counts = PairedCounts::default();
                    counts.increment(key(strand, current, adj));
                    assert_eq!(counts.get(key(strand, current, adj)), 1);

                    let others: u32 = [Strand::OT, Strand::OB]
                        .into_iter()
                        .flat_map(|s| {
                            Base::KNOWN
                                .into_iter()
                                .flat_map(move |c| Base::KNOWN.into_iter().map(move |a| (s, c, a)))
                        })
                        .filter(|&(s, c, a)| (s, c, a) != (strand, current, adj))
                        .map(|(s, c, a)| counts.get(key(s, c, a)))
                        .sum();
                    assert_eq!(others, 0, "{strand:?}/{current}/{adj} leaked into another slot");
                }
            }
        }
    }

    #[test]
    fn unknown_components_are_not_counted() {
        let mut counts = PairedCounts::default();
        counts.increment(key(Strand::OT, Base::Unknown, Base::C));
        counts.increment(key(Strand::OT, Base::A, Base::Unknown));

        assert_eq!(counts.get(key(Strand::OT, Base::Unknown, Base::C)), 0);
        assert_eq!(counts.get(key(Strand::OT, Base::A, Base::Unknown)), 0);
        let total: u32 = [Strand::OT, Strand::OB]
            .into_iter()
            .flat_map(|s| {
                Base::KNOWN
                    .into_iter()
                    .flat_map(move |c| Base::KNOWN.into_iter().map(move |a| (s, c, a)))
            })
            .map(|(s, c, a)| counts.get(key(s, c, a)))
            .sum();
        assert_eq!(total, 0, "an unknown base must count nothing at all");
    }

    #[test]
    fn increments_accumulate() {
        let mut counts = PairedCounts::default();
        for _ in 0..7 {
            counts.increment(key(Strand::OB, Base::G, Base::T));
        }
        assert_eq!(counts.get(key(Strand::OB, Base::G, Base::T)), 7);
    }

    /// `--max-coverage 0` loads every read, so a counter must count a depth
    /// past `u16::MAX` exactly.
    #[test]
    fn a_counter_counts_past_u16() {
        let mut counts = PairedCounts::default();
        let depth = u32::from(u16::MAX) + 10;
        for _ in 0..depth {
            counts.increment(key(Strand::OT, Base::C, Base::G));
        }
        assert_eq!(counts.get(key(Strand::OT, Base::C, Base::G)), depth);
    }
}
