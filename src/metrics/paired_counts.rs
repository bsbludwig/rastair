use seqair_types::{Base, Strand};

/// Key for looking up or incrementing a paired count.
///
/// `current` is the base at the position being queried; `adj` is the adjacent
/// base (before or after, depending on which `PairedCounts` field is accessed).
/// Unknown bases or strands are silently ignored on both `get` and `increment`.
pub struct ReadKey {
    pub strand: Strand,
    pub current: Base,
    pub adj: Base,
}

/// Counts of (`current_base`, `adjacent_base`) pairs by strand.
///
/// Array-backed for O(1) access. Indexed as `[strand][current_base][adj_base]`.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct PairedCounts {
    inner: [[[u32; 4]; 4]; 2],
}

impl PairedCounts {
    pub fn get(&self, key: ReadKey) -> u32 {
        let (Some(si), Some(ci), Some(ai)) =
            (strand_idx(key.strand), base_idx(key.current), base_idx(key.adj))
        else {
            return 0;
        };
        self.inner[si][ci][ai]
    }

    pub fn increment(&mut self, key: ReadKey) {
        let (Some(si), Some(ci), Some(ai)) =
            (strand_idx(key.strand), base_idx(key.current), base_idx(key.adj))
        else {
            return;
        };
        self.inner[si][ci][ai] += 1;
    }
}

fn base_idx(base: Base) -> Option<usize> {
    match base {
        Base::A => Some(0),
        Base::C => Some(1),
        Base::G => Some(2),
        Base::T => Some(3),
        Base::Unknown => None,
    }
}

fn strand_idx(strand: Strand) -> Option<usize> {
    match strand {
        Strand::OT => Some(0),
        Strand::OB => Some(1),
        Strand::Unknown => None,
    }
}
