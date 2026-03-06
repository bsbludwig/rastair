use seqair_types::{Base, Strand};

/// Counts of (current_base, adjacent_base) pairs by strand.
///
/// Array-backed for O(1) access. Indexed as `[strand][current_base][adj_base]`.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct PairedCounts {
    inner: [[[u32; 4]; 4]; 2],
}

impl PairedCounts {
    pub fn get(&self, strand: Strand, current: Base, adj: Base) -> u32 {
        let (Some(si), Some(ci), Some(ai)) = (strand_idx(strand), base_idx(current), base_idx(adj))
        else {
            return 0;
        };
        self.inner[si][ci][ai]
    }

    pub fn increment(&mut self, strand: Strand, current: Base, adj: Base) {
        let (Some(si), Some(ci), Some(ai)) = (strand_idx(strand), base_idx(current), base_idx(adj))
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
