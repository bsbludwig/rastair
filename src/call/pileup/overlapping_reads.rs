use crate::call::{pileup::SimpleRead, process::PileupMappingParams};
use rastair_types::SmallVec;
use rustc_hash::FxHashMap;
use std::{
    collections::hash_map::Entry,
    hash::{Hash, Hasher},
};

pub(crate) enum NameCollector {
    Skip,
    Collect(NameBuffer),
}

impl NameCollector {
    pub(crate) fn new(params: &PileupMappingParams) -> Self {
        if params.keep_overlapping_reads { Self::Skip } else { Self::Collect(NameBuffer::new()) }
    }
}

/// Reusable `HashMap` for dedup across pileup positions.
///
/// Maps each qname to the index of the first read seen with that name.
/// `prepare()` clears and pre-reserves capacity; the allocation persists
/// across positions so there is no per-position overhead.
pub(crate) struct NameBuffer(FxHashMap<NameKey, usize>);

impl NameBuffer {
    fn new() -> Self {
        Self(FxHashMap::default())
    }

    pub(super) fn prepare(&mut self, max_reads: usize) {
        self.0.clear();
        self.0.reserve(max_reads);
    }

    /// Record a read name. Returns the first-seen index if this is a duplicate.
    pub(super) fn see(&mut self, name: &[u8], this_idx: usize) -> Option<usize> {
        let key = NameKey { ptr: name.as_ptr(), len: name.len() };
        match self.0.entry(key) {
            Entry::Vacant(e) => {
                e.insert(this_idx);
                None
            }
            Entry::Occupied(e) => Some(*e.get()),
        }
    }
}

/// Raw qname pointer borrowed from htslib for one `from_hts` call.
///
/// Hashes by the last 4 bytes of the name for fast bucket distribution —
/// typical qnames encode tile/x/y coordinates there, so the suffix is highly
/// discriminating. Equality uses the full byte slice, so the Hash/Eq contract
/// holds: equal names always share a suffix.
///
/// # Safety
/// Only constructed from pointers that remain valid for the duration of the
/// enclosing `from_hts` call (htslib keeps pileup data alive until the next
/// `bam_plp_auto` step).
#[derive(Clone, Copy)]
struct NameKey {
    ptr: *const u8,
    len: usize,
}

impl NameKey {
    fn as_bytes(self) -> &'static [u8] {
        // Safety: caller must ensure that `ptr` is valid for `len` bytes and
        // that the data is not mutated
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
}

impl Hash for NameKey {
    /// Custom hashing that uses the last 8 byte of the name
    ///
    /// Read names often have common prefixes, so we focus on the suffix.
    fn hash<H: Hasher>(&self, state: &mut H) {
        const MAX_SUFFIX: usize = 8; // FxHashMap's default hash is 64-bit, so 8 bytes is enough to fill a hash value

        let bytes = self.as_bytes();
        if bytes.len() >= MAX_SUFFIX {
            state.write_u64(u64::from_ne_bytes(
                bytes[bytes.len() - MAX_SUFFIX..].try_into().unwrap(),
            ));
        } else {
            let n = bytes.len().min(MAX_SUFFIX);
            let mut suffix = [0u8; MAX_SUFFIX];
            suffix[MAX_SUFFIX - n..].copy_from_slice(&bytes[bytes.len() - n..]);
            state.write_u64(u64::from_ne_bytes(suffix));
        }
    }
}

impl PartialEq for NameKey {
    fn eq(&self, other: &Self) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl Eq for NameKey {}

/// Decide which of a duplicate pair to remove.
pub(super) fn resolve_pair(
    reads: &[SimpleRead],
    this_idx: usize,
    other_idx: usize,
    to_remove: &mut SmallVec<usize, 16>,
) {
    let this_read = &reads[this_idx];
    let other_read = &reads[other_idx];
    if this_read.base == other_read.base || this_read.second {
        to_remove.push(this_idx);
    } else {
        to_remove.push(other_idx);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{Base, default};

    fn sr(base: Base, second: bool) -> SimpleRead {
        SimpleRead { base, second, ..default() }
    }

    /// Simulate the inline dedup that `from_hts` performs during alignment collection.
    fn run_dedup(reads: Vec<(&[u8], SimpleRead)>) -> Vec<SimpleRead> {
        let mut buf = NameBuffer::new();
        buf.prepare(reads.len());
        let mut raw_reads: Vec<SimpleRead> = Vec::new();
        let mut to_remove: SmallVec<usize, 16> = SmallVec::new();
        for (name, read) in reads {
            let this_idx = raw_reads.len();
            raw_reads.push(read);
            if let Some(other_idx) = buf.see(name, this_idx) {
                resolve_pair(&raw_reads, this_idx, other_idx, &mut to_remove);
            }
        }
        to_remove.sort_unstable();
        for &idx in to_remove.iter().rev() {
            raw_reads.swap_remove(idx);
        }
        raw_reads
    }

    #[test]
    fn no_duplicates() {
        assert_eq!(
            run_dedup(vec![
                (b"read1".as_ref(), sr(Base::A, false)),
                (b"read2".as_ref(), sr(Base::C, false)),
                (b"read3".as_ref(), sr(Base::G, false)),
            ])
            .len(),
            3
        );
    }

    #[test]
    fn same_base_keeps_first() {
        let reads = run_dedup(vec![
            (b"read1".as_ref(), sr(Base::A, false)),
            (b"read1".as_ref(), sr(Base::A, true)),
        ]);
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].base, Base::A);
        assert!(!reads[0].second);
    }

    #[test]
    fn different_base_removes_second_in_pair() {
        let reads = run_dedup(vec![
            (b"read1".as_ref(), sr(Base::A, false)),
            (b"read1".as_ref(), sr(Base::C, true)),
        ]);
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].base, Base::A);
    }

    #[test]
    fn different_base_first_is_second() {
        let reads = run_dedup(vec![
            (b"read1".as_ref(), sr(Base::A, true)),
            (b"read1".as_ref(), sr(Base::C, false)),
        ]);
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].base, Base::C);
    }

    #[test]
    fn multiple_pairs() {
        assert_eq!(
            run_dedup(vec![
                (b"read1".as_ref(), sr(Base::A, false)),
                (b"read2".as_ref(), sr(Base::C, false)),
                (b"read1".as_ref(), sr(Base::A, true)),
                (b"read3".as_ref(), sr(Base::G, false)),
                (b"read2".as_ref(), sr(Base::C, true)),
            ])
            .len(),
            3
        );
    }

    #[test]
    fn same_suffix_different_name_kept() {
        // Both end in "_AAA" (same 4-byte suffix) but different full names
        assert_eq!(
            run_dedup(vec![
                (b"prefix_AAA".as_ref(), sr(Base::A, false)),
                (b"other__AAA".as_ref(), sr(Base::C, false)),
            ])
            .len(),
            2
        );
    }

    #[test]
    fn three_reads_same_name() {
        assert_eq!(
            run_dedup(vec![
                (b"read1".as_ref(), sr(Base::A, false)),
                (b"read1".as_ref(), sr(Base::A, true)),
                (b"read1".as_ref(), sr(Base::A, true)),
            ])
            .len(),
            1
        );
    }

    #[test]
    fn empty() {
        assert_eq!(run_dedup(vec![]).len(), 0);
    }

    #[test]
    fn realistic_qnames_not_deduped() {
        // /1 and /2 are distinct names — not duplicates
        assert_eq!(
            run_dedup(vec![
                (b"inst:1:fc:1:1234:5678:9012/1".as_ref(), sr(Base::A, false)),
                (b"inst:1:fc:1:1234:5678:9012/2".as_ref(), sr(Base::A, true)),
                (b"inst:1:fc:1:1234:9999:8888/1".as_ref(), sr(Base::C, false)),
            ])
            .len(),
            3
        );
    }
}
