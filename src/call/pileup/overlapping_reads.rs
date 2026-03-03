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
        if params.keep_overlapping_reads {
            Self::Skip
        } else {
            Self::Collect(NameBuffer::new(params.linear_dedup_threshold))
        }
    }
}

/// Reusable buffer for read-name deduplication across pileup positions.
///
/// Below `threshold` reads, uses a linear scan through parallel suffix/name
/// arrays — cache-friendly and SIMD-amenable at low depth. Above the
/// threshold, falls back to an `FxHashMap`.
///
/// `prepare()` picks the active mode and clears/reserves the appropriate inner
/// buffer; `see()` delegates without allocation. Both inner buffers persist
/// across positions so only the per-position clear+reserve happens on the hot path.
pub(crate) struct NameBuffer {
    threshold: usize,
    hash: HashInner,
    linear: LinearInner,
    active: DedupeMode,
}

#[derive(Clone, Copy)]
enum DedupeMode {
    Hash,
    Linear,
}

impl NameBuffer {
    pub(crate) fn new(threshold: usize) -> Self {
        Self {
            threshold,
            hash: HashInner(FxHashMap::default()),
            linear: LinearInner::new(),
            active: DedupeMode::Hash,
        }
    }

    pub(super) fn prepare(&mut self, max_reads: usize) {
        if max_reads <= self.threshold {
            self.active = DedupeMode::Linear;
            self.linear.prepare(max_reads);
        } else {
            self.active = DedupeMode::Hash;
            self.hash.prepare(max_reads);
        }
    }

    pub(super) fn see(&mut self, name: &[u8], this_idx: usize) -> Option<usize> {
        match self.active {
            DedupeMode::Hash => self.hash.see(name, this_idx),
            DedupeMode::Linear => self.linear.see(name),
        }
    }
}

/// Hashmap-based inner buffer (used above `threshold`).
///
/// Maps each qname to the index of the first read seen with that name.
/// `prepare()` clears and pre-reserves capacity; the allocation persists
/// across positions so there is no per-position overhead.
struct HashInner(FxHashMap<NameKey, usize>);

impl HashInner {
    fn prepare(&mut self, max_reads: usize) {
        self.0.clear();
        self.0.reserve(max_reads);
    }

    fn see(&mut self, name: &[u8], this_idx: usize) -> Option<usize> {
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

/// Linear-scan inner buffer (used at or below `threshold`).
///
/// Stores parallel arrays: `suffixes[i]` is the last-4-byte suffix of the
/// i-th read seen, `keys[i]` is its full name pointer. The vec index equals
/// the read's position in `raw_reads`, so no separate index field is needed.
///
/// On each `see()`, scans `suffixes` for a matching u32 — the tight loop over
/// a contiguous array of u32s is easily auto-vectorised by the compiler — then
/// verifies the full name only on a suffix hit.
///
/// Both vecs are always cleared and reserved to `capacity` on `prepare()` so
/// their allocations survive across positions and grow monotonically to the
/// high-water mark.
struct LinearInner {
    suffixes: Vec<u32>,
    keys: Vec<NameKey>,
}

impl LinearInner {
    fn new() -> Self {
        Self { suffixes: Vec::new(), keys: Vec::new() }
    }

    fn prepare(&mut self, capacity: usize) {
        self.suffixes.clear();
        self.suffixes.reserve(capacity);
        self.keys.clear();
        self.keys.reserve(capacity);
    }

    /// Record a read name and return the first-seen index if this is a duplicate.
    ///
    /// Always appends to both inner vecs so that the vec index stays in sync
    /// with the caller's `raw_reads` index.
    fn see(&mut self, name: &[u8]) -> Option<usize> {
        let suffix = name_suffix_u32(name);
        let key = NameKey { ptr: name.as_ptr(), len: name.len() };

        // Zip suffixes and keys so a suffix collision on a different full name
        // doesn't prevent finding the true duplicate at a later position.
        let duplicate = self
            .suffixes
            .iter()
            .zip(self.keys.iter())
            .position(|(&s, k)| s == suffix && k.as_bytes() == name);

        self.suffixes.push(suffix);
        self.keys.push(key);

        duplicate
    }
}

/// Extract the last 4 bytes of `name` as a native-endian u32.
///
/// Read names carry tile/x/y coordinates in their suffix; 4 bytes is enough
/// to discriminate most pairs while fitting 8 entries per 32-byte cache line.
fn name_suffix_u32(name: &[u8]) -> u32 {
    const SUFFIX: usize = 4;
    let mut buf = [0u8; SUFFIX];
    let n = name.len().min(SUFFIX);
    buf[SUFFIX - n..].copy_from_slice(&name[name.len() - n..]);
    u32::from_ne_bytes(buf)
}

/// Raw qname pointer borrowed from htslib for one `from_hts` call.
///
/// Hashes by the last 8 bytes of the name for fast bucket distribution —
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
    /// Custom hashing that uses the last 8 bytes of the name.
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
    use proptest::prelude::*;

    fn sr(base: Base, second: bool) -> SimpleRead {
        SimpleRead { base, second, ..default() }
    }

    /// Simulate the inline dedup that `from_hts` performs during alignment collection.
    fn run_dedup_with(reads: Vec<(&[u8], SimpleRead)>, threshold: usize) -> Vec<SimpleRead> {
        let mut buf = NameBuffer::new(threshold);
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

    fn run_dedup(reads: Vec<(&[u8], SimpleRead)>) -> Vec<SimpleRead> {
        run_dedup_with(reads, 0) // threshold=0 → always hash
    }

    fn run_dedup_linear(reads: Vec<(&[u8], SimpleRead)>) -> Vec<SimpleRead> {
        run_dedup_with(reads, usize::MAX) // threshold=MAX → always linear
    }

    #[test]
    fn no_duplicates() {
        let reads = vec![
            (b"read1".as_ref(), sr(Base::A, false)),
            (b"read2".as_ref(), sr(Base::C, false)),
            (b"read3".as_ref(), sr(Base::G, false)),
        ];
        assert_eq!(run_dedup(reads.clone()).len(), 3);
        assert_eq!(run_dedup_linear(reads).len(), 3);
    }

    #[test]
    fn same_base_keeps_first() {
        let reads = vec![
            (b"read1".as_ref(), sr(Base::A, false)),
            (b"read1".as_ref(), sr(Base::A, true)),
        ];
        for result in [run_dedup(reads.clone()), run_dedup_linear(reads)] {
            assert_eq!(result.len(), 1);
            assert_eq!(result[0].base, Base::A);
            assert!(!result[0].second);
        }
    }

    #[test]
    fn different_base_removes_second_in_pair() {
        let reads = vec![
            (b"read1".as_ref(), sr(Base::A, false)),
            (b"read1".as_ref(), sr(Base::C, true)),
        ];
        for result in [run_dedup(reads.clone()), run_dedup_linear(reads)] {
            assert_eq!(result.len(), 1);
            assert_eq!(result[0].base, Base::A);
        }
    }

    #[test]
    fn different_base_first_is_second() {
        let reads = vec![
            (b"read1".as_ref(), sr(Base::A, true)),
            (b"read1".as_ref(), sr(Base::C, false)),
        ];
        for result in [run_dedup(reads.clone()), run_dedup_linear(reads)] {
            assert_eq!(result.len(), 1);
            assert_eq!(result[0].base, Base::C);
        }
    }

    #[test]
    fn multiple_pairs() {
        let reads = vec![
            (b"read1".as_ref(), sr(Base::A, false)),
            (b"read2".as_ref(), sr(Base::C, false)),
            (b"read1".as_ref(), sr(Base::A, true)),
            (b"read3".as_ref(), sr(Base::G, false)),
            (b"read2".as_ref(), sr(Base::C, true)),
        ];
        assert_eq!(run_dedup(reads.clone()).len(), 3);
        assert_eq!(run_dedup_linear(reads).len(), 3);
    }

    #[test]
    fn same_suffix_different_name_kept() {
        // Both end in "_AAA" (same 4-byte suffix) but different full names.
        // This exercises the "suffix collision but different name" path in LinearInner.
        let reads = vec![
            (b"prefix_AAA".as_ref(), sr(Base::A, false)),
            (b"other__AAA".as_ref(), sr(Base::C, false)),
        ];
        assert_eq!(run_dedup(reads.clone()).len(), 2);
        assert_eq!(run_dedup_linear(reads).len(), 2);
    }

    #[test]
    fn suffix_collision_finds_true_duplicate() {
        // read0 and read1 share a 4-byte suffix but are different names.
        // read2 is a true duplicate of read1. The linear scan must not stop
        // at the false suffix hit on read0 and must find the match at read1.
        let reads = vec![
            (b"prefix_AAA".as_ref(), sr(Base::A, false)),
            (b"other__AAA".as_ref(), sr(Base::C, false)),
            (b"other__AAA".as_ref(), sr(Base::C, true)),
        ];
        assert_eq!(run_dedup(reads.clone()).len(), 2);
        assert_eq!(run_dedup_linear(reads).len(), 2);
    }

    #[test]
    fn three_reads_same_name() {
        let reads = vec![
            (b"read1".as_ref(), sr(Base::A, false)),
            (b"read1".as_ref(), sr(Base::A, true)),
            (b"read1".as_ref(), sr(Base::A, true)),
        ];
        assert_eq!(run_dedup(reads.clone()).len(), 1);
        assert_eq!(run_dedup_linear(reads).len(), 1);
    }

    #[test]
    fn empty() {
        assert_eq!(run_dedup(vec![]).len(), 0);
        assert_eq!(run_dedup_linear(vec![]).len(), 0);
    }

    #[test]
    fn realistic_qnames_not_deduped() {
        // /1 and /2 are distinct names — not duplicates
        let reads = vec![
            (b"inst:1:fc:1:1234:5678:9012/1".as_ref(), sr(Base::A, false)),
            (b"inst:1:fc:1:1234:5678:9012/2".as_ref(), sr(Base::A, true)),
            (b"inst:1:fc:1:1234:9999:8888/1".as_ref(), sr(Base::C, false)),
        ];
        assert_eq!(run_dedup(reads.clone()).len(), 3);
        assert_eq!(run_dedup_linear(reads).len(), 3);
    }

    /// Strategy: generate a small pool of distinct names, then produce a
    /// sequence of reads that samples from the pool (creating realistic
    /// duplicates). The pool size is 1–6; sequence length is 0–30.
    fn reads_strategy() -> impl Strategy<Value = Vec<(Vec<u8>, Base, bool)>> {
        prop::collection::vec(prop::collection::vec(any::<u8>(), 1usize..=30), 1usize..=6)
            .prop_flat_map(|pool| {
                let n = pool.len();
                prop::collection::vec(
                    (
                        0..n,
                        prop_oneof![Just(Base::A), Just(Base::C), Just(Base::G), Just(Base::T)],
                        any::<bool>(),
                    ),
                    0usize..=30,
                )
                .prop_map(move |entries| {
                    entries
                        .into_iter()
                        .map(|(i, base, second)| (pool[i].clone(), base, second))
                        .collect()
                })
            })
    }

    proptest! {
        /// Both dedup modes must produce the same multiset of surviving reads
        /// for any input sequence. Compared as sorted `(base as u8, second)`
        /// pairs since `SimpleRead` has no `PartialEq`.
        #[test]
        fn hash_and_linear_agree(entries in reads_strategy()) {
            let reads: Vec<(&[u8], SimpleRead)> = entries
                .iter()
                .map(|(name, base, second)| (name.as_slice(), sr(*base, *second)))
                .collect();

            let mut hash_out = run_dedup_with(reads.clone(), 0);
            let mut linear_out = run_dedup_with(reads, usize::MAX);

            // Sort by (base discriminant, second) for order-independent comparison.
            // swap_remove doesn't preserve insertion order, but both modes remove
            // the same logical reads so the multisets must be identical.
            let key = |r: &SimpleRead| (r.base as u8, r.second);
            hash_out.sort_unstable_by_key(key);
            linear_out.sort_unstable_by_key(key);

            prop_assert_eq!(hash_out.len(), linear_out.len());
            for (h, l) in hash_out.iter().zip(linear_out.iter()) {
                prop_assert_eq!((h.base, h.second), (l.base, l.second));
            }
        }
    }
}
