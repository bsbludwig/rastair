//! A small, efficient byte vector type that uses inline storage for small vectors
//
// Ideas:
// - view into a larger array? maybe actually just not use this at all and have
//   a view into a known big array per genome segment?
#![deny(missing_docs)]

use std::{
    cmp::Ordering,
    fmt,
    hash::{Hash, Hasher},
    ops::Deref,
    sync::Arc,
};

/// A structure that stores byte vectors efficiently, using inline storage for small vectors
/// and shared heap storage for larger ones. This is immutable but cheaply clonable.
///
/// # Examples
///
/// Creating a small vector (uses inline storage):
///
/// ```
/// # use smallarcvec::SmallByteVec;
/// let small = SmallByteVec::from([1, 2, 3]);
/// assert_eq!(small.len(), 3);
/// assert_eq!(&small[..], &[1, 2, 3]);
/// assert!(small.is_inline());
/// ```
///
/// Creating a larger vector (uses heap storage with Arc):
///
/// ```
/// # use smallarcvec::SmallByteVec;
/// let large_data: Vec<u8> = (0..30).collect();
/// let large = SmallByteVec::from(&large_data);
/// assert_eq!(large.len(), 30);
/// assert!(!large.is_inline());
/// ```
///
/// Using as a byte slice through Deref:
///
/// ```
/// # use smallarcvec::SmallByteVec;
/// let bytes = SmallByteVec::from([1, 2, 3, 4]);
///
/// // Can use slice methods directly
/// let first_two = &bytes[0..2];
/// assert_eq!(first_two, &[1, 2]);
///
/// // Can iterate over bytes
/// let sum: u32 = bytes.iter().map(|&b| b as u32).sum();
/// assert_eq!(sum, 10); // 1 + 2 + 3 + 4
/// ```
pub struct SmallByteVec {
    inner: SmallByteVecInner,
}

const INLINE_CAPACITY: usize = 22;

enum SmallByteVecInner {
    // Up to 22 bytes inline
    Inline {
        // Could get one more byte if we used an enum with 23 variants so the rest are niches
        // like <https://docs.rs/smol_str/0.3.2/src/smol_str/lib.rs.html#439>
        len: u8,
        data: [u8; INLINE_CAPACITY],
    },
    // For larger arrays, use Arc
    Heap {
        data: Arc<[u8]>,
    },
    // For static arrays
    Static(&'static [u8]),
}

/// Compile-time assertion that SmallByteVec is at most 24 bytes
const _: [(); 1] = [(); (std::mem::size_of::<SmallByteVec>() <= 24) as usize];

impl Clone for SmallByteVec {
    /// Clones the `SmallByteVec`, creating a new instance with the same data
    /// or, if stored on the heap, a new reference to the same data.
    ///
    /// This is a very cheap operation.
    fn clone(&self) -> Self {
        Self {
            inner: match &self.inner {
                SmallByteVecInner::Inline { len, data } => {
                    SmallByteVecInner::Inline { len: *len, data: *data }
                }
                SmallByteVecInner::Heap { data } => {
                    SmallByteVecInner::Heap { data: Arc::clone(data) }
                }
                SmallByteVecInner::Static(data) => SmallByteVecInner::Static(data),
            },
        }
    }
}

impl Deref for SmallByteVec {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        match &self.inner {
            SmallByteVecInner::Inline { len, data } => &data[..(*len as usize)],
            SmallByteVecInner::Heap { data } => data,
            SmallByteVecInner::Static(data) => data,
        }
    }
}

impl SmallByteVec {
    /// Construct a new `SmallByteVec` with no data.
    ///
    /// Since `SmallByteVec` is immutable, this is mainly for using it as a
    /// placeholder.
    pub fn new() -> Self {
        SmallByteVec::new_static(&[])
    }

    /// Construct a new `SmallByteVec` from a static slice of bytes.
    ///
    /// ```rust
    /// # use smallarcvec::SmallByteVec;
    /// let data: &'static [u8] = &[1, 2, 3, 4, 5];
    /// let static_vec = SmallByteVec::new_static(data);
    /// assert_eq!(static_vec.len(), 5);
    /// ```
    ///
    /// ```compile_fail
    /// # use smallarcvec::SmallByteVec;
    /// let data = vec![1, 2, 3, 4, 5];
    /// let static_vec = SmallByteVec::new_static(&data); // This will not compile
    /// ```
    pub fn new_static(data: &'static [u8]) -> Self {
        Self { inner: SmallByteVecInner::Static(data) }
    }

    /// Length of the byte vector.
    pub fn len(&self) -> usize {
        match &self.inner {
            SmallByteVecInner::Inline { len, .. } => *len as usize,
            SmallByteVecInner::Heap { data } => data.len(),
            SmallByteVecInner::Static(data) => data.len(),
        }
    }

    /// True if the byte vector is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// True if the byte vector is using inline storage.
    pub fn is_inline(&self) -> bool {
        matches!(self.inner, SmallByteVecInner::Inline { .. })
    }
}

impl Default for SmallByteVec {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: AsRef<[u8]>> From<T> for SmallByteVec {
    fn from(data: T) -> Self {
        let bytes = data.as_ref();
        if bytes.len() <= INLINE_CAPACITY {
            let mut arr = [0u8; INLINE_CAPACITY];
            arr[..bytes.len()].copy_from_slice(bytes);
            Self { inner: SmallByteVecInner::Inline { len: bytes.len() as u8, data: arr } }
        } else {
            Self { inner: SmallByteVecInner::Heap { data: Arc::from(bytes.to_vec()) } }
        }
    }
}

#[cfg(not(tarpaulin_include))]
impl fmt::Debug for SmallByteVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let data = &self[..];
        f.debug_struct("SmallByteVec")
            .field("len", &self.len())
            .field("data", &data)
            .field(
                "storage",
                &match &self.inner {
                    SmallByteVecInner::Inline { .. } => "inline",
                    SmallByteVecInner::Heap { .. } => "heap",
                    SmallByteVecInner::Static(_) => "static",
                },
            )
            .finish()
    }
}

impl PartialEq for SmallByteVec {
    fn eq(&self, other: &Self) -> bool {
        self.deref() == other.deref()
    }
}

impl Eq for SmallByteVec {}

impl Hash for SmallByteVec {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.deref().hash(state);
    }
}

impl PartialOrd for SmallByteVec {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SmallByteVec {
    fn cmp(&self, other: &Self) -> Ordering {
        self.deref().cmp(other.deref())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        let empty = SmallByteVec::from([]);
        assert_eq!(empty.len(), 0);
        assert!(empty.is_empty());
        assert_eq!(&empty[..], &[]);
    }

    #[test]
    fn test_inline() {
        let data = [1, 2, 3, 4, 5];
        let small = SmallByteVec::from(data);
        assert_eq!(small.len(), 5);
        assert!(!small.is_empty());
        assert!(small.is_inline());
        assert_eq!(&small[..], &data[..]);
    }

    #[test]
    fn test_heap() {
        let data: Vec<u8> = (0..30).collect();
        let large = SmallByteVec::from(&data);
        assert_eq!(large.len(), 30);
        assert!(!large.is_empty());
        assert!(!large.is_inline());
        assert_eq!(&large[..], &data[..]);
    }

    #[test]
    fn test_static() {
        let data: &'static [u8] = &[1, 2, 3, 4, 5];
        let static_vec = SmallByteVec::new_static(data);
        assert_eq!(static_vec.len(), 5);
        assert!(!static_vec.is_empty());
        assert!(!static_vec.is_inline());
        assert_eq!(&static_vec[..], data);
    }

    #[test]
    fn test_clone() {
        let data = [1, 2, 3, 4, 5];
        let small = SmallByteVec::from(data);
        let small_clone = small.clone();
        assert_eq!(&small[..], &small_clone[..]);

        let data: Vec<u8> = (0..30).collect();
        let large = SmallByteVec::from(&data);
        let large_clone = large.clone();
        assert_eq!(&large[..], &large_clone[..]);

        let data: &'static [u8] = &[1, 2, 3, 4, 5];
        let static_vec = SmallByteVec::new_static(data);
        let static_clone = static_vec.clone();
        assert_eq!(&static_vec[..], &static_clone[..]);
    }

    #[test]
    fn test_from_vec() {
        let data = vec![1, 2, 3, 4, 5];
        let small = SmallByteVec::from(data.clone());
        assert_eq!(&small[..], &data[..]);

        let data: Vec<u8> = (0..30).collect();
        let large = SmallByteVec::from(data.clone());
        assert_eq!(&large[..], &data[..]);
    }

    #[test]
    fn test_equality() {
        // Same content should be equal regardless of storage type
        let inline1 = SmallByteVec::from([1, 2, 3]);
        let inline2 = SmallByteVec::from([1, 2, 3]);
        let static_vec = SmallByteVec::new_static(&[1, 2, 3]);

        assert_eq!(inline1, inline2);
        assert_eq!(inline1, static_vec);

        // Different content should not be equal
        let different = SmallByteVec::from([1, 2, 4]);
        assert_ne!(inline1, different);

        // Longer vectors
        let long1: Vec<u8> = (0..30).collect();
        let long2: Vec<u8> = (0..30).collect();
        let heap1 = SmallByteVec::from(&long1);
        let heap2 = SmallByteVec::from(&long2);

        assert_eq!(heap1, heap2);

        // Different length vectors
        let short = SmallByteVec::from([0, 1]);
        let long = SmallByteVec::from([0, 1, 2]);
        assert_ne!(short, long);
    }

    #[test]
    fn test_ordering() {
        // Test ordering for same length but different content
        let a = SmallByteVec::from([1, 2, 3]);
        let b = SmallByteVec::from([1, 2, 4]);
        assert!(a < b);

        // Test ordering for different lengths
        let short = SmallByteVec::from([1, 2]);
        let long = SmallByteVec::from([1, 2, 0]);
        assert!(short < long);

        // Test ordering for empty and non-empty
        let empty = SmallByteVec::from([]);
        let non_empty = SmallByteVec::from([1]);
        assert!(empty < non_empty);

        // Test heap stored vectors
        let heap1 = SmallByteVec::from(&(0..25).collect::<Vec<u8>>());
        let heap2 = SmallByteVec::from(&(0..26).collect::<Vec<u8>>());
        assert!(heap1 < heap2);

        // Test mixed storage types
        let inline = SmallByteVec::from([1, 2, 3]);
        let static_vec = SmallByteVec::new_static(&[1, 2, 4]);
        assert!(inline < static_vec);
    }

    #[test]
    fn test_hashing() {
        let data = SmallByteVec::from([1, 2, 3]);
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        data.hash(&mut hasher);
        let hash1 = hasher.finish();

        let data = SmallByteVec::from([2, 3, 4]);
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        data.hash(&mut hasher);
        let hash2 = hasher.finish();
        assert!(hash1 != hash2);
    }

    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_any_vec_contents_preserved(data: Vec<u8>) {
            let small_vec = SmallByteVec::from(data.clone());
            prop_assert_eq!(&small_vec[..], &data[..]);
        }

        #[test]
        fn test_clone_preserves_contents(data: Vec<u8>) {
            let small_vec = SmallByteVec::from(data);
            let cloned = small_vec.clone();
            prop_assert_eq!(&small_vec[..], &cloned[..]);
        }

        #[test]
        fn test_len_matches_original(data: Vec<u8>) {
            let small_vec = SmallByteVec::from(data.clone());
            prop_assert_eq!(small_vec.len(), data.len());
        }

        #[test]
        fn test_is_empty_consistent(data: Vec<u8>) {
            let small_vec = SmallByteVec::from(data.clone());
            prop_assert_eq!(small_vec.is_empty(), data.is_empty());
        }

        #[test]
        fn test_edge_case_inline_capacity(bytes in prop::collection::vec(any::<u8>(), INLINE_CAPACITY)) {
            let small_vec = SmallByteVec::from(bytes.clone());
            prop_assert_eq!(&small_vec[..], &bytes[..]);
            // Should be using inline storage
            match &small_vec.inner {
                SmallByteVecInner::Inline { .. } => {},
                _ => prop_assert!(false, "Should be using inline storage for INLINE_CAPACITY bytes")
            }
        }

        #[test]
        fn test_edge_case_over_inline_capacity(bytes in prop::collection::vec(any::<u8>(), INLINE_CAPACITY + 1)) {
            let small_vec = SmallByteVec::from(bytes.clone());
            prop_assert_eq!(&small_vec[..], &bytes[..]);
            // Should be using heap storage
            match &small_vec.inner {
                SmallByteVecInner::Heap { .. } => {},
                _ => prop_assert!(false, "Should be using heap storage for more than INLINE_CAPACITY bytes")
            }
        }
    }
}
