//! Case counts for property tests that pin their own `test_cases`.
//!
//! A count compiled into `#[hegel::test(test_cases = N)]` wins over
//! `HEGEL_TEST_CASES`, so that variable does not reach the pinned tests. They
//! are written `test_cases = pinned::cases(N)` instead, and
//! `RASTAIR_PINNED_CASES` overrides them all.
//!
//! Shared with the integration tests via `#[path = "../src/pinned.rs"]`.

/// `n`, or the value of `RASTAIR_PINNED_CASES` when it is set.
pub fn cases(n: u64) -> u64 {
    match std::env::var("RASTAIR_PINNED_CASES") {
        Ok(v) if !v.is_empty() => {
            v.parse().expect("RASTAIR_PINNED_CASES must be a positive integer")
        }
        _ => n,
    }
}
