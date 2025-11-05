use std::{
    fmt,
    ops::{self, Deref},
};

/// The root mean square (RMS) of a set of values.
///
/// RMS is a statistical measure of the magnitude of a varying quantity.
///
/// # Examples
///
/// ```rust
/// # use rastair_types::RootMeanSquare;
/// let data = [1, 2, 3, 4, 5];
/// // Explicitly construct from an iterator
/// let rms = RootMeanSquare::from_iter(data);
/// // our using `collect`
/// let rms: RootMeanSquare = data.into_iter().collect();
/// // you can use the value as a float
/// assert_eq!(rms.round(), 3.0);
/// ```
#[derive(Clone, Copy, Default, serde::Serialize, serde::Deserialize)]
#[must_use]
pub struct RootMeanSquare(f64);

impl<T: Into<f64>> FromIterator<T> for RootMeanSquare {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut count = 0;
        let sum_of_squares: f64 = iter
            .into_iter()
            .map(|x| {
                let x: f64 = x.into();
                count += 1;
                x.powi(2)
            })
            .sum();
        if count == 0 {
            return RootMeanSquare(0.0);
        }
        let average_of_squares = sum_of_squares / f64::from(count);
        RootMeanSquare(average_of_squares.sqrt())
    }
}

impl Deref for RootMeanSquare {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl ops::Mul<f64> for RootMeanSquare {
    type Output = f64;

    fn mul(self, rhs: f64) -> Self::Output {
        self.0 * rhs
    }
}

impl ops::Mul<RootMeanSquare> for f64 {
    type Output = f64;

    fn mul(self, rhs: RootMeanSquare) -> Self::Output {
        self * rhs.0
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Debug for RootMeanSquare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RMS({:.2})", self.0)
    }
}

#[cfg_attr(coverage_nightly, coverage(off))]
impl fmt::Display for RootMeanSquare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// Extension trait for iterators that adds an `.rms()` method.
///
/// This trait provides a convenient way to calculate the root mean square
/// of an iterator's items without explicitly calling `collect()`.
///
/// # Examples
///
/// ```rust
/// # use rastair_types::{RootMeanSquare, RootMeanSquareExt};
/// let data = [1, 2, 3, 4, 5];
/// let rms = data.into_iter().rms();
/// assert_eq!(rms.round(), 3.0);
/// ```
pub trait RootMeanSquareExt: Iterator {
    /// Calculates the root mean square of the iterator's items.
    ///
    /// This is equivalent to collecting into a `RootMeanSquare`.
    fn rms<T>(self) -> RootMeanSquare
    where
        Self: Sized + Iterator<Item = T>,
        T: Into<f64>,
    {
        RootMeanSquare::from_iter(self)
    }
}

impl<I: Iterator> RootMeanSquareExt for I {}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::{collection::vec, prelude::*};

    proptest! {
        #[test]
        fn test_rms_constant_value(value: u8) {
            let data = vec![value; 100];
            let rms = RootMeanSquare::from_iter(data);
            prop_assert_eq!(*rms, f64::from(value));
        }

        #[test]
        fn test_rms_never_negative(data: Vec<u8>) {
            let rms = RootMeanSquare::from_iter(data);
            prop_assert!(*rms >= 0.0);
        }

        #[test]
        fn test_rms_zero_iff_all_zeros(data: Vec<u8>) {
            let rms = RootMeanSquare::from_iter(data.clone());
            if data.iter().all(|&x| x == 0) {
                prop_assert_eq!(rms.0, 0.0);
            } else if !data.is_empty() {
                prop_assert!(*rms > 0.0);
            }
        }

        #[test]
        fn test_rms_greater_than_or_equal_to_mean(data in vec(any::<u8>(), 1..300)) {
            let mean = data.iter().map(|&x| f64::from(x)).sum::<f64>() / data.len() as f64;
            let rms = RootMeanSquare::from_iter(data);
            prop_assert!(*rms >= mean);
        }

        #[test]
        fn test_rms_less_than_or_equal_to_max(data in vec(any::<u8>(), 1..300)) {
            let max = f64::from(*data.iter().max().unwrap());
            let rms = RootMeanSquare::from_iter(data);
            prop_assert!(*rms <= max);
        }
    }
}
