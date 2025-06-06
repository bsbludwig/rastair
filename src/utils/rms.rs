use std::{fmt, ops::Deref};

/// The root mean square (RMS) of a set of values.
///
/// RMS is a statistical measure of the magnitude of a varying quantity.
///
/// # Examples
///
/// ```rust
/// # use rastair2::utils::RootMeanSquare;
/// let data = [1, 2, 3, 4, 5];
/// // constructing the type calculates the value
/// let rms = RootMeanSquare::new(&data);
/// // you can use the value as a float
/// assert_eq!(rms.round(), 3.0);
/// ```
#[derive(Clone, Copy)]
pub struct RootMeanSquare(f64);

impl RootMeanSquare {
    pub fn new<T: Copy + Into<f64>>(data: &[T]) -> RootMeanSquare {
        if data.is_empty() {
            return RootMeanSquare(0.0);
        }
        let sum_of_squares: f64 = data
            .iter()
            .map(|x| {
                let x: f64 = (*x).into();
                x.powi(2)
            })
            .sum();
        let average_of_squares = sum_of_squares / data.len() as f64;
        RootMeanSquare(average_of_squares.sqrt())
    }
}

impl Deref for RootMeanSquare {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[cfg(not(tarpaulin_include))]
impl fmt::Debug for RootMeanSquare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RMS({:.2})", self.0)
    }
}

#[cfg(not(tarpaulin_include))]
impl fmt::Display for RootMeanSquare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::{collection::vec, prelude::*};

    proptest! {
        #[test]
        fn test_rms_constant_value(value: u8) {
            let data = vec![value; 100];
            let rms = RootMeanSquare::new(&data);
            prop_assert_eq!(*rms, f64::from(value));
        }

        #[test]
        fn test_rms_never_negative(data: Vec<u8>) {
            let rms = RootMeanSquare::new(&data);
            prop_assert!(*rms >= 0.0);
        }

        #[test]
        fn test_rms_zero_iff_all_zeros(data: Vec<u8>) {
            let rms = RootMeanSquare::new(&data);
            if data.iter().all(|&x| x == 0) {
                prop_assert_eq!(rms.0, 0.0);
            } else if !data.is_empty() {
                prop_assert!(*rms > 0.0);
            }
        }

        #[test]
        fn test_rms_greater_than_or_equal_to_mean(data in vec(any::<u8>(), 1..300)) {
            let mean = data.iter().map(|&x| f64::from(x)).sum::<f64>() / data.len() as f64;
            let rms = RootMeanSquare::new(&data);
            prop_assert!(*rms >= mean);
        }

        #[test]
        fn test_rms_less_than_or_equal_to_max(data in vec(any::<u8>(), 1..300)) {
            let max = f64::from(*data.iter().max().unwrap());
            let rms = RootMeanSquare::new(&data);
            prop_assert!(*rms <= max);
        }
    }
}
