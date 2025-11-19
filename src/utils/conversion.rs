use rastair_types::{Probability, RootMeanSquare};

/// Shortcut to get the default value of a type.
pub fn default<T: Default>() -> T {
    T::default()
}

pub trait IntoF64 {
    /// Dangerously convert this into f64!
    #[track_caller]
    fn f(self) -> f64;
}

macro_rules! impl_into_f64 {
    ($($t:ty),* $(,)?) => {
        $(
            impl IntoF64 for $t {
                #[track_caller]
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::cast_lossless, reason = "macro code")]
                fn f(self) -> f64 {
                    self as f64
                }
            }
        )*
    };
}

impl_into_f64!(u8, u16, u32, u64, u128, usize, i8, i16, i32, i64, i128, isize, f32);

impl IntoF64 for RootMeanSquare {
    #[track_caller]
    fn f(self) -> f64 {
        *self
    }
}

impl IntoF64 for Probability {
    #[track_caller]
    fn f(self) -> f64 {
        *self
    }
}
