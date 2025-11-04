use rastair_types::RootMeanSquare;

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
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
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
