use seqair_types::Base;

pub fn one_hot_encode_base(base: impl Into<Base>) -> (f64, f64, f64, f64) {
    match base.into() {
        Base::A => (1., 0., 0., 0.),
        Base::C => (0., 1., 0., 0.),
        Base::G => (0., 0., 1., 0.),
        Base::T => (0., 0., 0., 1.),
        Base::Unknown => (0., 0., 0., 0.),
    }
}

/// Safe division that returns 0.0 instead of NaN when denominator is 0
#[inline]
pub fn safe_div(numerator: f64, denominator: f64) -> f64 {
    if denominator == 0.0 { 0.0 } else { numerator / denominator }
}
