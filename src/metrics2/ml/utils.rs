use rastair_types::Base;
pub fn one_hot_encode_base(base: impl Into<Base>) -> (f64, f64, f64, f64) {
    match base.into() {
        Base::A => (1., 0., 0., 0.),
        Base::C => (0., 1., 0., 0.),
        Base::G => (0., 0., 1., 0.),
        Base::T => (0., 0., 0., 1.),
        Base::Unknown => (0., 0., 0., 0.),
    }
}
