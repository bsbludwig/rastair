use crate::{
    utils::Base,
    vcf::{ByStrand, Record},
};

pub fn one_hot_encode_base(base: impl Into<Base>) -> (f64, f64, f64, f64) {
    match base.into() {
        Base::A => (1., 0., 0., 0.),
        Base::C => (0., 1., 0., 0.),
        Base::G => (0., 0., 1., 0.),
        Base::T => (0., 0., 0., 1.),
        Base::Unknown => (0., 0., 0., 0.),
    }
}

pub fn get_strand_base_quality(record: &Record, base: Base) -> ByStrand<f64> {
    record
        .info
        .strand_specific_base_quality
        .iter()
        .find(|x| x.base == base)
        .copied()
        .unwrap_or_default()
}

pub fn get_strand_map_quality(record: &Record, base: Base) -> ByStrand<f64> {
    record
        .info
        .strand_specific_mapping_quality
        .iter()
        .find(|x| x.base == base)
        .copied()
        .unwrap_or_default()
}
