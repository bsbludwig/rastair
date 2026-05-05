use crate::metrics::MetricsForIndel;
use crate::metrics::ml::features::utils::one_hot_encode_base;
use crate::{metrics::PileupMetrics, utils::IntoF64 as _};
use color_eyre::Result;
use rastair_types::{Base, RootMeanSquare};

pub struct CommonIndelFeatures {
    mapq: RootMeanSquare,
    mapq0: f64,
    read_complexity: RootMeanSquare,
    position_in_read: RootMeanSquare,
    num_aligned_bases: RootMeanSquare,
    num_indels_in_read: RootMeanSquare,
    indel_len: f64,
    indel_complexity: f64,
    indel_base_count: [f64; 4],
}

pub struct InsertionFeatures {
    common: CommonIndelFeatures,
    insertion_baseq: RootMeanSquare,
}

pub struct DeletionFeatures {
    common: CommonIndelFeatures,
    ref_one_hot: [f64; 4],
}

pub fn insertion(
    current: &MetricsForIndel,
    buf: &mut [f64; InsertionFeatures::FEATURES],
) -> Result<()> {
    let features = InsertionFeatures::extract(current);
    let common_len = CommonIndelFeatures::FEATURES;
    features.common.write_to(&mut buf[..common_len]);
    buf[CommonIndelFeatures::FEATURES] = features.insertion_baseq.f();
    Ok(())
}

pub fn deletion(
    current: &MetricsForIndel,
    buf: &mut [f64; DeletionFeatures::FEATURES],
) -> Result<()> {
    let features = DeletionFeatures::extract(current);
    let common_len = CommonIndelFeatures::FEATURES;
    features.common.write_to(&mut buf[..common_len]);
    buf[CommonIndelFeatures::FEATURES..].copy_from_slice(&features.ref_one_hot);
    Ok(())
}

impl CommonIndelFeatures {
    pub const FEATURES: usize = size_of::<Self>() / size_of::<f64>();

    fn write_to(&self, buf: &mut [f64]) {
        buf[0] = self.mapq.f();
        buf[1] = self.mapq0;
        buf[2] = self.read_complexity.f();
        buf[3] = self.position_in_read.f();
        buf[4] = self.num_aligned_bases.f();
        buf[5] = self.num_indels_in_read.f();
        buf[6] = self.indel_len;
        buf[7] = self.indel_complexity;
        buf[8..12].copy_from_slice(&self.indel_base_count);
    }

    fn extract(current: &MetricsForIndel) -> CommonIndelFeatures {
        let indel = &current.indel;
        let PileupMetrics { pos_metrics: pos, ref_metrics: r, .. } = &current.metrics;
        let CountAndEntropy { counts, entropy } = CountAndEntropy::from_bases(indel.allele.bases());

        CommonIndelFeatures {
            mapq: r.mapq,
            mapq0: pos.mapq0.f(),
            read_complexity: RootMeanSquare::default(), // TODO: calculate RMS read complexity
            position_in_read: r.position_in_read,
            num_aligned_bases: r.num_aligned_bases,
            num_indels_in_read: r.num_indels,
            indel_len: indel.allele.len() as f64,
            indel_complexity: entropy,
            indel_base_count: [counts[0].f(), counts[1].f(), counts[2].f(), counts[3].f()],
        }
    }
}

impl InsertionFeatures {
    pub const FEATURES: usize = size_of::<Self>() / size_of::<f64>();

    fn extract(current: &MetricsForIndel) -> InsertionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let insertion_baseq = RootMeanSquare::default(); // TODO: calculate RMS base quality of inserted bases
        InsertionFeatures { common, insertion_baseq }
    }
}

impl DeletionFeatures {
    pub const FEATURES: usize = size_of::<Self>() / size_of::<f64>();

    fn extract(current: &MetricsForIndel) -> DeletionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let ref_one_hot = one_hot_encode_base(current.metrics.ref_base()).try_into().unwrap();
        DeletionFeatures { common, ref_one_hot }
    }
}

struct CountAndEntropy {
    counts: [u32; 4],
    entropy: f64,
}

impl CountAndEntropy {
    fn from_bases(bases: &[Base]) -> Self {
        let mut counts = [0; 4];
        for &b in bases {
            match b {
                Base::A => counts[0] += 1,
                Base::C => counts[1] += 1,
                Base::G => counts[2] += 1,
                Base::T => counts[3] += 1,
                _ => {}
            }
        }
        let count: u32 = counts.iter().sum();
        let entropy = if count > 0 {
            counts
                .iter()
                .filter(|&&c| c > 0)
                .map(|&c| {
                    let p = (c.f()) / (count.f());
                    -p * p.log2()
                })
                .sum()
        } else {
            0.0
        };
        Self { counts, entropy }
    }
}
