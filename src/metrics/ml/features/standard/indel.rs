use crate::metrics::MetricsForIndel;
use crate::metrics::ml::features::utils::one_hot_encode_base;
use crate::utils::IntoF64 as _;
use color_eyre::Result;
use rastair_types::Base;

use crate::call::pileup::indels::{IndelAllele, IndelObservation};

pub struct CommonIndelFeatures {
    indel_len: f64,
    indel_complexity: f64,
    indel_base_count: [f64; 4],
    indel_dominance: f64,
    mapq_rms: f64,
    mapq0_rate: f64,
    baseq_rms: f64,
    edge_dist_rms: f64,
    depth: f64,
    strand_balance: f64,
    ctx_before_2: [f64; 4],
    ctx_before_1: [f64; 4],
    ctx_after_1: [f64; 4],
    ctx_after_2: [f64; 4],
    homopolymer_run: f64,
}

pub struct InsertionFeatures {
    common: CommonIndelFeatures,
    insertion_baseq_rms: f64,
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
    buf[common_len] = features.insertion_baseq_rms;
    Ok(())
}

pub fn deletion(
    current: &MetricsForIndel,
    buf: &mut [f64; DeletionFeatures::FEATURES],
) -> Result<()> {
    let features = DeletionFeatures::extract(current);
    let common_len = CommonIndelFeatures::FEATURES;
    features.common.write_to(&mut buf[..common_len]);
    buf[common_len..].copy_from_slice(&features.ref_one_hot);
    Ok(())
}

impl CommonIndelFeatures {
    pub const FEATURES: usize = size_of::<Self>() / size_of::<f64>();

    fn write_to(&self, buf: &mut [f64]) {
        buf[0] = self.indel_len;
        buf[1] = self.indel_complexity;
        buf[2..6].copy_from_slice(&self.indel_base_count);
        buf[6] = self.indel_dominance;
        buf[7] = self.mapq_rms;
        buf[8] = self.mapq0_rate;
        buf[9] = self.baseq_rms;
        buf[10] = self.edge_dist_rms;
        buf[11] = self.depth;
        buf[12] = self.strand_balance;
        buf[13..17].copy_from_slice(&self.ctx_before_2);
        buf[17..21].copy_from_slice(&self.ctx_before_1);
        buf[21..25].copy_from_slice(&self.ctx_after_1);
        buf[25..29].copy_from_slice(&self.ctx_after_2);
        buf[29] = self.homopolymer_run;
    }

    fn extract(current: &MetricsForIndel) -> CommonIndelFeatures {
        let indel = &current.indel;
        let pileup = &current.metrics.pileup;
        let observations = &pileup.indel_observations;
        let allele = &indel.allele;

        let agg = compute_aggregates(observations, allele);
        let count_and_entropy = CountAndEntropy::from_bases(allele.bases());
        let dominance = compute_dominance(&current.metrics.indels.alleles, allele);

        let ctx = &pileup.context;
        let (b2a, b2c, b2g, b2t) = one_hot_encode_base(ctx.before_2);
        let (b1a, b1c, b1g, b1t) = one_hot_encode_base(ctx.before_1);
        let (a1a, a1c, a1g, a1t) = one_hot_encode_base(ctx.after_1);
        let (a2a, a2c, a2g, a2t) = one_hot_encode_base(ctx.after_2);

        CommonIndelFeatures {
            indel_len: allele.len() as f64,
            indel_complexity: count_and_entropy.entropy,
            indel_base_count: [
                count_and_entropy.counts[0].f(),
                count_and_entropy.counts[1].f(),
                count_and_entropy.counts[2].f(),
                count_and_entropy.counts[3].f(),
            ],
            indel_dominance: dominance,
            mapq_rms: agg.mapq_rms,
            mapq0_rate: if agg.total > 0 { agg.mapq0_count.f() / agg.total.f() } else { 0.0 },
            baseq_rms: agg.baseq_rms,
            edge_dist_rms: agg.edge_dist_rms,
            depth: agg.total.f(),
            strand_balance: strand_balance(agg.fwd_count, agg.rev_count),
            ctx_before_2: [b2a, b2c, b2g, b2t],
            ctx_before_1: [b1a, b1c, b1g, b1t],
            ctx_after_1: [a1a, a1c, a1g, a1t],
            ctx_after_2: [a2a, a2c, a2g, a2t],
            homopolymer_run: pileup.homopolymer_run.f(),
        }
    }
}

impl InsertionFeatures {
    pub const FEATURES: usize = size_of::<Self>() / size_of::<f64>();

    fn extract(current: &MetricsForIndel) -> InsertionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let observations = &current.metrics.pileup.indel_observations;
        let allele = &current.indel.allele;
        let insertion_baseq_rms = insertion_baseq_rms(observations, allele);
        InsertionFeatures { common, insertion_baseq_rms }
    }
}

impl DeletionFeatures {
    pub const FEATURES: usize = size_of::<Self>() / size_of::<f64>();

    fn extract(current: &MetricsForIndel) -> DeletionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let ref_one_hot = one_hot_encode_base(current.metrics.pileup.reference_base).into();
        DeletionFeatures { common, ref_one_hot }
    }
}

struct Aggregates {
    mapq_rms: f64,
    mapq0_count: u32,
    baseq_rms: f64,
    edge_dist_rms: f64,
    fwd_count: u32,
    rev_count: u32,
    total: u32,
}

fn compute_aggregates(observations: &[IndelObservation], allele: &IndelAllele) -> Aggregates {
    let mut mapq_sq_sum: f64 = 0.0;
    let mut mapq0_count: u32 = 0;
    let mut baseq_sq_sum: f64 = 0.0;
    let mut edge_sq_sum: f64 = 0.0;
    let mut fwd_count: u32 = 0;
    let mut rev_count: u32 = 0;
    let mut total: u32 = 0;

    for obs in observations {
        if &obs.allele != allele {
            continue;
        }
        total += 1;

        let mq = obs.mapq.f();
        mapq_sq_sum += mq * mq;
        if obs.mapq == 0 {
            mapq0_count += 1;
        }

        let bq = obs.base_qual.f();
        baseq_sq_sum += bq * bq;

        let edge = (obs.pos_in_read as f64).min((obs.read_length - obs.pos_in_read) as f64);
        edge_sq_sum += edge * edge;

        if obs.reverse {
            rev_count += 1;
        } else {
            fwd_count += 1;
        }
    }

    let mapq_rms = if total > 0 { (mapq_sq_sum / total.f()).sqrt() } else { 0.0 };
    let baseq_rms = if total > 0 { (baseq_sq_sum / total.f()).sqrt() } else { 0.0 };
    let edge_dist_rms = if total > 0 { (edge_sq_sum / total.f()).sqrt() } else { 0.0 };

    Aggregates { mapq_rms, mapq0_count, baseq_rms, edge_dist_rms, fwd_count, rev_count, total }
}

fn insertion_baseq_rms(observations: &[IndelObservation], allele: &IndelAllele) -> f64 {
    let mut sq_sum: f64 = 0.0;
    let mut count: u32 = 0;
    for obs in observations {
        if &obs.allele != allele {
            continue;
        }
        for &q in &obs.insertion_base_quals {
            let qf = q.f();
            sq_sum += qf * qf;
            count += 1;
        }
    }
    if count > 0 { (sq_sum / count.f()).sqrt() } else { 0.0 }
}

fn strand_balance(fwd: u32, rev: u32) -> f64 {
    if fwd == 0 || rev == 0 {
        return 0.0;
    }
    let min = fwd.min(rev).f();
    let max = fwd.max(rev).f();
    min / max
}

fn compute_dominance(
    alleles: &[crate::call::pileup::indels::IndelAlleleCounts],
    target: &IndelAllele,
) -> f64 {
    let total: u32 = alleles.iter().map(|a| a.total()).sum();
    if total == 0 {
        return 0.0;
    }
    let max_count = alleles.iter().map(|a| a.total()).max().unwrap_or(0);
    let this_count = alleles.iter().find(|a| &a.allele == target).map(|a| a.total()).unwrap_or(0);
    if max_count > 0 { this_count.f() / max_count.f() } else { 0.0 }
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
