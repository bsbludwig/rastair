use crate::metrics::MetricsForIndel;
use crate::metrics::ml::features::define_features;
use crate::utils::IntoF64 as _;
use rastair_types::Base;

use crate::call::pileup::indels::{IndelAllele, IndelObservation};

define_features! {
    /// Features shared between insertions and deletions.
    pub struct CommonIndelFeatures {
        /// Length of the indel allele in bases.
        scalar indel_len;
        /// Shannon entropy of the allele's base composition (A/C/G/T distribution).
        scalar indel_complexity;
        /// Ratio of this allele's read count to the most frequent allele at this
        /// position. Near 1.0 when this is the dominant allele; low when many
        /// alleles compete.
        scalar indel_dominance;
        /// RMS mapping quality of reads supporting this specific indel allele.
        scalar mapq_rms;
        /// `mapq_rms` divided by the RMS mapping quality of all reads at this
        /// position. Below 1.0 means the indel-supporting reads map worse than the
        /// surrounding pileup (suspicious).
        scalar relative_mapq;
        /// RMS base quality at the pileup anchor position for reads supporting this
        /// allele.
        scalar baseq_rms;
        /// RMS of `min(pos_in_read, read_length - pos_in_read)` across
        /// indel-supporting reads. Low values mean the indel is near read
        /// boundaries (less reliable).
        scalar edge_dist_rms;
        /// Fraction of reads at this position that support this indel allele (VAF).
        /// Coverage-invariant, unlike a raw supporting-read count.
        scalar allele_fraction;
        /// Strand bias: (OT - OB) / (OT + OB). -1.0 = all OB, 0.0 = balanced, +1.0
        /// = all OT.
        scalar strand_bias;
        /// Length of the longest homopolymer run on the reference spanning this
        /// position.
        scalar homopolymer_run;
        /// Fraction of reads covering this position that have a soft-clip in their
        /// CIGAR.
        scalar soft_clip_rate;
        /// Length of the longest dinucleotide repeat (e.g. ATAT) on the reference
        /// spanning this position.
        scalar dinucleotide_run;
        /// Fraction of reads supporting this indel allele that have a homopolymer or
        /// dinucleotide repeat at their read ends (potential alignment artifact).
        scalar repeat_fraction;
    }
}

define_features! {
    pub struct InsertionFeatures {
        /// Features shared with deletions.
        flatten common: CommonIndelFeatures;
        /// RMS base quality of the inserted bases across all reads supporting this
        /// insertion allele.
        scalar insertion_baseq_rms;
        /// Relative baseq vs. entire pileup
        scalar relative_insertion_baseq;
    }
}

define_features! {
    pub struct DeletionFeatures {
        /// Features shared with insertions.
        flatten common: CommonIndelFeatures;
        /// RMS base quality of the first base after the deletion across all reads
        /// supporting this allele. Low quality suggests the deletion boundary is
        /// uncertain. (The anchor base quality is already in [`CommonIndelFeatures`].)
        scalar post_del_baseq_rms;
        /// `post_del_baseq_rms` divided by the RMS base quality of all reads at this
        /// position. Below 1.0 means the deletion boundary is supported by
        /// worse-quality bases than the surrounding pileup.
        scalar relative_post_del_baseq;
    }
}

impl CommonIndelFeatures {
    fn extract(current: &MetricsForIndel) -> CommonIndelFeatures {
        let indel = &current.indel;
        let pileup = &current.metrics.pileup;
        let observations = &pileup.indel_observations;
        let allele = &indel.allele;

        let agg = compute_aggregates(observations, allele);
        let dominance = compute_dominance(&current.metrics.indels.alleles, allele);

        let total_reads = pileup.reads.len().max(1) as f64;
        let soft_clip_rate = pileup.soft_clip_count.f() / total_reads;

        let pos_depth = current.metrics.pos_metrics.depth;
        let pos_mapq = *current.metrics.pos_metrics.mapq;

        CommonIndelFeatures {
            indel_len: allele.len() as f64,
            indel_complexity: allele_entropy(allele.bases()),
            indel_dominance: dominance,
            mapq_rms: agg.mapq_rms,
            relative_mapq: if pos_mapq > 0.0 { agg.mapq_rms / pos_mapq } else { 0.0 },
            baseq_rms: agg.baseq_rms,
            edge_dist_rms: agg.edge_dist_rms,
            allele_fraction: if pos_depth > 0 { agg.total.f() / pos_depth.f() } else { 0.0 },
            strand_bias: strand_bias(agg.ot_count, agg.ob_count),
            homopolymer_run: pileup.homopolymer_run.f(),
            soft_clip_rate,
            dinucleotide_run: pileup.dinucleotide_run.f(),
            repeat_fraction: if agg.total > 0 { agg.repeat_count.f() / agg.total.f() } else { 0.0 },
        }
    }
}

impl InsertionFeatures {
    pub fn extract(current: &MetricsForIndel) -> InsertionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let observations = &current.metrics.pileup.indel_observations;
        let allele = &current.indel.allele;
        let insertion_baseq_rms = insertion_baseq_rms(observations, allele);
        let relative_insertion_baseq = if *current.metrics.pos_metrics.baseq > 0.0 {
            insertion_baseq_rms / *current.metrics.pos_metrics.baseq
        } else {
            0.0
        };
        InsertionFeatures { common, insertion_baseq_rms, relative_insertion_baseq }
    }
}

impl DeletionFeatures {
    pub fn extract(current: &MetricsForIndel) -> DeletionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let observations = &current.metrics.pileup.indel_observations;
        let allele = &current.indel.allele;
        let post_del_baseq_rms = post_del_baseq_rms(observations, allele);
        let relative_post_del_baseq = if *current.metrics.pos_metrics.baseq > 0.0 {
            post_del_baseq_rms / *current.metrics.pos_metrics.baseq
        } else {
            0.0
        };
        DeletionFeatures { common, post_del_baseq_rms, relative_post_del_baseq }
    }
}

struct Aggregates {
    mapq_rms: f64,
    baseq_rms: f64,
    edge_dist_rms: f64,
    ot_count: u32,
    ob_count: u32,
    total: u32,
    repeat_count: u32,
}

fn compute_aggregates(observations: &[IndelObservation], allele: &IndelAllele) -> Aggregates {
    let mut mapq_sq_sum: f64 = 0.0;
    let mut baseq_sq_sum: f64 = 0.0;
    let mut edge_sq_sum: f64 = 0.0;
    let mut ot_count: u32 = 0;
    let mut ob_count: u32 = 0;
    let mut total: u32 = 0;
    let mut repeat_count: u32 = 0;

    for obs in observations {
        if &obs.allele != allele {
            continue;
        }
        total += 1;

        let mq = obs.mapq.f();
        mapq_sq_sum += mq * mq;

        let bq = obs.base_qual.f();
        baseq_sq_sum += bq * bq;

        let edge = (obs.pos_in_read as f64).min((obs.read_length - obs.pos_in_read) as f64);
        edge_sq_sum += edge * edge;

        match obs.strand {
            rastair_types::Strand::OT => ot_count += 1,
            rastair_types::Strand::OB => ob_count += 1,
            rastair_types::Strand::Unknown => {}
        }

        if obs.has_repeat {
            repeat_count += 1;
        }
    }

    let mapq_rms = if total > 0 { (mapq_sq_sum / total.f()).sqrt() } else { 0.0 };
    let baseq_rms = if total > 0 { (baseq_sq_sum / total.f()).sqrt() } else { 0.0 };
    let edge_dist_rms = if total > 0 { (edge_sq_sum / total.f()).sqrt() } else { 0.0 };

    Aggregates { mapq_rms, baseq_rms, edge_dist_rms, ot_count, ob_count, total, repeat_count }
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

fn post_del_baseq_rms(observations: &[IndelObservation], allele: &IndelAllele) -> f64 {
    let mut sq_sum: f64 = 0.0;
    let mut count: u32 = 0;
    for obs in observations {
        if &obs.allele != allele {
            continue;
        }
        if obs.post_del_base_qual > 0 {
            let q = obs.post_del_base_qual.f();
            sq_sum += q * q;
            count += 1;
        }
    }
    if count > 0 { (sq_sum / count.f()).sqrt() } else { 0.0 }
}

fn strand_bias(ot: u32, ob: u32) -> f64 {
    let total = ot + ob;
    if total == 0 {
        return 0.0;
    }
    (ot.f() - ob.f()) / total.f()
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

/// Shannon entropy (bits) of the A/C/G/T composition of the allele's bases.
fn allele_entropy(bases: &[Base]) -> f64 {
    let mut counts = [0u32; 4];
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
    if count == 0 {
        return 0.0;
    }
    counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c.f() / count.f();
            -p * p.log2()
        })
        .sum()
}
