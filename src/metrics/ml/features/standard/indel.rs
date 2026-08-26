use crate::metrics::MetricsForIndel;
use crate::metrics::ml::features::define_features;
use crate::utils::IntoF32 as _;
use seqair_types::Base;

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
        /// 1.0 if this indel is a tandem-repeat slippage event — its bases are
        /// whole copies of the repeat unit of the reference tract it sits in
        /// (homopolymer or short tandem repeat). The dominant indel artifact mode.
        scalar indel_in_repeat;
        /// Length in bases of the reference tandem-repeat tract the indel slips
        /// within (0 when not in a repeat). Allele-aware generalisation of
        /// `homopolymer_run`: longer tracts slip more often.
        scalar repeat_tract_length;
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
        let metrics = &current.metrics;
        let d = metrics.indel_data.as_ref().expect("indel data present for indel call position");
        let observations = &d.observations;
        let allele = &indel.allele;

        let agg = compute_aggregates(observations, allele);
        let dominance = compute_dominance(&d.counts.alleles, allele);

        let total_reads = metrics.pos_metrics.depth.max(1).f();
        let soft_clip_rate = d.soft_clip_count.f() / total_reads;

        let pos_depth = metrics.pos_metrics.depth;
        let pos_mapq = metrics.pos_metrics.mapq.f();

        let repeat = RepeatContext::detect(allele, &d.ref_window, d.ref_anchor as usize);

        CommonIndelFeatures {
            indel_len: allele.len().f(),
            indel_complexity: allele_entropy(allele.bases()),
            indel_dominance: dominance,
            mapq_rms: agg.mapq_rms,
            relative_mapq: if pos_mapq > 0.0 { agg.mapq_rms / pos_mapq } else { 0.0 },
            baseq_rms: agg.baseq_rms,
            edge_dist_rms: agg.edge_dist_rms,
            allele_fraction: if pos_depth > 0 { agg.total.f() / pos_depth.f() } else { 0.0 },
            strand_bias: strand_bias(agg.ot_count, agg.ob_count),
            indel_in_repeat: if repeat.in_repeat { 1.0 } else { 0.0 },
            repeat_tract_length: repeat.tract_length.f(),
            homopolymer_run: d.homopolymer_run.f(),
            soft_clip_rate,
            dinucleotide_run: d.dinucleotide_run.f(),
            repeat_fraction: if agg.total > 0 { agg.repeat_count.f() / agg.total.f() } else { 0.0 },
        }
    }
}

impl InsertionFeatures {
    pub fn extract(current: &MetricsForIndel) -> InsertionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let d = current
            .metrics
            .indel_data
            .as_ref()
            .expect("indel data present for indel call position");
        let observations = &d.observations;
        let allele = &current.indel.allele;
        let insertion_baseq_rms = insertion_baseq_rms(observations, allele);
        let pos_baseq = current.metrics.pos_metrics.baseq.f();
        let relative_insertion_baseq =
            if pos_baseq > 0.0 { insertion_baseq_rms / pos_baseq } else { 0.0 };
        InsertionFeatures { common, insertion_baseq_rms, relative_insertion_baseq }
    }
}

impl DeletionFeatures {
    pub fn extract(current: &MetricsForIndel) -> DeletionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let d = current
            .metrics
            .indel_data
            .as_ref()
            .expect("indel data present for indel call position");
        let observations = &d.observations;
        let allele = &current.indel.allele;
        let post_del_baseq_rms = post_del_baseq_rms(observations, allele);
        let pos_baseq = current.metrics.pos_metrics.baseq.f();
        let relative_post_del_baseq =
            if pos_baseq > 0.0 { post_del_baseq_rms / pos_baseq } else { 0.0 };
        DeletionFeatures { common, post_del_baseq_rms, relative_post_del_baseq }
    }
}

struct Aggregates {
    mapq_rms: f32,
    baseq_rms: f32,
    edge_dist_rms: f32,
    ot_count: u32,
    ob_count: u32,
    total: u32,
    repeat_count: u32,
}

fn compute_aggregates(observations: &[IndelObservation], allele: &IndelAllele) -> Aggregates {
    let mut mapq_sq_sum: f32 = 0.0;
    let mut baseq_sq_sum: f32 = 0.0;
    let mut edge_sq_sum: f32 = 0.0;
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

        let edge = (obs.pos_in_read.f()).min((obs.read_length - obs.pos_in_read).f());
        edge_sq_sum += edge * edge;

        match obs.strand {
            seqair_types::Strand::OT => ot_count += 1,
            seqair_types::Strand::OB => ob_count += 1,
            seqair_types::Strand::Unknown => {}
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

fn insertion_baseq_rms(observations: &[IndelObservation], allele: &IndelAllele) -> f32 {
    let mut sq_sum: f32 = 0.0;
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

fn post_del_baseq_rms(observations: &[IndelObservation], allele: &IndelAllele) -> f32 {
    let mut sq_sum: f32 = 0.0;
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

fn strand_bias(ot: u32, ob: u32) -> f32 {
    let total = ot + ob;
    if total == 0 {
        return 0.0;
    }
    (ot.f() - ob.f()) / total.f()
}

fn compute_dominance(
    alleles: &[crate::call::pileup::indels::IndelAlleleCounts],
    target: &IndelAllele,
) -> f32 {
    let total: u32 = alleles.iter().map(|a| a.total()).sum();
    if total == 0 {
        return 0.0;
    }
    let max_count = alleles.iter().map(|a| a.total()).max().unwrap_or(0);
    let this_count = alleles.iter().find(|a| &a.allele == target).map(|a| a.total()).unwrap_or(0);
    if max_count > 0 { this_count.f() / max_count.f() } else { 0.0 }
}

/// Shannon entropy (bits) of the A/C/G/T composition of the allele's bases.
fn allele_entropy(bases: &[Base]) -> f32 {
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

/// Tandem-repeat slippage descriptor for an indel allele, derived from the
/// local reference window.
struct RepeatContext {
    /// The indel's bases are whole copies of the repeat unit of the reference
    /// tract it sits in (homopolymer or short tandem repeat).
    in_repeat: bool,
    /// Length in bases of that reference tandem-repeat tract (0 if not in one).
    tract_length: u32,
}

impl RepeatContext {
    fn none() -> Self {
        Self { in_repeat: false, tract_length: 0 }
    }

    /// `window[anchor]` is the anchor base; the inserted/deleted bases of
    /// `allele` follow it (start at `anchor + 1` in reference coordinates).
    fn detect(allele: &IndelAllele, window: &[Base], anchor: usize) -> Self {
        let bases = allele.bases();
        let unit = repeat_unit(bases);
        if bases.is_empty() || unit.is_empty() || anchor >= window.len() {
            return Self::none();
        }

        // Reference immediately downstream of the anchor — where inserted bases
        // land, and what deleted bases are drawn from.
        let downstream = window.get(anchor + 1..).unwrap_or(&[]);
        let mut copies = count_unit_copies(downstream, unit);

        // Homopolymer tracts are phase-free, so extend across the anchor and
        // upstream to capture the full run length regardless of indel alignment.
        if unit.len() == 1 {
            let base = unit.first().copied();
            let mut i = anchor as isize;
            while i >= 0 && window.get(i as usize).copied() == base {
                copies += 1;
                i -= 1;
            }
        }

        let units_in_indel = bases.len() / unit.len();
        let in_repeat = match allele {
            // An insertion that duplicates an already-present unit expands the tract.
            IndelAllele::Insertion(_) => copies >= 1,
            // Deleted bases are reference; it's slippage when they are part of a
            // ≥2-copy array (extra units flank, or the deleted span itself repeats).
            IndelAllele::Deletion(_) => copies > units_in_indel || units_in_indel >= 2,
        };

        Self { in_repeat, tract_length: u32::try_from(copies * unit.len()).unwrap_or(u32::MAX) }
    }
}

/// The minimal repeating unit of `bases` (`AAAA`→`A`, `ATAT`→`AT`, `ACG`→`ACG`).
/// Returns the whole slice when it has no shorter period.
fn repeat_unit(bases: &[Base]) -> &[Base] {
    let n = bases.len();
    for p in 1..=n {
        if n.is_multiple_of(p) && (p..n).all(|i| bases.get(i) == bases.get(i - p)) {
            return bases.get(..p).unwrap_or(bases);
        }
    }
    bases
}

/// How many consecutive copies of `unit` appear at the start of `seq`.
fn count_unit_copies(seq: &[Base], unit: &[Base]) -> usize {
    if unit.is_empty() {
        return 0;
    }
    let mut copies = 0;
    while seq.get(copies * unit.len()..(copies + 1) * unit.len()) == Some(unit) {
        copies += 1;
    }
    copies
}

#[cfg(test)]
mod tests {
    use super::*;
    use seqair_types::SmallVec;

    fn bases(s: &str) -> SmallVec<Base, 4> {
        s.bytes().map(Base::from).collect()
    }
    fn window(s: &str) -> Vec<Base> {
        s.bytes().map(Base::from).collect()
    }

    #[test]
    fn repeat_unit_minimal_period() {
        assert_eq!(repeat_unit(&bases("AAAA")), &bases("A")[..]);
        assert_eq!(repeat_unit(&bases("ATAT")), &bases("AT")[..]);
        assert_eq!(repeat_unit(&bases("ACG")), &bases("ACG")[..]);
        // period 2 does not divide length 5, so there is no shorter period
        assert_eq!(repeat_unit(&bases("ATATA")), &bases("ATATA")[..]);
    }

    #[test]
    fn count_copies_counts_leading_units_only() {
        assert_eq!(count_unit_copies(&window("AAAAC"), &bases("A")), 4);
        assert_eq!(count_unit_copies(&window("ATATATG"), &bases("AT")), 3);
        assert_eq!(count_unit_copies(&window("CAAA"), &bases("A")), 0);
    }

    #[test]
    fn homopolymer_insertion_is_slippage() {
        // C[anchor=0] A A A A A — insert A
        let win = window("CAAAAA");
        let r = RepeatContext::detect(&IndelAllele::Insertion(bases("A")), &win, 0);
        assert!(r.in_repeat);
        assert_eq!(r.tract_length, 5);
    }

    #[test]
    fn homopolymer_deletion_is_slippage() {
        let win = window("CAAAAA");
        let r = RepeatContext::detect(&IndelAllele::Deletion(bases("A")), &win, 0);
        assert!(r.in_repeat);
        assert_eq!(r.tract_length, 5);
    }

    #[test]
    fn dinucleotide_insertion_is_slippage() {
        // C[anchor=0] A T A T A T — insert AT
        let win = window("CATATAT");
        let r = RepeatContext::detect(&IndelAllele::Insertion(bases("AT")), &win, 0);
        assert!(r.in_repeat);
        assert_eq!(r.tract_length, 6);
    }

    #[test]
    fn non_repeat_insertion_is_not_slippage() {
        let win = window("CGGGC");
        let r = RepeatContext::detect(&IndelAllele::Insertion(bases("A")), &win, 0);
        assert!(!r.in_repeat);
        assert_eq!(r.tract_length, 0);
    }

    #[test]
    fn single_copy_deletion_is_not_slippage() {
        // delete the lone A in C[anchor]ACGT — one copy, no flanking repeat
        let win = window("CACGT");
        let r = RepeatContext::detect(&IndelAllele::Deletion(bases("A")), &win, 0);
        assert!(!r.in_repeat);
    }

    #[test]
    fn homopolymer_run_extends_across_anchor() {
        // A A [anchor=2]A A A C — insert A; tract spans upstream + downstream
        let win = window("AAAAAC");
        let r = RepeatContext::detect(&IndelAllele::Insertion(bases("A")), &win, 2);
        assert!(r.in_repeat);
        assert_eq!(r.tract_length, 5);
    }
}
