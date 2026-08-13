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
        /// This allele's share of all indel-supporting fragments at the position.
        /// 1.0 when it is the only indel allele; 0.5 when two are tied.
        ///
        /// Deliberately not a ratio to the *most frequent* allele, which is what
        /// this used to be: that scores the winner exactly 1.0 by construction, so
        /// it was constant across the single-allele majority of loci and read the
        /// same for a clean 10-vs-0 as for an ambiguous 10-vs-9.
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
        /// `-log10(p)` of the two-sided exact binomial strand-bias test, against the
        /// strand mix of the *rest* of the locus. 0.0 = exactly as skewed as the
        /// surrounding coverage; higher = more surprising. Clamped to 30.
        ///
        /// Supplies the depth `strand_bias` above lacks: that one is a bare ratio, so
        /// a 2/0 split scores a maximal 1.0 identically to a 20/0 split, though the
        /// first is a coin flip and the second is not. This also conditions on the
        /// locus' own strand skew, so genuinely one-strand coverage does not read as
        /// allele bias.
        ///
        /// Not a gate. As a hard filter at alpha 0.05 this test rejected 4,288 true
        /// chr12 indels to remove 252 false ones, because its null is false for TAPS
        /// — OT and OB reads present different sequence after C→T conversion, so
        /// genuine indel support is strand-asymmetric. It is still real evidence when
        /// weighed against everything else here, which is the point of moving it from
        /// `--indel-strand-bias-alpha` (now defaulted off) into the feature vector.
        scalar strand_bias_log_p;
        /// `ln(1 + depth)` at this position. The vector carried `allele_fraction`
        /// but nothing to scale it by, so a VAF of 0.5 from 2-of-4 fragments was
        /// indistinguishable from 20-of-40. Log-compressed because the difference
        /// between 4x and 40x matters and 400x versus 4000x does not.
        scalar log_depth;
        /// `ln(1 + supporting fragments)`, the other half of `allele_fraction`.
        scalar log_alt_count;
        /// Mean number of indels in the CIGAR of the reads supporting this allele,
        /// the read itself included. Alignments that need several indels to fit are
        /// where spurious ones are placed; a genuine indel usually sits in a read
        /// that needs only it.
        scalar read_indel_burden;
        /// Mean `matching_bases / read_length` over the supporting reads: how much
        /// of each read actually aligned. Low means the support comes from reads
        /// that fit the reference poorly, independent of their mapping quality.
        scalar read_match_fraction;
        /// Fraction of supporting fragments flagged noisy — soft-clipped, or
        /// carrying a terminal tandem repeat.
        ///
        /// The pre-existing `soft_clip_rate` is taken over *every* read at the
        /// position, so it describes the neighbourhood rather than the evidence.
        /// This one is the allele's own, which is what separates a slippage
        /// artifact from a real indel in the same repeat tract.
        scalar noisy_fraction;
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

        let total_reads = pileup.reads.len().max(1).f();
        let soft_clip_rate = pileup.soft_clip_count.f() / total_reads;

        let pos_depth = current.metrics.pos_metrics.depth;
        let pos_mapq = current.metrics.pos_metrics.mapq.f();

        let repeat = RepeatContext::detect(
            allele,
            &pileup.indel_ref_window,
            pileup.indel_ref_anchor as usize,
        );

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
            homopolymer_run: pileup.homopolymer_run.f(),
            soft_clip_rate,
            dinucleotide_run: pileup.dinucleotide_run.f(),
            repeat_fraction: if agg.total > 0 { agg.repeat_count.f() / agg.total.f() } else { 0.0 },
            strand_bias_log_p: strand_bias_log_p(&current.metrics.indels, allele),
            log_depth: pos_depth.f().ln_1p(),
            log_alt_count: agg.total.f().ln_1p(),
            read_indel_burden: mean(agg.indel_burden_sum, agg.total),
            read_match_fraction: mean(agg.match_fraction_sum, agg.total),
            noisy_fraction: if agg.total > 0 { agg.noisy_count.f() / agg.total.f() } else { 0.0 },
        }
    }
}

/// Mean of a sum accumulated over `count` supporting reads; 0.0 when there are
/// none, matching how the RMS aggregates report an empty allele.
fn mean(sum: f32, count: u32) -> f32 {
    if count > 0 { sum / count.f() } else { 0.0 }
}

impl InsertionFeatures {
    pub fn extract(current: &MetricsForIndel) -> InsertionFeatures {
        let common = CommonIndelFeatures::extract(current);
        let observations = &current.metrics.pileup.indel_observations;
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
        let observations = &current.metrics.pileup.indel_observations;
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
    noisy_count: u32,
    indel_burden_sum: f32,
    match_fraction_sum: f32,
}

fn compute_aggregates(observations: &[IndelObservation], allele: &IndelAllele) -> Aggregates {
    let mut mapq_sq_sum: f32 = 0.0;
    let mut baseq_sq_sum: f32 = 0.0;
    let mut edge_sq_sum: f32 = 0.0;
    let mut ot_count: u32 = 0;
    let mut ob_count: u32 = 0;
    let mut total: u32 = 0;
    let mut repeat_count: u32 = 0;
    let mut noisy_count: u32 = 0;
    let mut indel_burden_sum: f32 = 0.0;
    let mut match_fraction_sum: f32 = 0.0;

    for obs in observations {
        if &obs.allele != allele {
            continue;
        }
        total += 1;

        let mq = obs.mapq.f();
        mapq_sq_sum += mq * mq;

        let bq = obs.base_qual.f();
        baseq_sq_sum += bq * bq;

        // `saturating_sub`: these two come from different htslib accessors (query
        // position vs stored SEQ length), so nothing in the type system stops
        // `pos_in_read` exceeding `read_length`. Plain `-` wraps in release and
        // turns one malformed record into an edge distance of ~4e9, which drags
        // the RMS for the whole allele.
        let edge = obs.pos_in_read.min(obs.read_length.saturating_sub(obs.pos_in_read)).f();
        edge_sq_sum += edge * edge;

        indel_burden_sum += obs.num_indels_in_read.f();
        if obs.read_length > 0 {
            match_fraction_sum += obs.matching_bases.f() / obs.read_length.f();
        }
        if obs.noisy {
            noisy_count += 1;
        }

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

    Aggregates {
        mapq_rms,
        baseq_rms,
        edge_dist_rms,
        ot_count,
        ob_count,
        total,
        repeat_count,
        noisy_count,
        indel_burden_sum,
        match_fraction_sum,
    }
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

/// `-log10` of this allele's strand-bias p-value, or 0.0 when the allele is not
/// among the locus' counted alleles (nothing to test).
///
/// Clamped to 30, which is far past any p a real pileup produces and keeps the
/// value finite for `p = 0` underflow. The counts come from
/// [`IndelCounts`](crate::call::pileup::indels::IndelCounts) rather than from
/// `compute_aggregates` so the test and its null are drawn from the same
/// fragment-level bookkeeping.
fn strand_bias_log_p(
    indels: &crate::call::pileup::indels::IndelCounts,
    allele: &IndelAllele,
) -> f32 {
    let Some(counts) = indels.alleles.iter().find(|a| &a.allele == allele) else {
        return 0.0;
    };
    let p = counts.strand_bias_p_value(indels.null_ot_fraction(counts));
    (-p.max(1e-30).log10() as f32).clamp(0.0, 30.0)
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
    let this_count = alleles.iter().find(|a| &a.allele == target).map(|a| a.total()).unwrap_or(0);
    this_count.f() / total.f()
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

    use crate::call::pileup::indels::IndelAlleleCounts;

    fn allele_counts(seq: &str, total: u32) -> IndelAlleleCounts {
        let mut counts = IndelAlleleCounts::new(IndelAllele::Insertion(bases(seq)));
        counts.ot = total;
        counts
    }

    /// Dominance is a share of all indel support, not a ratio to the winner. The
    /// old form scored whichever allele led exactly 1.0 by construction, so it was
    /// constant wherever only one allele existed — the common case — and could not
    /// tell a clean call from a coin flip between two alleles.
    #[test]
    fn dominance_is_a_share_not_a_ratio_to_the_winner() {
        let contested = [allele_counts("A", 10), allele_counts("T", 9)];
        let lead = compute_dominance(&contested, &IndelAllele::Insertion(bases("A")));
        assert!((lead - 10.0 / 19.0).abs() < 1e-6, "leading allele of a near-tie: {lead}");

        let alone = [allele_counts("A", 10)];
        assert_eq!(compute_dominance(&alone, &IndelAllele::Insertion(bases("A"))), 1.0);

        // An allele that is not counted at this locus has no share of it.
        assert_eq!(compute_dominance(&alone, &IndelAllele::Insertion(bases("G"))), 0.0);
        assert_eq!(compute_dominance(&[], &IndelAllele::Insertion(bases("A"))), 0.0);
    }

    /// `pos_in_read` and `read_length` come from different htslib accessors, so a
    /// malformed record can put the position past the end. The subtraction must
    /// saturate: wrapping turns one record into an edge distance of ~4e9 and
    /// destroys the RMS for the whole allele.
    #[test]
    fn edge_distance_saturates_past_the_read_end() {
        let mut obs = IndelObservation {
            allele: IndelAllele::Insertion(bases("A")),
            strand: seqair_types::Strand::OT,
            reverse: false,
            pos_in_read: 150,
            read_length: 100,
            mapq: 60,
            base_qual: 30,
            matching_bases: 100,
            num_indels_in_read: 1,
            insertion_base_quals: SmallVec::new(),
            post_del_base_qual: 0,
            has_repeat: false,
            noisy: false,
        };
        let agg = compute_aggregates(std::slice::from_ref(&obs), &obs.allele.clone());
        assert_eq!(agg.edge_dist_rms, 0.0, "position past the end clamps to distance 0");

        obs.pos_in_read = 10;
        let agg = compute_aggregates(std::slice::from_ref(&obs), &obs.allele.clone());
        assert_eq!(agg.edge_dist_rms, 10.0, "min(10, 100 - 10)");
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
