use crate::utils::{Base, Strand};
use probability::prelude::{Binomial, Discrete as _};
use seqair_types::SmallVec;

/// How long a tandem repeat at either read terminus has to be before the read is
/// treated as slippage-prone, expressed in whole repeat units per period.
///
/// Units rather than bases, so the threshold means the same thing for every
/// period: three units of a dinucleotide is six bases, three units of a
/// homopolymer is three. Getting this wrong is expensive in a way that is easy to
/// miss — a shared *base* window of 3 makes the period-2 check reduce to
/// `seq[0] == seq[2] || seq[len - 3] == seq[len - 1]`, which fires on 43.75% of
/// random reads and so waters down whatever the flag feeds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TerminalRepeatLimits {
    /// Minimum homopolymer length, in bases.
    pub homopolymer_units: usize,
    /// Minimum dinucleotide repeat length, in 2 bp units.
    pub dinucleotide_units: usize,
}

impl Default for TerminalRepeatLimits {
    /// Chosen so a flagged read is unusual rather than typical: on random
    /// sequence a 4 bp terminal homopolymer occurs ~3% of the time and a 3 unit
    /// (6 bp) dinucleotide repeat ~0.8%, against 43.75% for the 3 bp shared
    /// window these replaced.
    fn default() -> Self {
        Self { homopolymer_units: 4, dinucleotide_units: 3 }
    }
}

/// A specific indel allele observed in reads.
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum IndelAllele {
    /// Bases inserted after the anchor position (not including the anchor).
    Insertion(SmallVec<Base, 4>),
    /// Reference bases deleted after the anchor position.
    Deletion(SmallVec<Base, 4>),
}

impl IndelAllele {
    pub fn len(&self) -> usize {
        match self {
            Self::Insertion(bases) | Self::Deletion(bases) => bases.len(),
        }
    }

    pub fn bases(&self) -> &[Base] {
        match self {
            Self::Insertion(bases) | Self::Deletion(bases) => bases,
        }
    }

    pub fn is_insertion(&self) -> bool {
        matches!(self, Self::Insertion(_))
    }

    pub fn is_deletion(&self) -> bool {
        matches!(self, Self::Deletion(_))
    }
}

/// A single read's indel observation at a pileup position.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndelObservation {
    pub allele: IndelAllele,
    pub strand: Strand,
    pub reverse: bool,
    pub pos_in_read: u32,
    pub read_length: u32,
    pub mapq: u8,
    pub base_qual: u8,
    pub matching_bases: u32,
    pub num_indels_in_read: u32,
    pub insertion_base_quals: SmallVec<u8, 4>,
    pub post_del_base_qual: u8,
    pub has_repeat: bool,
    /// Terminal tandem repeat or soft-clip: this fragment's alignment is the kind
    /// that slips. Excluded from the hard-filter pathway's genotyping counts on
    /// both sides of the ratio — see [`IndelCounts::noisy_ref_count`].
    pub noisy: bool,
}

/// Aggregated indel counts at a position, ready for calling.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct IndelCounts {
    /// Each unique indel allele with its OT/OB strand counts.
    pub alleles: SmallVec<IndelAlleleCounts, 2>,
    /// Fragments at this position on the original-top strand, supporting any
    /// allele or none. Together with [`Self::ob_depth`] this provides the null
    /// for the strand-bias test; see [`IndelCounts::null_ot_fraction`].
    pub ot_depth: u32,
    /// Fragments at this position on the original-bottom strand.
    pub ob_depth: u32,
    /// Fragments at this position that do NOT support an indel.
    ///
    /// Together with [`Self::total_indel_reads`] this partitions the fragments in
    /// `Pileup::reads`: both sides are counted after the same read filters,
    /// coverage cap and overlap deduplication.
    pub ref_count: u32,
    /// Fragments with a terminal homopolymer/dinucleotide repeat, subtracted from
    /// depth on the ML pathway.
    pub depth_offset: u32,
    /// Non-supporting fragments that look noisy — a terminal tandem repeat or a
    /// soft-clip. The reference-side half of the hard-filter pathway's noise
    /// exclusion; the alt-side half is [`IndelAlleleCounts::noisy`].
    ///
    /// Both halves are dropped together by [`Self::clean_depth`], so the noise
    /// exclusion cannot move the alt/depth ratio on its own. Subtracting only
    /// this side would inflate VAF by exactly the noise rate — which is the whole
    /// reason the two are tracked separately from the one-sided
    /// [`Self::depth_offset`] the ML pathway still uses.
    pub noisy_ref_count: u32,
}

impl IndelCounts {
    pub fn total_indel_reads(&self) -> u32 {
        self.alleles.iter().map(|a| a.total()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }

    /// Fragments at this position with neither a terminal tandem repeat nor a
    /// soft-clip, over all alleles and the reference alike.
    ///
    /// The denominator the hard-filter pathway genotypes against. Every allele's
    /// [`IndelAlleleCounts::clean_total`] is a subset of it by construction, so
    /// the binomial always sees `alt <= depth`.
    pub fn clean_depth(&self) -> u32 {
        let clean_ref = self.ref_count.saturating_sub(self.noisy_ref_count);
        clean_ref + self.alleles.iter().map(|a| a.clean_total()).sum::<u32>()
    }

    /// Expected OT share for `allele` under the null that it is drawn from the
    /// same strand mix as the rest of the locus.
    ///
    /// Taken from the fragments *not* supporting this allele, so that an allele
    /// cannot dilute the null with its own skew — this is the binomial analogue
    /// of conditioning on the margins in a 2x2 strand-bias table. A locus whose
    /// coverage is genuinely strand-skewed therefore does not make every allele
    /// on it look biased.
    ///
    /// Smoothed with one fragment's worth of prior mass, so the estimate cannot
    /// reach exactly 0 or 1 when the background fragments all happen to fall on
    /// one strand — routine at low depth, not an edge case, and a degenerate null
    /// makes a single fragment on the other strand infinitely surprising.
    ///
    /// The prior is centred on the locus' *whole* strand mix rather than on 0.5,
    /// which is what makes the estimate behave as the background thins out. A
    /// homozygous indel has no background at all, and a flat 0.5 prior would then
    /// judge it against balanced coverage it never had: at a locus covered only on
    /// OT, every hom-alt allele is an n/0 split and gets rejected for a skew that
    /// belongs to the coverage, not to the allele. Falling back to the locus mix
    /// instead states the right null — an allele that *is* the locus cannot be
    /// shown to be skewed relative to it. With a substantial background the prior
    /// is swamped and this is the background fraction either way.
    pub fn null_ot_fraction(&self, allele: &IndelAlleleCounts) -> f64 {
        let locus_ot = f64::from(self.ot_depth);
        let locus_ob = f64::from(self.ob_depth);
        let locus_fraction = (locus_ot + 0.5) / (locus_ot + locus_ob + 1.0);

        let ot = f64::from(self.ot_depth.saturating_sub(allele.ot));
        let ob = f64::from(self.ob_depth.saturating_sub(allele.ob));
        (ot + locus_fraction) / (ot + ob + 1.0)
    }
}

/// Fragments supporting one indel allele, split by the strand of the original
/// DNA duplex they came from.
///
/// The split is by OT/OB, not by the alignment's reverse flag. Both mates of a
/// fragment share an OT/OB assignment but have opposite reverse flags, so under
/// per-fragment deduplication (the default) a reverse-flag split would collapse
/// to whichever mate happened to survive — inside the mate-overlap window every
/// fragment would look single-stranded and no allele would ever look strand
/// balanced. OT/OB is invariant under that deduplication, and is the axis
/// strand-specific artifacts actually fall along in TAPS data.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IndelAlleleCounts {
    pub allele: IndelAllele,
    /// Supporting fragments on the original-top strand.
    pub ot: u32,
    /// Supporting fragments on the original-bottom strand.
    pub ob: u32,
    /// Supporting fragments whose strand could not be determined. Counted in
    /// [`Self::total`] but never evidence for or against strand bias.
    ///
    /// Unreachable on the `call` path as it stands: both orientation modes
    /// (`flags` and `--guess-read-orientation`) always resolve to OT or OB, and
    /// `Pileup::from_hts` drops an alignment whose orientation is `None`. Kept
    /// because [`Strand`] is a shared type whose third variant other entry points
    /// (`call-reads`, `.mpk` written by another build) do produce, and because
    /// "no orientation" must not silently read as "OT" in a strand-bias test.
    pub unknown_strand: u32,
    /// Supporting fragments excluded from genotyping as noise — see
    /// [`IndelCounts::noisy_ref_count`] for the reference-side counterpart.
    ///
    /// Counted in [`Self::total`] (it is real support, and the strand-bias test
    /// uses all of it) but not in [`Self::clean_total`].
    pub noisy: u32,
}

impl IndelAlleleCounts {
    pub fn new(allele: IndelAllele) -> Self {
        Self { allele, ot: 0, ob: 0, unknown_strand: 0, noisy: 0 }
    }

    /// Records one supporting fragment.
    pub fn add(&mut self, strand: Strand, noisy: bool) {
        match strand {
            Strand::OT => self.ot += 1,
            Strand::OB => self.ob += 1,
            Strand::Unknown => self.unknown_strand += 1,
        }
        if noisy {
            self.noisy += 1;
        }
    }

    pub fn total(&self) -> u32 {
        self.ot + self.ob + self.unknown_strand
    }

    /// Supporting fragments that are not noise: the numerator matching
    /// [`IndelCounts::clean_depth`].
    pub fn clean_total(&self) -> u32 {
        self.total().saturating_sub(self.noisy)
    }

    /// Two-sided exact binomial p-value for this allele's OT/OB split against
    /// `null_ot_fraction`: the probability of a split at least this lopsided if
    /// the supporting fragments were drawn from the locus' own strand mix.
    ///
    /// Low means strand biased. Replaces an earlier `ot > 0 && ob > 0` rule,
    /// which at the default `min_indel_ao` of 2 rejected a 2/0 split — an event
    /// with probability 0.5 under the null, and so no evidence of anything.
    /// Measured on real data that rule tracked the chance rate almost exactly
    /// below AO≈4; see `.claude/notes/indel-strand-concordance-vs-fragment-dedup.md`.
    ///
    /// [`Self::unknown_strand`] fragments are excluded: they are evidence
    /// neither for nor against bias. An allele supported only by fragments of
    /// undetermined orientation therefore scores 1.0 (not biased) rather than
    /// being rejected — absence of strand information is not evidence of skew.
    pub fn strand_bias_p_value(&self, null_ot_fraction: f64) -> f64 {
        two_sided_binomial_p(self.ot, self.ot + self.ob, null_ot_fraction)
    }
}

/// Two-sided exact binomial test of `successes` out of `trials` against `p`.
///
/// Uses the min-tail doubling convention rather than summing all outcomes at
/// most as likely as the observed one: the two agree exactly at p=0.5 (the
/// dominant case here, since sequencing coverage is near strand balanced) and
/// doubling is the more conservative of the two off-centre, which is the safer
/// direction for a filter that gates variant output.
fn two_sided_binomial_p(successes: u32, trials: u32, p: f64) -> f64 {
    if trials == 0 {
        return 1.0;
    }
    // `Binomial::new` panics outside the open interval, and a NaN would compare
    // false against the alpha threshold and silently disable the filter.
    // Callers go through `IndelCounts::null_ot_fraction`, which is smoothed away
    // from the boundaries; this only has to keep a bad caller from panicking.
    const BOUND: f64 = 1e-9;
    let p = if p.is_finite() { p.clamp(BOUND, 1.0 - BOUND) } else { 0.5 };
    let dist = Binomial::new(trials as usize, p);
    let k = successes as usize;

    let lower: f64 = (0..=k).map(|i| dist.mass(i)).sum();
    let upper: f64 = (k..=trials as usize).map(|i| dist.mass(i)).sum();

    (2.0 * lower.min(upper)).min(1.0)
}

#[cfg(test)]
mod strand_bias_tests {
    use super::*;

    fn insertion(ot: u32, ob: u32) -> IndelAlleleCounts {
        let mut counts = IndelAlleleCounts::new(IndelAllele::Insertion(SmallVec::new()));
        counts.ot = ot;
        counts.ob = ob;
        counts
    }

    /// The whole point of the change: at the default `min_indel_ao` of 2, a 2/0
    /// split is a coin flip and must not read as evidence of bias. The values
    /// below are exact — `2 * 0.5^n` — so this doubles as a check on the
    /// statistic itself.
    #[test]
    fn one_sided_support_p_value_halves_with_each_fragment() {
        for (ot, expected) in [(1, 1.0), (2, 0.5), (3, 0.25), (4, 0.125), (6, 0.03125)] {
            let p = insertion(ot, 0).strand_bias_p_value(0.5);
            assert!((p - expected).abs() < 1e-9, "{ot}/0 split: expected p={expected}, got {p}");
        }
    }

    #[test]
    fn balanced_support_is_never_biased() {
        for n in 1..=20 {
            assert_eq!(insertion(n, n).strand_bias_p_value(0.5), 1.0);
        }
    }

    #[test]
    fn p_value_is_symmetric_between_strands() {
        for (ot, ob) in [(7, 1), (9, 2), (12, 0)] {
            assert_eq!(
                insertion(ot, ob).strand_bias_p_value(0.5),
                insertion(ob, ot).strand_bias_p_value(0.5),
                "{ot}/{ob} and {ob}/{ot} are equally lopsided"
            );
        }
    }

    /// Against a skewed null, support matching that skew is unremarkable while
    /// the same split against a balanced null is significant.
    #[test]
    fn null_fraction_shifts_what_counts_as_biased() {
        let allele = insertion(8, 0);
        assert!(allele.strand_bias_p_value(0.5) < 0.05);
        assert!(allele.strand_bias_p_value(0.9) > 0.05);
    }

    /// Unknown-orientation support carries no strand information, so an allele
    /// made entirely of it cannot be shown to be biased.
    #[test]
    fn support_without_orientation_is_not_biased() {
        let mut allele = insertion(0, 0);
        allele.unknown_strand = 12;
        assert_eq!(allele.total(), 12);
        assert_eq!(allele.strand_bias_p_value(0.5), 1.0);
    }

    /// A degenerate null must not produce a NaN that silently compares false
    /// against the alpha threshold.
    #[test]
    fn degenerate_null_stays_finite() {
        for null in [0.0, 1.0, -1.0, 2.0, f64::NAN] {
            let p = insertion(4, 1).strand_bias_p_value(null);
            assert!((0.0..=1.0).contains(&p), "null {null} produced p={p}");
        }
    }
}
