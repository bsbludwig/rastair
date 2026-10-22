//! Indel records from the pair-HMM's verdict: no model, no binomial genotyper,
//! no strand gate.
//!
//! [`call_on_pair_hmm`] writes the records of the events anchored at one
//! position and no other. Two called events at one anchor make one
//! [`IndelRecords::Compound`] `1/2` record (always PASS, QUAL the weaker
//! margin's); otherwise one biallelic record per event, its genotype from the
//! called pair's copies. A scored event not called stays a `0/0` call, written
//! only under `--all` with `low_indel_qual` (read off the genotype, not the
//! margin). Emission and `phasing::problem::indel_site` match on
//! [`IndelRecords`] rather than re-deriving it from the genotypes.

use super::GenotypeTag;
use crate::call::{
    pair_hmm::{EventIdx, HaplotypeScores, MAX_EVENTS},
    pileup::indels::IndelAllele,
};
use better_default::Default;
use seqair_types::{Phred, SmallVec};
use std::num::NonZeroU8;
use tracing::debug;

/// CLI parameters for indel calling.
#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize, Default)]
pub struct IndelParams {
    /// Enable experimental indel calling
    ///
    /// Indels are called by a TAPS-aware pair-HMM.
    #[arg(long)]
    #[arg(help_heading = crate::utils::cli::sections::INDELS)]
    pub experimental_indels: bool,

    /// Minimum alternate observations for an indel allele to be scored
    ///
    /// Reads (one per fragment) whose alignment carries the allele; at least
    /// 5% of the column's depth is required in any case. Of the
    /// alleles reaching it in a cluster of neighbouring columns, the three
    /// most reads carry are scored, at most two per column, and a scored
    /// allele is called only if the pair-HMM finds it present.
    #[arg(long, default_value_t = 2)]
    #[arg(help_heading = crate::utils::cli::sections::INDELS)]
    #[default(2)]
    pub min_indel_ao: u32,
}

/// Result of indel calling at a single position for one allele.
#[derive(Debug, Clone)]
pub struct IndelCall {
    pub allele: IndelAllele,
    pub genotype: GenotypeTag,
    /// The pair-HMM's margin of presence for the allele; its `QUAL`.
    pub quality: Phred,
    /// The genotype's quality; its `GC`.
    pub genotype_quality: Phred,
    /// The event the allele is in the cluster's `HaplotypeScores`.
    pub event: EventIdx,
}

impl IndelCall {
    /// Whether the genotype holds the allele; a call that does not is written
    /// only under `--all`, filtered.
    #[must_use]
    pub fn is_called(&self) -> bool {
        self.genotype != GenotypeTag::HomRef
    }
}

/// A column's indel records.
#[derive(Debug, Clone)]
pub enum IndelRecords {
    /// Both alleles called, `1/2`: one multi-allelic record, as two biallelic
    /// ones would each carry a genotype referring to an allele they do not
    /// declare.
    Compound([IndelCall; 2]),
    /// One biallelic record per scored allele, in record order.
    Separate(SmallVec<IndelCall, MAX_EVENTS>),
}

impl IndelRecords {
    #[must_use]
    pub fn calls(&self) -> &[IndelCall] {
        match self {
            Self::Compound(calls) => calls,
            Self::Separate(calls) => calls,
        }
    }
}

/// Call the indels anchored at `pos` on the pair-HMM's scores of their
/// cluster.
///
/// An event is called when its margin reaches
/// [`Margin::PRESENT`](crate::call::pair_hmm::Margin::PRESENT). Two called
/// events at one anchor make a compound het, genotyped by the cluster's called
/// pair ([`crate::call::pair_hmm::Verdict::called_pair`]), as is any other
/// called event, its `GC` the event's own GQ
/// ([`crate::call::pair_hmm::Verdict::dosage`]), as SNVs have. A scored event that is not
/// called stays as a `0/0` call. There is no model and no strand gate; on
/// chr12 no strand test did better than none. An allele the pair-HMM did not
/// score, past the cluster's leading events or without a read to score it on,
/// is not called: nothing speaks for it.
pub fn call_on_pair_hmm(scores: &HaplotypeScores, pos: u32) -> IndelRecords {
    const ONE: NonZeroU8 = NonZeroU8::new(1).expect("nonzero");
    const TWO: NonZeroU8 = NonZeroU8::new(2).expect("nonzero");
    let verdict = &scores.verdict;
    let call = |event: EventIdx, genotype: GenotypeTag, quality: Phred, genotype_quality| {
        let allele = scores.event(event)?.allele.as_indel()?.clone();
        Some(IndelCall { allele, genotype, quality, genotype_quality, event })
    };

    // Events are in record order, so a compound het's two alleles are ALT 1
    // and 2 as they stand.
    let here: SmallVec<EventIdx, MAX_EVENTS> = scores.indel_events_at(pos).collect();
    // The called pair, not the dosage: each event's own dosage can say two
    // copies, and two events at one anchor cannot both be on both chromosomes.
    let compound = || {
        let [first, second] = *here.as_ref() else { return None };
        if verdict.copies_of(first) == 0 || verdict.copies_of(second) == 0 {
            return None;
        }
        // One record standing for both alleles is as sure as the less certain.
        let quality = verdict.margin(first)?.min(verdict.margin(second)?).qual();
        let genotype = GenotypeTag::alt_het(ONE, TWO);
        Some([
            call(first, genotype, quality, verdict.quality)?,
            call(second, genotype, quality, verdict.quality)?,
        ])
    };
    if let Some(calls) = compound() {
        return IndelRecords::Compound(calls);
    }

    let separate = here.iter().filter_map(|&event| {
        let quality = verdict.margin(event)?.qual();
        let genotype = match verdict.copies_of(event) {
            0 => return call(event, GenotypeTag::hom_ref(), quality, verdict.quality),
            1 => GenotypeTag::ref_het(ONE),
            _ => GenotypeTag::hom_alt(ONE),
        };
        let genotype_quality = verdict.dosage(event).map_or_else(
            || {
                debug!(pos, ?event, "no dosage for a called indel, GC is the pair's");
                verdict.quality
            },
            |dosage| dosage.gq,
        );
        call(event, genotype, quality, genotype_quality)
    });
    IndelRecords::Separate(separate.collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::pair_hmm::fixtures::*;
    use seqair_types::Base;

    const ALT_1: NonZeroU8 = NonZeroU8::new(1).expect("nonzero");
    const ALT_2: NonZeroU8 = NonZeroU8::new(2).expect("nonzero");

    fn deletion(seq: &str) -> IndelAllele {
        IndelAllele::Deletion(seq.bytes().map(Base::from).collect())
    }

    /// `alleles` scored in record order, haplotype `i + 1` carrying
    /// `alleles[i]`, over `reads` of each kind.
    fn scores(alleles: &[&str], reads: &[([f32; 3], usize)]) -> HaplotypeScores {
        HaplotypeScores::at_column(
            alleles.iter().map(|a| deletion(a)).collect(),
            rows(reads).into(),
        )
        .expect("a verdict")
    }

    fn on_pair_hmm(scores: &HaplotypeScores) -> Vec<(IndelAllele, GenotypeTag)> {
        call_on_pair_hmm(scores, 0).calls().iter().map(|c| (c.allele.clone(), c.genotype)).collect()
    }

    /// A het, a homozygote and the reference, each as the best pair has it,
    /// with `QUAL` from the margin and a called allele's `GC` from its
    /// dosage; under the reference
    /// the allele stays as a `0/0` call for `--all` to show.
    #[test]
    fn a_called_allele_takes_its_copies_in_the_best_pair() {
        for (reads, expected) in [
            (
                vec![(REF_READ, 4), (ALT1_READ, 4)],
                vec![(deletion("A"), GenotypeTag::RefHet(ALT_1))],
            ),
            (vec![(ALT1_READ, 8)], vec![(deletion("A"), GenotypeTag::HomAlt(ALT_1))]),
            (vec![(REF_READ, 8)], vec![(deletion("A"), GenotypeTag::HomRef)]),
        ] {
            let scores = scores(&["A"], &reads);
            assert_eq!(on_pair_hmm(&scores), expected);
            let verdict = &scores.verdict;
            for call in call_on_pair_hmm(&scores, 0).calls() {
                assert_eq!(call.event, EventIdx::new(0).expect("event"));
                let margin = verdict.margin(call.event).expect("a margin");
                assert_eq!(*call.quality, *margin.qual());
                let gc = match verdict.dosage(call.event) {
                    Some(dosage) if call.is_called() => dosage.gq,
                    _ => verdict.quality,
                };
                assert_eq!(*call.genotype_quality, *gc);
            }
        }
    }

    /// Two called alleles in the best pair are one compound record: both
    /// calls `1/2`, in the alleles' order, at the weaker allele's `QUAL`.
    #[test]
    fn two_called_alleles_make_a_compound_het() {
        let scores = scores(&["AA", "A"], &[(ALT1_READ, 9), (ALT2_READ, 4)]);
        let IndelRecords::Compound(calls) = call_on_pair_hmm(&scores, 0) else {
            panic!("a compound record");
        };
        let compound = GenotypeTag::AltHet(ALT_1, ALT_2);
        let alleles: Vec<_> = calls.iter().map(|c| (c.allele.clone(), c.genotype)).collect();
        assert_eq!(alleles, [(deletion("AA"), compound), (deletion("A"), compound)]);
        let [first, second] = [0, 1]
            .map(|i| scores.verdict.margin(EventIdx::new(i).expect("event")).expect("margin"));
        assert!(calls.iter().all(|c| *c.quality == *first.min(second).qual()));
    }

    /// The best pair holds the second allele, but too weakly to call it, so
    /// the second is a `0/0` and the first takes the zygosity the reads give
    /// it without the second: a homozygote when no read carries the
    /// reference, a het when one does.
    #[test]
    fn an_uncalled_half_of_a_compound_leaves_the_other_alone() {
        let hom = scores(&["A", "AA"], &[([-9.0, -1.0, -9.0], 2), ([-9.0, -1.5, -1.0], 4)]);
        assert_eq!(
            on_pair_hmm(&hom),
            [(deletion("A"), GenotypeTag::HomAlt(ALT_1)), (deletion("AA"), GenotypeTag::HomRef)]
        );
        let het = scores(
            &["A", "AA"],
            &[([-1.0, -9.0, -2.0], 1), ([-9.0, -1.0, -9.0], 2), ([-9.0, -1.5, -1.0], 2)],
        );
        assert_eq!(
            on_pair_hmm(&het),
            [(deletion("A"), GenotypeTag::RefHet(ALT_1)), (deletion("AA"), GenotypeTag::HomRef)]
        );
    }

    /// Eight reads carry an insertion at 5 and a deletion at 20; one lacks
    /// the deletion and, ending at the insertion, barely prefers the
    /// reference there. The pair puts that read on the second chromosome, so
    /// the insertion is a het, and its `GC` is the dosage's.
    #[test]
    fn a_linked_indel_takes_the_pairs_dosage() {
        use crate::call::pair_hmm::{Cluster, Event, ReadScores};
        const NO: f32 = f32::NEG_INFINITY;
        let events = [
            Event { anchor: 5, allele: IndelAllele::Insertion([Base::T].into()).into() },
            Event { anchor: 20, allele: deletion("A").into() },
        ];
        let cluster = Cluster::new(events.into_iter().collect()).expect("a cluster");
        // By haplotype: reference, insertion, deletion, both.
        let both = ReadScores::new(None, 60, [-9.0, -4.5, -4.5, 0.0, NO, NO, NO, NO]);
        let other = ReadScores::new(None, 60, [0.0, -0.3, -4.5, -4.5, NO, NO, NO, NO]);
        let reads: Vec<_> = std::iter::repeat_n(both, 8).chain(std::iter::once(other)).collect();
        let scores = HaplotypeScores::new(cluster, reads.into()).expect("scores");
        let insertion = EventIdx::new(0).expect("event");
        let dosage = scores.verdict.dosage(insertion).expect("present");

        let records = call_on_pair_hmm(&scores, 5);
        let [call] = records.calls() else { panic!("one record") };
        assert_eq!(call.genotype, GenotypeTag::RefHet(ALT_1));
        assert_eq!(dosage.copies, 1);
        assert_eq!(*call.genotype_quality, *dosage.gq);
    }

    const A: u32 = 10;
    const B: u32 = 20;

    /// Two events at anchor [`A`] and one at [`B`], over `n` reads each that
    /// only the haplotype carrying the events `haps` lists explains (as bitmasks
    /// over the three events).
    fn cluster(haps: &[(usize, usize)]) -> HaplotypeScores {
        use crate::call::pair_hmm::{Cluster, Event, MAX_HAPLOTYPES, ReadScores};
        let at = |anchor, seq| Event { anchor, allele: deletion(seq).into() };
        let events = [at(A, "A"), at(A, "AA"), at(B, "C")].into_iter().collect();
        let reads: Vec<ReadScores> = haps
            .iter()
            .flat_map(|&(hap, n)| {
                let mut log10 = [-9.0; MAX_HAPLOTYPES];
                log10[hap] = -1.0;
                std::iter::repeat_n(ReadScores::new(None, 60, log10), n)
            })
            .collect();
        let cluster = Cluster::new(events).expect("a cluster");
        HaplotypeScores::new(cluster, reads.into()).expect("a verdict")
    }

    fn calls_at(scores: &HaplotypeScores, pos: u32) -> (bool, Vec<(IndelAllele, GenotypeTag)>) {
        let records = call_on_pair_hmm(scores, pos);
        let calls = records.calls().iter().map(|c| (c.allele.clone(), c.genotype)).collect();
        (matches!(records, IndelRecords::Compound(_)), calls)
    }

    /// A column's records hold the events anchored at it and no other: the
    /// two at [`A`] on different chromosomes are one compound record, the one
    /// at [`B`], on the first chromosome with an event at [`A`], is its own
    /// `0/1`.
    #[test]
    fn a_column_calls_the_events_of_its_cluster_anchored_at_it() {
        let scores = cluster(&[(0b101, 6), (0b010, 6)]);
        let compound = GenotypeTag::AltHet(ALT_1, ALT_2);
        assert_eq!(
            calls_at(&scores, A),
            (true, vec![(deletion("A"), compound), (deletion("AA"), compound)])
        );
        assert_eq!(
            calls_at(&scores, B),
            (false, vec![(deletion("C"), GenotypeTag::RefHet(ALT_1))])
        );
        assert_eq!(calls_at(&scores, B + 1), (false, vec![]));
    }

    /// An event of the cluster the called pair does not carry is a `0/0`
    /// call at its own column, for `--all` to show.
    #[test]
    fn an_absent_event_of_a_cluster_is_a_hom_ref_call() {
        let scores = cluster(&[(0b001, 6), (0b010, 6)]);
        assert_eq!(calls_at(&scores, B), (false, vec![(deletion("C"), GenotypeTag::HomRef)]));
        assert!(calls_at(&scores, A).0);
    }
}
