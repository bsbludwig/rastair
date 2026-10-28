//! Indel records from the pair-HMM's verdict: no model, no binomial genotyper,
//! no strand gate.
//!
//! [`call_on_pair_hmm`] writes the records of the events anchored at one
//! position and no other. Two events at one anchor the called pair carries
//! make one [`IndelRecords::Compound`] record (always `1/2` and PASS, QUAL the
//! weaker margin's); otherwise one biallelic record per event. A scored event the
//! called pair leaves out stays a `0/0` call, written only under `--all` with
//! `low_indel_qual` (read off the genotype, not the margin). The called pair
//! decides only which events are called; a called record's `GT`, and every
//! record's `GQ` and `PL`, come from its own genotype likelihoods, as the
//! **Genotypes** section of `call::pair_hmm` describes. Emission and
//! `phasing::problem::indel_site` match on [`IndelRecords`] rather than
//! re-deriving it from the genotypes.

use super::GenotypeTag;
use crate::call::{
    pair_hmm::{EventIdx, GenotypeLikelihoods, HaplotypeScores, MAX_EVENTS},
    pileup::{from_seqair::CandidateFloor, indels::IndelAllele},
};
use better_default::Default;
use seqair_types::{Phred, Pos0, SmallVec};
use std::num::NonZeroU8;

/// CLI parameters for indel calling.
#[derive(Debug, Clone, clap::Args, serde::Serialize, serde::Deserialize, Default)]
pub struct IndelParams {
    /// Minimum alternate observations for an indel allele to be scored
    ///
    /// Reads (one per fragment) whose alignment carries the allele; at least
    /// 5% of the column's depth is required in any case. Of the
    /// alleles reaching it in a cluster of neighbouring columns, the three
    /// most reads carry are scored, at most two per column, and a scored
    /// allele is called only if the pair-HMM finds it present.
    #[arg(long, default_value_t = CandidateFloor { min_reads: 2 })]
    #[arg(help_heading = crate::utils::cli::sections::INDELS)]
    #[default(CandidateFloor { min_reads: 2 })]
    pub min_indel_ao: CandidateFloor,
}

/// Result of indel calling at a single position for one allele.
#[derive(Debug, Clone)]
pub struct IndelCall {
    pub allele: IndelAllele,
    pub genotype: GenotypeTag,
    /// The event the allele is in the cluster's `HaplotypeScores`.
    pub event: EventIdx,
}

/// What a record says about its genotype: one per record, whether it holds one
/// allele or a compound het's two.
#[derive(Debug, Clone)]
pub struct RecordScore {
    /// The pair-HMM's margin of presence for the allele (the weaker of a
    /// compound's two); the record's `QUAL`.
    pub quality: Phred,
    /// The record's genotype likelihoods, over its alleles in ALT order.
    pub likelihoods: GenotypeLikelihoods,
}

impl RecordScore {
    /// The `GQ` of `genotype`.
    pub fn genotype_quality(&self, genotype: GenotypeTag) -> Phred {
        let [a, b] = genotype.alleles();
        self.likelihoods.quality(a, b)
    }
}

/// A column's indel records.
#[derive(Debug, Clone)]
pub enum IndelRecords {
    /// Both alleles called, always `1/2`: one multi-allelic record, as two
    /// biallelic ones would each carry a genotype referring to an allele
    /// they do not declare.
    Compound { calls: [IndelCall; 2], score: RecordScore },
    /// One biallelic record per scored allele, in record order.
    Separate(SmallVec<BiallelicRecord, MAX_EVENTS>),
}

/// One allele's record.
#[derive(Debug, Clone)]
pub struct BiallelicRecord {
    pub call: IndelCall,
    pub score: RecordScore,
}

impl IndelRecords {
    /// Every call with its record's score; a compound's two share theirs.
    pub fn scored_calls(&self) -> impl Iterator<Item = (&IndelCall, &RecordScore)> {
        let (compound, separate) = match self {
            Self::Compound { calls, score } => (Some((calls, score)), None),
            Self::Separate(calls) => (None, Some(calls)),
        };
        compound
            .into_iter()
            .flat_map(|(calls, score)| calls.iter().map(move |call| (call, score)))
            .chain(separate.into_iter().flatten().map(|record| (&record.call, &record.score)))
    }
}

/// Call the indels anchored at `pos` on the pair-HMM's scores of their
/// cluster.
///
/// An event is called when the cluster's called pair
/// ([`crate::call::pair_hmm::Verdict::called_pair`]) carries it; its record's
/// genotype is the likeliest of its likelihoods carrying an alt. Two called
/// events at one anchor make one compound record. A scored event that is not
/// called stays as a `0/0` call.
/// There is no model and no strand gate; on chr12 no strand test did better
/// than none. An allele the pair-HMM did not
/// score, past the cluster's leading events or without a read to score it on,
/// is not called: nothing speaks for it.
pub fn call_on_pair_hmm(scores: &HaplotypeScores, pos: Pos0) -> IndelRecords {
    const ONE: NonZeroU8 = NonZeroU8::new(1).expect("nonzero");
    const TWO: NonZeroU8 = NonZeroU8::new(2).expect("nonzero");
    let verdict = &scores.verdict();
    let call = |event: EventIdx, genotype: GenotypeTag| {
        let allele = scores.event(event)?.allele.as_indel()?.clone();
        Some(IndelCall { allele, genotype, event })
    };
    // The called pair is always among a record's pairs, so every record has
    // likelihoods, and a called event's carry an alt.
    let score = |quality: Phred, record: &[EventIdx]| {
        let likelihoods = verdict.genotype_likelihoods(record, |other| verdict.is_called(other))?;
        Some(RecordScore { quality, likelihoods })
    };

    // Events are in record order, so a compound het's two alleles are ALT 1
    // and 2 as they stand.
    let here: SmallVec<EventIdx, MAX_EVENTS> = scores.indel_events_at(pos).collect();
    // The called pair, not each event's own best copies: each alone can say two
    // copies, and two events at one anchor cannot both be on both chromosomes.
    let compound = || {
        let [first, second] = *here.as_ref() else { return None };
        if verdict.copies_of(first) == 0 || verdict.copies_of(second) == 0 {
            return None;
        }
        // One record standing for both alleles is as sure as the less certain.
        let quality = verdict.margin(first)?.min(verdict.margin(second)?).qual();
        let score = score(quality, &[first, second])?;
        // Two indels at one anchor cannot share a chromosome, so the called
        // pair is `1/2`, and it is also the likeliest genotype of the record's
        // own likelihoods (`a_compound_het_is_always_one_two`).
        let genotype = GenotypeTag::alt_het(ONE, TWO);
        Some(IndelRecords::Compound {
            calls: [call(first, genotype)?, call(second, genotype)?],
            score,
        })
    };
    if let Some(records) = compound() {
        return records;
    }

    let separate = here.iter().filter_map(|&event| {
        let score = score(verdict.margin(event)?.qual(), &[event])?;
        let genotype = match verdict.copies_of(event) {
            0 => GenotypeTag::HomRef,
            _ => score.likelihoods.likeliest_with_an_alt()?,
        };
        Some(BiallelicRecord { call: call(event, genotype)?, score })
    });
    IndelRecords::Separate(separate.collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::pair_hmm::fixtures::*;
    use crate::call::test_helpers::pos0;
    use crate::utils::MappingQuality;
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
        call_on_pair_hmm(scores, pos0(0))
            .scored_calls()
            .map(|(c, _)| (c.allele.clone(), c.genotype))
            .collect()
    }

    /// A het, a homozygote and the reference, each as the best pair has it,
    /// with `QUAL` from the margin and `GQ` and `PL` from the record's
    /// likelihoods; under the reference the allele stays as a `0/0` call for
    /// `--all` to show.
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
            let verdict = &scores.verdict();
            for (call, score) in call_on_pair_hmm(&scores, pos0(0)).scored_calls() {
                assert_eq!(call.event, EventIdx::new(0).expect("event"));
                let margin = verdict.margin(call.event).expect("a margin");
                assert_eq!(*score.quality, *margin.qual());
                // The written genotype is the likeliest of the three, and is
                // sure on eight reads.
                let gq = score.genotype_quality(call.genotype);
                let pl = score.likelihoods.phred();
                let written = match call.genotype {
                    GenotypeTag::HomRef => 0,
                    GenotypeTag::RefHet(_) => 1,
                    _ => 2,
                };
                assert_eq!(pl.get(written), Some(&0), "{pl:?}");
                assert!(*gq > 10.0, "{gq:?}");
            }
        }
    }

    /// Two called alleles in the best pair are one compound record: both
    /// calls `1/2`, in the alleles' order, at the weaker allele's `QUAL`.
    #[test]
    fn two_called_alleles_make_a_compound_het() {
        let scores = scores(&["AA", "A"], &[(ALT1_READ, 9), (ALT2_READ, 4)]);
        let IndelRecords::Compound { calls, score } = call_on_pair_hmm(&scores, pos0(0)) else {
            panic!("a compound record");
        };
        let compound = GenotypeTag::AltHet(ALT_1, ALT_2);
        let alleles: Vec<_> = calls.iter().map(|c| (c.allele.clone(), c.genotype)).collect();
        assert_eq!(alleles, [(deletion("AA"), compound), (deletion("A"), compound)]);
        let [first, second] = [0, 1]
            .map(|i| scores.verdict().margin(EventIdx::new(i).expect("event")).expect("margin"));
        assert_eq!(*score.quality, *first.min(second).qual());
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
    /// the insertion is a het.
    #[test]
    fn a_linked_indel_takes_the_pairs_copies() {
        use crate::call::pair_hmm::{Cluster, Event, ReadScores};
        const NO: f32 = f32::NEG_INFINITY;
        let events = [
            Event { anchor: pos0(5), allele: IndelAllele::Insertion([Base::T].into()).into() },
            Event { anchor: pos0(20), allele: deletion("A").into() },
        ];
        let cluster = Cluster::new(events.into_iter().collect()).expect("a cluster");
        // By haplotype: reference, insertion, deletion, both.
        let both = ReadScores::new(
            None,
            MappingQuality::from_byte(60),
            [-9.0, -4.5, -4.5, 0.0, NO, NO, NO, NO],
        );
        let other = ReadScores::new(
            None,
            MappingQuality::from_byte(60),
            [0.0, -0.3, -4.5, -4.5, NO, NO, NO, NO],
        );
        let reads: Vec<_> = std::iter::repeat_n(both, 8).chain(std::iter::once(other)).collect();
        let scores = HaplotypeScores::new(cluster, reads.into()).expect("scores");
        let insertion = EventIdx::new(0).expect("event");

        let records = call_on_pair_hmm(&scores, pos0(5));
        let [(call, _)] = records.scored_calls().collect::<Vec<_>>()[..] else {
            panic!("one record")
        };
        assert_eq!(call.genotype, GenotypeTag::RefHet(ALT_1));
        assert_eq!(scores.verdict().copies_of(insertion), 1);
    }

    const A: u32 = 10;
    const B: u32 = 20;

    /// Two events at anchor [`A`] and one at [`B`], over `n` reads each that
    /// only the haplotype carrying the events `haps` lists explains (as bitmasks
    /// over the three events).
    fn cluster(haps: &[(usize, usize)]) -> HaplotypeScores {
        use crate::call::pair_hmm::{Cluster, Event, MAX_HAPLOTYPES, ReadScores};
        let at = |anchor, seq| Event { anchor: pos0(anchor), allele: deletion(seq).into() };
        let events = [at(A, "A"), at(A, "AA"), at(B, "C")].into_iter().collect();
        let reads: Vec<ReadScores> = haps
            .iter()
            .flat_map(|&(hap, n)| {
                let mut log10 = [-9.0; MAX_HAPLOTYPES];
                log10[hap] = -1.0;
                std::iter::repeat_n(ReadScores::new(None, MappingQuality::from_byte(60), log10), n)
            })
            .collect();
        let cluster = Cluster::new(events).expect("a cluster");
        HaplotypeScores::new(cluster, reads.into()).expect("a verdict")
    }

    fn calls_at(scores: &HaplotypeScores, pos: u32) -> (bool, Vec<(IndelAllele, GenotypeTag)>) {
        let records = call_on_pair_hmm(scores, pos0(pos));
        let calls = records.scored_calls().map(|(c, _)| (c.allele.clone(), c.genotype)).collect();
        (matches!(records, IndelRecords::Compound { .. }), calls)
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

    /// Every event at an anchor gets a record, and two the called pair carries
    /// are one compound record whose own likelihoods make `1/2` the likeliest
    /// genotype with an alt. Two indels at [`A`] cannot share a chromosome, so
    /// the called pair is `1/2` there; a genotype without one of them is at
    /// most a few pairs (the event at [`B`] placed either way) each a margin
    /// of `Margin::PRESENT` below the best pair, and when the best pair holds
    /// an uncalled event at [`B`] the called pair is less than that margin
    /// below it.
    #[hegel::test]
    fn a_compound_het_is_always_one_two(tc: hegel::TestCase) {
        use crate::call::pair_hmm::{Cluster, Event, MAX_HAPLOTYPES, ReadScores};
        use hegel::generators as gs;

        let at = |anchor, seq| Event { anchor: pos0(anchor), allele: deletion(seq).into() };
        let events = [at(A, "A"), at(A, "AA"), at(B, "C")].into_iter().collect();
        let cluster = Cluster::new(events).expect("a cluster");
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(16));
        let reads: Vec<ReadScores> = (0..n)
            .map(|_| {
                let log10: [f32; MAX_HAPLOTYPES] =
                    tc.draw(gs::arrays(gs::floats::<f32>().min_value(-12.0).max_value(0.0)));
                ReadScores::new(None, MappingQuality::from_byte(60), log10)
            })
            .collect();
        let Some(scores) = HaplotypeScores::new(cluster, reads.into()) else { return };
        let verdict = scores.verdict();

        let at_a: Vec<_> = scores.indel_events_at(pos0(A)).collect();
        let both_called = at_a.iter().all(|&event| verdict.copies_of(event) > 0);
        match call_on_pair_hmm(&scores, pos0(A)) {
            IndelRecords::Compound { calls, score } => {
                assert!(both_called, "{verdict:?}");
                let compound = GenotypeTag::AltHet(ALT_1, ALT_2);
                assert!(calls.iter().all(|call| call.genotype == compound), "{calls:?}");
                assert_eq!(score.likelihoods.likeliest_with_an_alt(), Some(compound));
            }
            IndelRecords::Separate(records) => {
                assert!(!both_called, "{verdict:?}");
                assert_eq!(records.len(), at_a.len(), "{verdict:?}");
            }
        }
        assert_eq!(call_on_pair_hmm(&scores, pos0(B)).scored_calls().count(), 1);
    }

    /// A record's `GT` is the genotype its own `PL` gives 0 among those
    /// carrying the allele, as an SNV record's is
    /// (`genotype::tests::the_genotype_is_the_likeliest_of_its_likelihoods`).
    /// Two deletions at [`A`] and [`B`], each read over the first, the second
    /// or both: a read over one scores the haplotypes with and without the
    /// other alike, so a het of either has a cis and a trans pair that tie,
    /// and the homozygote one.
    #[hegel::test]
    fn an_indel_genotype_is_the_likeliest_of_its_likelihoods(tc: hegel::TestCase) {
        use crate::call::pair_hmm::{Cluster, Event, ReadScores, fixtures::NO};
        use hegel::generators as gs;

        let events = [
            Event { anchor: pos0(A), allele: deletion("A").into() },
            Event { anchor: pos0(B), allele: deletion("C").into() },
        ];
        let cluster = Cluster::new(events.into_iter().collect()).expect("a cluster");
        let n = tc.draw(gs::integers::<usize>().min_value(1).max_value(16));
        let reads: Vec<ReadScores> = (0..n)
            .map(|_| {
                let log10 = || gs::floats::<f32>().min_value(-12.0).max_value(0.0);
                let [without, with] = [tc.draw(log10()), tc.draw(log10())];
                // By haplotype: reference, first, second, both.
                let row = match tc.draw(gs::integers::<u8>().max_value(2)) {
                    0 => [without, with, without, with],
                    1 => [without, without, with, with],
                    _ => [without, tc.draw(log10()), tc.draw(log10()), with],
                };
                let [a, b, c, d] = row;
                ReadScores::new(None, MappingQuality::from_byte(60), [a, b, c, d, NO, NO, NO, NO])
            })
            .collect();
        let Some(scores) = HaplotypeScores::new(cluster, reads.into()) else { return };

        for pos in [A, B] {
            for (call, score) in call_on_pair_hmm(&scores, pos0(pos)).scored_calls() {
                if call.genotype == GenotypeTag::HomRef {
                    continue;
                }
                let pl = score.likelihoods.phred();
                let [a, b] = call.genotype.alleles();
                let written = pl.get(GenotypeLikelihoods::index(a, b)).copied();
                let least = pl.iter().skip(1).copied().min();
                assert_eq!(written, least, "{:?} at {pos}, PL {pl:?}", call.genotype);
            }
        }
    }
}
