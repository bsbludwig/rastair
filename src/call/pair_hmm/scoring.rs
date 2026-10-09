//! The compair side: haplotypes built from indel alleles, reads with their
//! gap model, and one batched alignment per column.

use super::{Cluster, Event, HapIdx, HaplotypeScores, MAX_EVENTS, MAX_HAPLOTYPES, ReadScores};
use crate::call::{phasing::FragmentId, pileup::indels::IndelAllele};
use compair::{
    Band, Betas, ConversionModel, Haplotype, Log10Likelihood, PcrIndelModel, Probability, Read,
    TapsEmission, Workspace,
};
use seqair_types::{Base, BaseQuality, SmallVec, Strand};
use std::borrow::Cow;

/// GATK's conservative PCR indel model, but with the gap-open quality in a
/// tandem repeat floored at Q20 rather than Q10: the library's PCR stutter
/// rate in a long tract.
const PCR_INDEL_MODEL: PcrIndelModel = PcrIndelModel {
    gap_open: BaseQuality::from_byte(45),
    floor: BaseQuality::from_byte(20),
    rate: 3.0,
    max_unit: 8,
    max_repeats: 20,
};

/// Gap extension everywhere (GATK's).
const GAP_EXTEND_QUAL: BaseQuality = BaseQuality::from_byte(10);

/// What a base without a quality (`0xff` in BAM) is scored at.
const MISSING_QUAL: BaseQuality = BaseQuality::from_byte(20);

/// The methylation level every CpG is scored at: one half makes a converted
/// and an unconverted read of a CpG cost the same, so the level cancels out
/// of the genotype wherever the haplotypes share the CpG. Where they do not
/// (an indel that creates or destroys one, or a deletion spelled across a
/// methylated CpG) the column's own level would be the right one, and compair
/// can carry it on the haplotype, but on chr12 at 26x and 73x the column
/// estimate, raw or shrunk towards one half with up to sixteen reads of
/// prior, measured below the uniform level on every row.
const UNKNOWN_BETA: Probability = Probability::new_panicky(0.5);

/// `haplotype` with `allele` after the base at `anchor`; `None` when a
/// deletion runs past the window.
fn apply(allele: &IndelAllele, haplotype: &Haplotype, anchor: usize) -> Option<Haplotype> {
    match allele {
        IndelAllele::Insertion(inserted) => haplotype.with_insertion(anchor, inserted),
        IndelAllele::Deletion(deleted) => haplotype.with_deletion(anchor, deleted.len()),
    }
}

/// One read as the pair-HMM sees it, placed on the window.
#[derive(Debug, Clone, PartialEq)]
pub struct HaplotypeRead {
    read: Read,
    /// Window index the read's first base is expected at.
    offset: i32,
    fragment: Option<FragmentId>,
    mapq: u8,
}

impl HaplotypeRead {
    /// An error for a read without a bisulfite strand, which the emission
    /// cannot score, or whose bases and qualities disagree in length.
    pub fn new(
        bases: &[Base],
        quals: &[BaseQuality],
        offset: i32,
        strand: Strand,
        fragment: Option<FragmentId>,
        mapq: u8,
    ) -> Result<Self, compair::Error> {
        let quals: Cow<'_, [BaseQuality]> = if quals.contains(&BaseQuality::UNAVAILABLE) {
            quals
                .iter()
                .map(|&q| if q == BaseQuality::UNAVAILABLE { MISSING_QUAL } else { q })
                .collect()
        } else {
            Cow::Borrowed(quals)
        };
        let read =
            Read::with_pcr_indel_model(bases, &quals, &PCR_INDEL_MODEL, GAP_EXTEND_QUAL, strand)?;
        Ok(Self { read, offset, fragment, mapq })
    }
}

/// The kernel's buffers, kept for a whole region: they are sized by band width
/// and haplotype length, which barely move from one column to the next.
#[derive(Debug, Default)]
pub struct HaplotypeScratch {
    workspace: Workspace,
    scores: Vec<Log10Likelihood>,
}

/// Score `reads` against the `reference` window, which starts at genomic
/// `window_start`, and against it with every compatible subset of `events`
/// applied, each at its own genomic anchor. An event that does not fit the
/// window is left out. `None` when no event fits, the events break
/// [`Cluster::new`]'s limits, or no read can be scored.
pub fn score(
    reference: &[Base],
    window_start: u64,
    events: impl IntoIterator<Item = Event>,
    reads: &[HaplotypeRead],
    scratch: &mut HaplotypeScratch,
) -> Option<HaplotypeScores> {
    let reference = Haplotype::new(reference);
    let window_anchor =
        |event: &Event| usize::try_from(u64::from(event.anchor).checked_sub(window_start)?).ok();
    let kept: SmallVec<Event, MAX_EVENTS> = events
        .into_iter()
        .take(MAX_EVENTS)
        .filter(|event| {
            let fits = window_anchor(event)
                .is_some_and(|anchor| apply(&event.allele, &reference, anchor).is_some());
            if !fits {
                tracing::debug!(
                    ?event,
                    window_start,
                    window = reference.len(),
                    "Indel event runs past the haplotype window; not scored"
                );
            }
            fits
        })
        .collect();
    if kept.is_empty() {
        return None;
    }
    let cluster = Cluster::new(kept)?;
    // Later anchors first, so an earlier event's window index still holds
    // after the edit.
    let mut haplotypes: SmallVec<Haplotype, MAX_HAPLOTYPES> = SmallVec::new();
    let mut longest = 0usize;
    for &subset in cluster.haplotypes() {
        let mut carried: SmallVec<&Event, MAX_EVENTS> =
            subset.events().filter_map(|e| cluster.events().get(e.position())).collect();
        carried.sort_by_key(|event| std::cmp::Reverse(event.anchor));
        let mut haplotype = reference.clone();
        for event in &carried {
            haplotype = apply(&event.allele, &haplotype, window_anchor(event)?)?;
        }
        longest = longest.max(carried.iter().map(|event| event.allele.len()).sum());
        haplotypes.push(haplotype);
    }
    let width = band_width(longest);

    let pairs = reads
        .iter()
        .map(|read| Ok((&read.read, Band::new(width, read.offset)?)))
        .collect::<Result<Vec<_>, compair::Error>>()
        .inspect_err(|error| {
            tracing::warn!(
                %error,
                width,
                "The pair-HMM band was refused; this cluster's indels are not called. This is a \
                 bug in rastair, please report it."
            );
        })
        .ok()?;
    if pairs.is_empty() {
        return None;
    }

    let emission = TapsEmission::new(ConversionModel::taps_default(), Betas::Uniform(UNKNOWN_BETA));
    let haplotypes: SmallVec<&Haplotype, MAX_HAPLOTYPES> = haplotypes.iter().collect();
    scratch.workspace.align_reads(&haplotypes, &pairs, &emission, &mut scratch.scores);
    let reads = rows(&scratch.scores, cluster.haplotypes(), reads)?;
    HaplotypeScores::new(cluster, reads)
}

/// A band wide enough for a read carrying the longest indel, which sits that
/// many diagonals off the reference. One width for every read of a column, so
/// the kernel can pack their pairs together.
fn band_width(longest_indel: usize) -> u32 {
    let width = longest_indel
        .checked_mul(2)
        .and_then(|w| w.checked_add(32))
        .and_then(|w| u32::try_from(w).ok())
        .unwrap_or(Band::MAX_WIDTH);
    width.max(Band::DEFAULT_WIDTH).next_multiple_of(8).min(Band::MAX_WIDTH)
}

/// One row per read from compair's read-major buffer (read `r` against the
/// `h`th of `subsets` at `r * subsets.len() + h`), without the reads outside
/// the band on every haplotype. `None` when the buffer is not one score per
/// read and haplotype: it has then lost track of which read is which.
#[expect(clippy::cast_possible_truncation, reason = "likelihoods are stored as f32")]
fn rows(
    scores: &[Log10Likelihood],
    subsets: &[HapIdx],
    reads: &[HaplotypeRead],
) -> Option<Box<[ReadScores]>> {
    let haplotypes = subsets.len();
    if haplotypes == 0
        || haplotypes > MAX_HAPLOTYPES
        || Some(scores.len()) != reads.len().checked_mul(haplotypes)
    {
        tracing::warn!(
            reads = reads.len(),
            haplotypes,
            scores = scores.len(),
            "The pair-HMM returned a score buffer of the wrong length; this cluster's indels \
             are not called. This is a bug in rastair, please report it."
        );
        return None;
    }
    Some(
        scores
            .chunks_exact(haplotypes)
            .zip(reads)
            .filter(|(row, _)| row.iter().any(|l| l.get().is_finite()))
            .map(|(row, read)| {
                let mut log10 = [f32::NEG_INFINITY; MAX_HAPLOTYPES];
                for (subset, l) in subsets.iter().zip(row) {
                    if let Some(slot) = log10.get_mut(subset.index()) {
                        *slot = l.get() as f32;
                    }
                }
                ReadScores { fragment: read.fragment, mapq: read.mapq, log10 }
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::pair_hmm::{EventIdx, HapIdx, HapPair};

    const WINDOW: &[u8] = b"ACGTTAGCATCGGATCCATGACTGATCGGCTAGCTTACGGATCCTAGGCTAAGCTTTTTTTTTTTTTGACTAGCATCGGATCGATCAGCTAGCTGATCGATGCTAGCTAGTC";

    /// The base before the window's thirteen-base T run.
    fn anchor() -> usize {
        WINDOW.windows(13).position(|w| w == b"TTTTTTTTTTTTT").expect("the run") - 1
    }

    fn deletion(n: usize) -> IndelAllele {
        IndelAllele::Deletion(std::iter::repeat_n(Base::T, n).collect())
    }

    fn insertion(bases: &str) -> IndelAllele {
        IndelAllele::Insertion(bases.bytes().map(Base::from).collect())
    }

    /// A read copied off `haplotype` from `start`, `len` long, at Q30, its own
    /// fragment.
    fn read_of(haplotype: &[u8], start: usize, len: usize, strand: Strand) -> HaplotypeRead {
        let slice = &haplotype[start..(start + len).min(haplotype.len())];
        let bases: Vec<Base> = slice.iter().map(|&b| Base::from(b)).collect();
        let quals = vec![BaseQuality::from_byte(30); slice.len()];
        let offset = i32::try_from(start).expect("fits");
        let fragment = FragmentId::new(u64::try_from(start).expect("fits") + 1);
        HaplotypeRead::new(&bases, &quals, offset, strand, fragment, 60).expect("a read")
    }

    /// The window with `n` bases of the T run deleted.
    fn alt_window(n: usize) -> Vec<u8> {
        let a = anchor();
        let mut w = WINDOW[..=a].to_vec();
        w.extend_from_slice(&WINDOW[a + 1 + n..]);
        w
    }

    fn indels(alleles: &[IndelAllele], reads: &[HaplotypeRead]) -> Option<HaplotypeScores> {
        let window: Vec<Base> = WINDOW.iter().map(|&b| Base::from(b)).collect();
        let anchor = u32::try_from(anchor()).expect("fits");
        let events = alleles.iter().map(|allele| Event { anchor, allele: allele.clone() });
        score(&window, 0, events, reads, &mut HaplotypeScratch::default())
    }

    fn best(scores: &HaplotypeScores) -> (usize, usize) {
        let (low, high) = scores.verdict.best.haplotypes();
        (low.index(), high.index())
    }

    #[test]
    fn a_missing_quality_is_scored_at_q20() {
        let read = HaplotypeRead::new(
            &[Base::A, Base::C],
            &[BaseQuality::UNAVAILABLE, BaseQuality::from_byte(30)],
            0,
            Strand::OT,
            None,
            60,
        )
        .expect("a read");
        assert_eq!(read.read.base_quals(), [MISSING_QUAL, BaseQuality::from_byte(30)]);
        assert!(HaplotypeRead::new(&[Base::A], &[], 0, Strand::OT, None, 60).is_err());
    }

    /// A reference-length stutter read at a homozygous deletion in a long run
    /// costs a gap at the stutter floor, which is less than the half-mixture
    /// penalty every read pays under a het.
    #[test]
    fn a_stutter_read_does_not_flip_a_homozygote() {
        let alt = alt_window(1);
        let mut reads: Vec<_> = (0..11).map(|i| read_of(&alt, i * 3, 70, Strand::OT)).collect();
        reads.push(read_of(WINDOW, anchor() - 40, 80, Strand::OB));
        let scores = indels(&[deletion(1)], &reads).expect("scores");
        assert_eq!(best(&scores), (1, 1));
    }

    /// Every row names its read's fragment and mapping quality, and has no
    /// score for the haplotype the column does not have.
    #[test]
    fn reference_reads_are_homozygous_reference() {
        let reads: Vec<_> = (0..8).map(|i| read_of(WINDOW, i * 3, 60, Strand::OT)).collect();
        let scores = indels(&[deletion(1)], &reads).expect("scores");
        assert_eq!(best(&scores), (0, 0));
        assert_eq!(scores.reads.len(), 8);
        let fragments: Vec<_> = scores.reads.iter().map(|r| r.fragment).collect();
        let expected: Vec<_> = (0..8u64).map(|i| FragmentId::new(i * 3 + 1)).collect();
        assert_eq!(fragments, expected);
        assert!(
            scores
                .reads
                .iter()
                .all(|r| r.mapq == 60 && r.log10[2..].iter().all(|l| *l == f32::NEG_INFINITY))
        );
        let alt = HapIdx::of_event(EventIdx::new(0).expect("the one scored event"));
        let preference: f32 =
            scores.reads.iter().filter_map(|r| r.log10_ratio(HapIdx::REFERENCE, alt)).sum();
        assert!(preference > 0.0, "the reads prefer the reference: {preference}");
    }

    #[test]
    fn a_het_is_a_het() {
        let alt = alt_window(1);
        let mut reads: Vec<_> = (0..6).map(|i| read_of(WINDOW, i * 4, 70, Strand::OT)).collect();
        reads.extend((0..6).map(|i| read_of(&alt, i * 4, 70, Strand::OB)));
        let scores = indels(&[deletion(1)], &reads).expect("scores");
        assert_eq!(
            scores.verdict.called_pair(),
            HapPair::new(HapIdx::REFERENCE, HapIdx::of_event(EventIdx::new(0).expect("event")))
        );
    }

    /// Reads that end inside the tract carry no indel and the column counts
    /// them for the reference; here 4 of 12 do, and the pair-HMM scores them
    /// alike under both haplotypes so they drop out.
    #[test]
    fn reads_ending_in_the_tract_do_not_vote_for_the_reference() {
        let alt = alt_window(1);
        let a = anchor();
        let mut reads: Vec<_> = (0..8).map(|i| read_of(&alt, i * 3, 70, Strand::OT)).collect();
        reads.extend((0..4).map(|i| read_of(WINDOW, i * 2, a + 6 - i * 2, Strand::OB)));
        let scores = indels(&[deletion(1)], &reads).expect("scores");
        assert_eq!(best(&scores), (1, 1));
    }

    #[test]
    fn two_alleles_can_make_a_compound_het() {
        let (one, two) = (alt_window(1), alt_window(2));
        let mut reads: Vec<_> = (0..6).map(|i| read_of(&one, i * 4, 70, Strand::OT)).collect();
        reads.extend((0..6).map(|i| read_of(&two, i * 4, 70, Strand::OB)));
        let scores = indels(&[deletion(1), deletion(2)], &reads).expect("scores");
        assert_eq!(best(&scores), (1, 2));
        let [first, second] = [0, 1].map(|i| EventIdx::new(i).expect("event"));
        assert!(
            scores
                .cluster()
                .haplotypes()
                .iter()
                .all(|hap| !(hap.carries(first) && hap.carries(second))),
            "two events at one anchor are never one haplotype: {:?}",
            scores.cluster().haplotypes()
        );
    }

    /// Two deletions 13 bp apart that every read carries together, as at
    /// chr12:50,733,241 and 50,733,254: the haplotype carrying both explains
    /// the reads, and each event is homozygous. With one event per haplotype
    /// each would lose the reads to the other's mismatches.
    #[test]
    fn two_deletions_on_one_chromosome_are_one_haplotype() {
        const FIRST: usize = 30;
        const SECOND: usize = FIRST + 13;
        let mut both = WINDOW[..=FIRST].to_vec();
        both.extend_from_slice(&WINDOW[FIRST + 3..=SECOND]);
        both.extend_from_slice(&WINDOW[SECOND + 2..]);
        let reads: Vec<_> = (0..10).map(|i| read_of(&both, i * 3, 70, Strand::OT)).collect();
        let window: Vec<Base> = WINDOW.iter().map(|&b| Base::from(b)).collect();
        let at = |anchor: usize, deleted: &[u8]| Event {
            anchor: u32::try_from(anchor).expect("fits"),
            allele: IndelAllele::Deletion(deleted.iter().map(|&b| Base::from(b)).collect()),
        };
        let events =
            [at(FIRST, &WINDOW[FIRST + 1..FIRST + 3]), at(SECOND, &WINDOW[SECOND + 1..SECOND + 2])];
        let scores =
            score(&window, 0, events, &reads, &mut HaplotypeScratch::default()).expect("scores");
        let [first, second] = [0, 1].map(|i| EventIdx::new(i).expect("event"));
        let both = HapIdx::of_event(first).index() | HapIdx::of_event(second).index();
        assert_eq!(best(&scores), (both, both));
        assert_eq!(scores.verdict.called_pair(), scores.verdict.best);
        assert_eq!(scores.verdict.copies_of(first), 2);
        assert_eq!(scores.verdict.copies_of(second), 2);
    }

    #[test]
    fn a_deletion_past_the_window_is_not_a_haplotype() {
        let reads = vec![read_of(WINDOW, 0, 60, Strand::OT)];
        assert!(indels(&[deletion(500)], &reads).is_none());
        let scores = indels(&[deletion(500), deletion(1)], &reads).expect("the one that fits");
        assert_eq!(scores.cluster().events().len(), 1);
        assert_eq!(scores.cluster().events()[0].allele, deletion(1));
    }

    /// Two insertions equal under a C->T collapse are two haplotypes, and the
    /// emission separates them: reads carrying `TTAC` on the bottom strand,
    /// where no conversion turns a C into a T, pick the unconverted allele.
    #[test]
    fn the_emission_chooses_between_two_conversion_twins() {
        let a = anchor();
        let mut window = WINDOW[..=a].to_vec();
        window.extend_from_slice(b"TTAC");
        window.extend_from_slice(&WINDOW[a + 1..]);
        let reads: Vec<_> = (0..10).map(|i| read_of(&window, i * 3, 70, Strand::OB)).collect();
        let scores = indels(&[insertion("TCAC"), insertion("TTAC")], &reads).expect("scores");
        assert_eq!(best(&scores), (2, 2), "the unconverted twin");
    }

    /// The emission has no strand to score such a read on.
    #[test]
    fn a_read_without_a_strand_is_refused() {
        let read = HaplotypeRead::new(
            &[Base::A],
            &[BaseQuality::from_byte(30)],
            0,
            Strand::Unknown,
            None,
            60,
        );
        assert_eq!(read, Err(compair::Error::UnknownStrand));
    }

    /// A score buffer that is not one score per read and haplotype has lost
    /// track of which read is which, so it gives no scores.
    #[test]
    fn a_short_score_buffer_gives_no_scores() {
        let scores = [Log10Likelihood::new(-1.0); 3];
        let reads = [read_of(WINDOW, 0, 60, Strand::OT), read_of(WINDOW, 3, 60, Strand::OT)];
        let subsets = [HapIdx::REFERENCE, HapIdx::of_event(EventIdx::new(0).expect("event"))];
        assert!(rows(&scores, &subsets, &reads).is_none());
        assert_eq!(rows(&[scores[0]; 4], &subsets, &reads).map(|r| r.len()), Some(2));
    }
}
