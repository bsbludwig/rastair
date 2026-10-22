//! The compair side: haplotypes built from indel alleles, reads with their
//! gap model, and one batched alignment per column.

use super::{
    Allele, Cluster, Event, HapIdx, HaplotypeScores, MAX_EVENTS, MAX_HAPLOTYPES, ReadScores,
};
use crate::call::{phasing::FragmentId, pileup::indels::IndelAllele};
use compair::{
    Band, Betas, ConversionModel, CpgRole, Haplotype, Log10Likelihood, PcrIndelModel, Probability,
    Read, TapsEmission, Workspace,
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
/// prior, measured below the uniform level on every row. A substitution's
/// `CpG`s are the exception, profiled per pair ([`with_profiled_cpgs`]).
const UNKNOWN_BETA: Probability = Probability::new_panicky(0.5);

/// `haplotype` with `allele` at `anchor`: a base in place of the anchor's,
/// or an indel after it. `None` when it runs past the window.
fn apply(allele: &Allele, haplotype: &Haplotype, anchor: usize) -> Option<Haplotype> {
    match allele {
        Allele::Base(base) => with_substitution(haplotype, anchor, *base),
        Allele::Indel(IndelAllele::Insertion(inserted)) => {
            haplotype.with_insertion(anchor, inserted)
        }
        Allele::Indel(IndelAllele::Deletion(deleted)) => {
            haplotype.with_deletion(anchor, deleted.len())
        }
    }
}

/// `haplotype` with `base` at `at`. compair recomputes the `CpG` roles, so a
/// substitution creates or destroys a `CpG` on this haplotype alone. Built
/// here rather than in compair: levels go on a haplotype only once all its
/// events are applied, so there are none an edit would have to carry over.
fn with_substitution(haplotype: &Haplotype, at: usize, base: Base) -> Option<Haplotype> {
    let mut bases = haplotype.bases().to_vec();
    *bases.get_mut(at)? = base;
    Some(Haplotype::new(bases))
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
    reads: Vec<HaplotypeRead>,
    scratch: &mut HaplotypeScratch,
) -> Option<HaplotypeScores> {
    prepare(reference, window_start, events, reads)?.align(scratch)
}

/// Which of a haplotype's scores a kernel haplotype gives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Levels {
    /// No divergent `CpG`: both.
    Shared,
    /// Divergent `CpG`s at the level profiled for the pair's only carrier.
    Alone,
    /// Divergent `CpG`s at the level profiled for both chromosomes.
    Doubled,
}

/// A cluster ready for the kernel: everything [`score`] needs before the
/// alignment and after it.
#[derive(Debug)]
pub struct Prepared {
    cluster: Cluster,
    /// What the kernel aligns against: each cluster haplotype once, or twice
    /// when it carries a divergent `CpG`.
    haplotypes: Vec<Haplotype>,
    /// By kernel haplotype: whose score it is, and which.
    slots: Vec<(HapIdx, Levels)>,
    reads: Vec<HaplotypeRead>,
    bands: Vec<Band>,
}

/// [`score`] up to the alignment.
pub fn prepare(
    reference: &[Base],
    window_start: u64,
    events: impl IntoIterator<Item = Event>,
    reads: Vec<HaplotypeRead>,
) -> Option<Prepared> {
    prepare_with(reference, window_start, events, reads, true)
}

/// [`prepare`], with the divergent `CpG`s' levels profiled per pair when
/// `profile` holds; without, every `CpG` is at [`UNKNOWN_BETA`], for tests
/// to compare against.
fn prepare_with(
    reference: &[Base],
    window_start: u64,
    events: impl IntoIterator<Item = Event>,
    reads: Vec<HaplotypeRead>,
    profile: bool,
) -> Option<Prepared> {
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
                    "Event runs past the haplotype window; not scored"
                );
            }
            fits
        })
        .collect();
    if kept.is_empty() {
        return None;
    }
    let mut cluster = Cluster::new(kept)?;
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
        longest = longest.max(carried.iter().map(|event| event.allele.shift()).sum());
        haplotypes.push(haplotype);
    }
    let width = band_width(longest);
    // shortcut: indel-made CpGs stay at β ½; extend once SNV-only clusters measure well.
    let substitutions_only =
        cluster.events().iter().all(|event| matches!(event.allele, Allele::Base(_)));
    let (haplotypes, slots) = if profile && substitutions_only {
        let anchors: SmallVec<usize, MAX_EVENTS> =
            cluster.events().iter().filter_map(window_anchor).collect();
        let profiled =
            with_profiled_cpgs(cluster.haplotypes(), &reference, &haplotypes, &anchors, &reads)?;
        cluster.divergent_cpgs = profiled.masks;
        (profiled.kernel, profiled.slots)
    } else {
        let slots = cluster.haplotypes().iter().map(|&hap| (hap, Levels::Shared)).collect();
        (haplotypes.into(), slots)
    };

    let bands = reads
        .iter()
        .map(|read| Band::new(width, read.offset))
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
    if bands.is_empty() {
        return None;
    }
    Some(Prepared { cluster, haplotypes, slots, reads, bands })
}

/// An all-substitution cluster's kernel haplotypes with their levels.
struct Profiled {
    /// What the kernel aligns against.
    kernel: Vec<Haplotype>,
    /// By kernel haplotype: whose score it is, and which.
    slots: Vec<(HapIdx, Levels)>,
    /// [`Cluster::divergent_cpgs`].
    masks: [u32; MAX_HAPLOTYPES],
}

/// The kernel haplotypes of an all-substitution cluster, whose haplotypes
/// (`haplotypes`, one per entry of `subsets`) share the window's coordinates:
/// each haplotype carrying a divergent `CpG` (one some haplotypes have at that
/// window index and others not) twice, its levels there profiled from the
/// reads as the pair's only carrier and as both chromosomes, so a hom-alt
/// de-novo `CpG` is not charged for being mostly methylated.
fn with_profiled_cpgs(
    subsets: &[HapIdx],
    reference: &Haplotype,
    haplotypes: &[Haplotype],
    anchors: &[usize],
    reads: &[HaplotypeRead],
) -> Option<Profiled> {
    let role_at = |haplotype: &Haplotype, q: usize| haplotype.site(q).map(|site| site.cpg);
    // A substitution changes the roles of its own base and its two neighbours.
    let mut divergent: SmallVec<(usize, CpgRole), 18> = SmallVec::new();
    for q in anchors.iter().flat_map(|&a| a.saturating_sub(1)..=a.saturating_add(1)) {
        for role in [CpgRole::TopC, CpgRole::BottomG] {
            let carriers = haplotypes.iter().filter(|h| role_at(h, q) == Some(role)).count();
            if carriers > 0 && carriers < haplotypes.len() && !divergent.contains(&(q, role)) {
                divergent.push((q, role));
            }
        }
    }

    let mut kernel = Vec::with_capacity(haplotypes.len() * 2);
    let mut slots = Vec::with_capacity(haplotypes.len() * 2);
    let mut masks = [0; MAX_HAPLOTYPES];
    for (&hap, haplotype) in subsets.iter().zip(haplotypes) {
        let mut mask = 0u32;
        let mut alone = vec![None; haplotype.len()];
        let mut doubled = vec![None; haplotype.len()];
        for (bit, &(q, role)) in divergent.iter().enumerate() {
            if role_at(haplotype, q) != Some(role) {
                continue;
            }
            mask |= 1 << bit;
            let counts = Conversions::at(reads, q, role)?;
            let partner = Partner::at(haplotypes, reference, q, role, &counts);
            if let Some(slot) = alone.get_mut(q) {
                *slot = Probability::new(counts.alone(partner)).ok();
            }
            if let Some(slot) = doubled.get_mut(q) {
                *slot = Probability::new(counts.doubled()).ok();
            }
        }
        if let Some(slot) = masks.get_mut(hap.index()) {
            *slot = mask;
        }
        if mask == 0 {
            kernel.push(haplotype.clone());
            slots.push((hap, Levels::Shared));
            continue;
        }
        for (betas, levels) in [(alone, Levels::Alone), (doubled, Levels::Doubled)] {
            let with = haplotype
                .clone()
                .with_betas(betas)
                .inspect_err(|error| {
                    tracing::warn!(
                        %error,
                        "Methylation levels do not fit their haplotype; this cluster is not \
                         scored. This is a bug in rastair, please report it."
                    );
                })
                .ok()?;
            kernel.push(with);
            slots.push((hap, levels));
        }
    }
    Some(Profiled { kernel, slots, masks })
}

/// What the partner of a divergent `CpG`'s only carrier has at its position,
/// as the converting strand reads it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Partner {
    /// The converted base (a de-novo `T>C` on OT): every unconverted read is
    /// the carrier's.
    Converted,
    /// The unconverted base, outside a `CpG`: every converted read is the
    /// carrier's.
    Unconverted,
    /// Neither: only the carrier's reads show either.
    Other,
}

impl Partner {
    /// The partner's base at window index `q` for the haplotypes with `role`
    /// there. shortcut: a pair-independent guess, the base every haplotype
    /// without the `CpG` agrees on, else the reference's.
    fn at(
        haplotypes: &[Haplotype],
        reference: &Haplotype,
        q: usize,
        role: CpgRole,
        counts: &Conversions,
    ) -> Self {
        let mut others = haplotypes
            .iter()
            .filter(|h| h.site(q).map(|site| site.cpg) != Some(role))
            .map(|h| h.bases().get(q));
        let first = others.next().flatten();
        let base = if others.all(|b| b == first) { first } else { reference.bases().get(q) };
        match base {
            Some(&b) if b == counts.converted_base => Self::Converted,
            Some(&b) if b == counts.unconverted_base => Self::Unconverted,
            _ => Self::Other,
        }
    }
}

/// A divergent `CpG`'s reads on the strand that converts it.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Conversions {
    converted: f64,
    unconverted: f64,
    converted_base: Base,
    unconverted_base: Base,
}

impl Conversions {
    /// Reads covering window index `q` on the strand that converts `role`,
    /// a read's own indels ignored. `None` for [`CpgRole::None`].
    fn at(reads: &[HaplotypeRead], q: usize, role: CpgRole) -> Option<Self> {
        let (strand, converted_base, unconverted_base) = match role {
            CpgRole::TopC => (Strand::OT, Base::T, Base::C),
            CpgRole::BottomG => (Strand::OB, Base::A, Base::G),
            CpgRole::None => return None,
        };
        let q = i64::try_from(q).ok()?;
        let mut counts =
            Self { converted: 0.0, unconverted: 0.0, converted_base, unconverted_base };
        for read in reads.iter().filter(|read| read.read.strand() == strand) {
            let base = usize::try_from(q - i64::from(read.offset))
                .ok()
                .and_then(|i| read.read.bases().get(i));
            match base {
                Some(&b) if b == converted_base => counts.converted += 1.0,
                Some(&b) if b == unconverted_base => counts.unconverted += 1.0,
                _ => {}
            }
        }
        Some(counts)
    }

    /// The level with both chromosomes carrying the `CpG`: the MAP under a
    /// Beta(2, 2) prior.
    fn doubled(self) -> f64 {
        (self.converted + 1.0) / (self.converted + self.unconverted + 2.0)
    }

    /// The level with the carrier alone: half the reads are its, and what the
    /// partner's base can produce is taken out.
    fn alone(self, partner: Partner) -> f64 {
        let n = self.converted + self.unconverted;
        let (converted, n) = match partner {
            Partner::Converted => ((n / 2.0 - self.unconverted).max(0.0), n / 2.0),
            Partner::Unconverted => (self.converted.min(n / 2.0), n / 2.0),
            Partner::Other => (self.converted, n),
        };
        (converted + 1.0) / (n + 2.0)
    }
}

/// Uniform for haplotypes without levels, which is every haplotype but those
/// Uniform for haplotypes without levels, which is every haplotype but those
/// [`with_profiled_cpgs`] built. The CPU kernel, the GPU plan and its CPU
/// rescore all take this one.
fn emission() -> TapsEmission<'static> {
    TapsEmission::new(ConversionModel::taps_default(), Betas::OfHaplotype(UNKNOWN_BETA))
}

impl Prepared {
    fn pairs(&self) -> Vec<(&Read, Band)> {
        self.reads.iter().map(|read| &read.read).zip(self.bands.iter().copied()).collect()
    }

    /// [`score`] from the alignment on, on the CPU.
    pub fn align(self, scratch: &mut HaplotypeScratch) -> Option<HaplotypeScores> {
        scratch.workspace.align_reads(
            &self.haplotypes,
            &self.pairs(),
            &emission(),
            &mut scratch.scores,
        );
        self.finish(&scratch.scores)
    }

    fn finish(self, scores: &[Log10Likelihood]) -> Option<HaplotypeScores> {
        let reads = rows(scores, &self.slots, &self.reads)?;
        HaplotypeScores::new(self.cluster, reads)
    }

    /// Adds this cluster's pairs to a launch, read-major as [`Self::align`]
    /// lays them out; the range is where its scores land.
    pub fn push(
        &self,
        pairs: &mut compair::gpu::GpuPairs,
    ) -> Result<std::ops::Range<usize>, compair::gpu::GpuError> {
        pairs.push_reads(&self.haplotypes, &self.pairs(), &emission())
    }

    /// [`Self::align`]'s result from this cluster's scores of a launch,
    /// rescoring on the CPU the pairs the GPU's `f32` cannot vouch for.
    pub fn finish_gpu(self, scores: &[compair::gpu::GpuScore]) -> Option<HaplotypeScores> {
        let emission = &emission();
        let rows = scores.chunks_exact(self.haplotypes.len().max(1));
        let likelihoods: Vec<Log10Likelihood> = rows
            .zip(self.reads.iter().zip(&self.bands))
            .flat_map(|(row, (read, &band))| {
                row.iter().zip(&self.haplotypes).map(move |(score, haplotype)| {
                    score.or_rescore(haplotype, &read.read, emission, band)
                })
            })
            .collect();
        self.finish(&likelihoods)
    }
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
/// `h`th of `slots` at `r * slots.len() + h`), without the reads outside
/// the band on every haplotype. `None` when the buffer is not one score per
/// read and haplotype: it has then lost track of which read is which.
#[expect(clippy::cast_possible_truncation, reason = "likelihoods are stored as f32")]
fn rows(
    scores: &[Log10Likelihood],
    slots: &[(HapIdx, Levels)],
    reads: &[HaplotypeRead],
) -> Option<Box<[ReadScores]>> {
    let haplotypes = slots.len();
    if haplotypes == 0
        || haplotypes > 2 * MAX_HAPLOTYPES
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
                let mut scores =
                    ReadScores::new(read.fragment, read.mapq, [f32::NEG_INFINITY; MAX_HAPLOTYPES]);
                for (&(hap, levels), l) in slots.iter().zip(row) {
                    let l = l.get() as f32;
                    let set = |log10: &mut [f32; MAX_HAPLOTYPES]| {
                        if let Some(slot) = log10.get_mut(hap.index()) {
                            *slot = l;
                        }
                    };
                    match levels {
                        Levels::Shared => {
                            set(&mut scores.log10_alone);
                            set(&mut scores.log10_doubled);
                        }
                        Levels::Alone => set(&mut scores.log10_alone),
                        Levels::Doubled => set(&mut scores.log10_doubled),
                    }
                }
                scores
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
        let events = alleles.iter().map(|allele| Event { anchor, allele: allele.clone().into() });
        score(&window, 0, events, reads.to_vec(), &mut HaplotypeScratch::default())
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
                .all(|r| r.mapq == 60 && r.log10_alone[2..].iter().all(|l| *l == f32::NEG_INFINITY))
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
            allele: IndelAllele::Deletion(deleted.iter().map(|&b| Base::from(b)).collect()).into(),
        };
        let events =
            [at(FIRST, &WINDOW[FIRST + 1..FIRST + 3]), at(SECOND, &WINDOW[SECOND + 1..SECOND + 2])];
        let scores = score(&window, 0, events, reads.clone(), &mut HaplotypeScratch::default())
            .expect("scores");
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
        assert_eq!(scores.cluster().events()[0].allele, deletion(1).into());
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

    /// A substitution is the window with one base changed, and reads carrying
    /// it on one chromosome make a het; with every read on the reference the
    /// reference pair wins.
    #[test]
    fn a_substitution_is_a_one_base_haplotype() {
        const AT: usize = 40;
        let mut alt = WINDOW.to_vec();
        alt[AT] = if WINDOW[AT] == b'A' { b'G' } else { b'A' };
        let snv = Event {
            anchor: u32::try_from(AT).expect("fits"),
            allele: Allele::Base(Base::from(alt[AT])),
        };
        let window: Vec<Base> = WINDOW.iter().map(|&b| Base::from(b)).collect();
        let reference = Haplotype::new(window.as_slice());
        let edited = apply(&snv.allele, &reference, AT).expect("inside the window");
        assert_eq!(edited.len(), reference.len());
        let changed: Vec<usize> =
            (0..edited.len()).filter(|&i| edited.bases()[i] != reference.bases()[i]).collect();
        assert_eq!(changed, [AT]);
        assert!(apply(&snv.allele, &reference, WINDOW.len()).is_none());

        let mut reads: Vec<_> = (0..6).map(|i| read_of(WINDOW, i * 4, 70, Strand::OT)).collect();
        reads.extend((0..6).map(|i| read_of(&alt, 2 + i * 4, 70, Strand::OB)));
        let scores = score(&window, 0, [snv.clone()], reads, &mut HaplotypeScratch::default())
            .expect("scores");
        let event = EventIdx::new(0).expect("event");
        assert_eq!(
            scores.verdict.called_pair(),
            HapPair::new(HapIdx::REFERENCE, HapIdx::of_event(event))
        );
        let [p_ref, p_het, p_hom] = scores.copy_posteriors(event);
        assert!(p_het > 0.99 && p_ref < 0.01 && p_hom < 0.01, "{p_ref} {p_het} {p_hom}");

        let reads: Vec<_> = (0..8).map(|i| read_of(WINDOW, i * 4, 70, Strand::OT)).collect();
        let scores =
            score(&window, 0, [snv], reads, &mut HaplotypeScratch::default()).expect("scores");
        assert_eq!(
            scores.verdict.called_pair(),
            HapPair::new(HapIdx::REFERENCE, HapIdx::REFERENCE)
        );
        assert!(!scores.verdict.margin(event).expect("scored").is_present());
    }

    /// A `C>T` that destroys a `CpG`: an OT `T` there is what methylation
    /// makes of the reference `C`, so OT reads alone barely move the margin,
    /// while the same reads on OB, where TAPS converts no `C`, call it.
    #[test]
    fn a_cpg_substitution_is_scored_by_the_emission() {
        let cg = WINDOW.windows(2).position(|w| w == b"CG").expect("a CpG");
        let mut alt = WINDOW.to_vec();
        alt[cg] = b'T';
        let snv = Event { anchor: u32::try_from(cg).expect("fits"), allele: Allele::Base(Base::T) };
        let window: Vec<Base> = WINDOW.iter().map(|&b| Base::from(b)).collect();
        let margin = |strand| {
            let reads: Vec<_> =
                (0..6).map(|i| read_of(&alt, cg.saturating_sub(10 + i * 3), 60, strand)).collect();
            let scores = score(&window, 0, [snv.clone()], reads, &mut HaplotypeScratch::default())
                .expect("scores");
            scores.verdict.margin(EventIdx::new(0).expect("event")).expect("scored").log10()
        };
        let (ot, ob) = (margin(Strand::OT), margin(Strand::OB));
        assert!(ob > 5.0, "unconverted strand: {ob}");
        assert!(ot < ob / 2.0, "a converted read is half explained by methylation: {ot} vs {ob}");
    }

    /// `bases` over the window with the base at `at` replaced by `alt`.
    fn with_base(at: usize, alt: u8) -> Vec<u8> {
        let mut edited = WINDOW.to_vec();
        edited[at] = alt;
        edited
    }

    /// A read of `haplotype` around `at` showing `shown` there, at Q30, its
    /// own fragment.
    fn read_showing(
        haplotype: &[u8],
        at: usize,
        shown: u8,
        i: usize,
        strand: Strand,
    ) -> HaplotypeRead {
        let mut bases = haplotype.to_vec();
        bases[at] = shown;
        let fragment = 1000 + u64::try_from(i).expect("fits");
        let mut read = read_of(&bases, at - 20 + i % 10, 50, strand);
        read.fragment = FragmentId::new(fragment);
        read
    }

    /// The substitution of `alt` at `at` scored over `reads`, with the
    /// divergent `CpG`s profiled or at one half.
    fn snv(at: usize, alt: u8, reads: &[HaplotypeRead], profile: bool) -> HaplotypeScores {
        let window: Vec<Base> = WINDOW.iter().map(|&b| Base::from(b)).collect();
        let event = Event {
            anchor: u32::try_from(at).expect("fits"),
            allele: Allele::Base(Base::from(alt)),
        };
        prepare_with(&window, 0, [event], reads.to_vec(), profile)
            .expect("prepared")
            .align(&mut HaplotypeScratch::default())
            .expect("scores")
    }

    /// `log10` of how much likelier the reads are under the homozygote of the
    /// pair's carrier than under the het.
    fn hom_over_het(scores: &HaplotypeScores, hom: HapIdx) -> f64 {
        let het = HapPair::new(HapIdx::REFERENCE, HapIdx::of_event(EventIdx::new(0).expect("e")));
        let under = |pair| -> f64 {
            scores.reads.iter().map(|read| read.log10_under(scores.cluster(), pair)).sum()
        };
        under(HapPair::new(hom, hom)) - under(het)
    }

    /// A `T` before a `G` in the window, not after a `C`.
    fn tg() -> usize {
        (30..WINDOW.len() - 1)
            .find(|&i| WINDOW[i] == b'T' && WINDOW[i + 1] == b'G' && WINDOW[i - 1] != b'C')
            .expect("a TG")
    }

    /// A homozygous `T>C` making a de-novo `CpG` that is 90 % methylated: on
    /// OT nine of ten reads show the converted `T`, which the reference
    /// explains for free, and a few OB reads show the `C`. At one half every
    /// converted read costs the homozygote 0.3 log10 and the het a third of
    /// that, so the het wins; the level profiled over both chromosomes makes
    /// the converted reads nearly free and the homozygote wins.
    #[test]
    fn a_methylated_hom_alt_denovo_cpg_is_homozygous() {
        let at = tg();
        let alt = with_base(at, b'C');
        let mut reads: Vec<_> = (0..10)
            .map(|i| read_showing(&alt, at, if i < 9 { b'T' } else { b'C' }, i, Strand::OT))
            .collect();
        reads.extend((10..13).map(|i| read_showing(&alt, at, b'C', i, Strand::OB)));
        let hom = HapIdx::of_event(EventIdx::new(0).expect("event"));
        let profiled = snv(at, b'C', &reads, true);
        let uniform = snv(at, b'C', &reads, false);
        assert_eq!(profiled.verdict.best, HapPair::new(hom, hom));
        assert_eq!(uniform.verdict.best, HapPair::new(HapIdx::REFERENCE, hom));
        let (p, u) = (hom_over_het(&profiled, hom), hom_over_het(&uniform, hom));
        assert!(p > 0.5 && u < 0.0, "profiled {p}, uniform {u}");

        // With ten OB reads both call the homozygote, the profile by more.
        reads.extend((13..20).map(|i| read_showing(&alt, at, b'C', i, Strand::OB)));
        let (p, u) = (
            hom_over_het(&snv(at, b'C', &reads, true), hom),
            hom_over_het(&snv(at, b'C', &reads, false), hom),
        );
        assert!(p > u + 1.0 && u > 0.0, "profiled {p}, uniform {u}");
    }

    /// An unmethylated het `T>C`: the alt's OT reads show the `C`, so the
    /// profile has nothing to explain away.
    #[test]
    fn an_unmethylated_denovo_het_stays_het() {
        let at = tg();
        let alt = with_base(at, b'C');
        let mut reads = Vec::new();
        for (i, strand) in (0..20).zip([Strand::OT, Strand::OB].into_iter().flat_map(|s| [s; 10])) {
            let (haplotype, shown) =
                if i % 2 == 0 { (WINDOW, b'T') } else { (alt.as_slice(), b'C') };
            reads.push(read_showing(haplotype, at, shown, i, strand));
        }
        let het = HapPair::new(HapIdx::REFERENCE, HapIdx::of_event(EventIdx::new(0).expect("e")));
        for profile in [true, false] {
            assert_eq!(snv(at, b'C', &reads, profile).verdict.called_pair(), het, "{profile}");
        }
    }

    /// A `C>T` destroying a reference `CpG` that is 90 % methylated, on a
    /// homozygous reference: the mirror of the de-novo case.
    #[test]
    fn a_methylated_reference_cpg_is_homozygous_reference() {
        let cg = (30..WINDOW.len() - 1)
            .find(|&i| WINDOW[i] == b'C' && WINDOW[i + 1] == b'G')
            .expect("a CpG");
        let mut reads: Vec<_> = (0..10)
            .map(|i| read_showing(WINDOW, cg, if i < 9 { b'T' } else { b'C' }, i, Strand::OT))
            .collect();
        reads.extend((10..13).map(|i| read_showing(WINDOW, cg, b'C', i, Strand::OB)));
        let reference = HapPair::new(HapIdx::REFERENCE, HapIdx::REFERENCE);
        assert_eq!(snv(cg, b'T', &reads, true).verdict.best, reference);
        let uniform = snv(cg, b'T', &reads, false).verdict.best;
        assert_ne!(uniform, reference, "at one half the converted reads make a het");
    }

    /// A de-novo `CpG` made from a `T` leaves the partner the converted base,
    /// one made from an `A` neither, and a `C` whose `G` is destroyed keeps
    /// the unconverted base.
    #[test]
    fn the_partner_is_what_the_haplotypes_without_the_cpg_read() {
        let hap = |bases: &[u8]| Haplotype::from_ascii(bases);
        let counts = |role| Conversions::at(&[], 1, role).expect("a CpG role");
        let partner = |reference: &[u8], alt: &[u8], q, role| {
            let (reference, alt) = (hap(reference), hap(alt));
            Partner::at(&[reference.clone(), alt], &reference, q, role, &counts(role))
        };
        assert_eq!(partner(b"ATGA", b"ACGA", 1, CpgRole::TopC), Partner::Converted);
        assert_eq!(partner(b"AAGA", b"ACGA", 1, CpgRole::TopC), Partner::Other);
        assert_eq!(partner(b"ACGA", b"ACTA", 1, CpgRole::TopC), Partner::Unconverted);
        assert_eq!(partner(b"ACAA", b"ACGA", 2, CpgRole::BottomG), Partner::Converted);
        assert_eq!(partner(b"ACGA", b"ATGA", 2, CpgRole::BottomG), Partner::Unconverted);
    }

    #[test]
    fn the_levels_take_out_what_the_partner_can_produce() {
        let counts = Conversions {
            converted: 9.0,
            unconverted: 1.0,
            converted_base: Base::T,
            unconverted_base: Base::C,
        };
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12;
        assert!(close(counts.doubled(), 10.0 / 12.0));
        assert!(close(counts.alone(Partner::Converted), 5.0 / 7.0));
        assert!(close(counts.alone(Partner::Unconverted), 6.0 / 7.0));
        assert!(close(counts.alone(Partner::Other), 10.0 / 12.0));
        let all_unconverted = Conversions { converted: 0.0, unconverted: 10.0, ..counts };
        assert!(close(all_unconverted.alone(Partner::Converted), 1.0 / 7.0));
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
        let slots = [HapIdx::REFERENCE, HapIdx::of_event(EventIdx::new(0).expect("event"))]
            .map(|hap| (hap, Levels::Shared));
        assert!(rows(&scores, &slots, &reads).is_none());
        assert_eq!(rows(&[scores[0]; 4], &slots, &reads).map(|r| r.len()), Some(2));
    }
}
