#[cfg(test)]
use super::likelihood::log_sum;
use super::likelihood::{BaseMeasurement, MolecularBase, consensus, log_add};
use super::{Fragment, Haplotype, LikelihoodParams, Reference, Site};
use crate::metrics::methylation::CpgSide;
use color_eyre::eyre::{OptionExt as _, Result, ensure};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use seqair_types::{Base, BaseQuality, Probability, Strand};
use std::collections::{BTreeMap, BTreeSet};

/// An unanchored replacement starting at a zero-based reference coordinate.
/// Exactly one of `deleted` and `inserted` must be nonempty.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Indel {
    pub start: u32,
    pub deleted: Vec<Base>,
    pub inserted: Vec<Base>,
}

impl Ord for Indel {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.start
            .cmp(&other.start)
            .then_with(|| {
                self.deleted
                    .iter()
                    .map(|base| *base as u8)
                    .cmp(other.deleted.iter().map(|base| *base as u8))
            })
            .then_with(|| {
                self.inserted
                    .iter()
                    .map(|base| *base as u8)
                    .cmp(other.inserted.iter().map(|base| *base as u8))
            })
    }
}

impl PartialOrd for Indel {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Indel {
    pub fn normalized(&self, reference: &Reference) -> Result<Self> {
        ensure!(
            self.deleted.is_empty() != self.inserted.is_empty(),
            "Expected a simple insertion or deletion"
        );
        ensure!(
            self.deleted.iter().chain(&self.inserted).all(|base| base.known_index().is_some()),
            "Unknown indel base"
        );
        ensure!(reference.base_at(self.start).is_some(), "Indel outside reference");
        for (offset, &base) in self.deleted.iter().enumerate() {
            let pos = self
                .start
                .checked_add(u32::try_from(offset)?)
                .ok_or_eyre("Indel coordinate overflow")?;
            ensure!(reference.base_at(pos) == Some(base), "Deletion reference mismatch");
        }
        let mut normalized = self.clone();
        let sequence = if normalized.deleted.is_empty() {
            &mut normalized.inserted
        } else {
            &mut normalized.deleted
        };
        while let Some(previous) = normalized.start.checked_sub(1) {
            if reference.base_at(previous) != sequence.last().copied() {
                break;
            }
            sequence.rotate_right(1);
            normalized.start = previous;
        }
        Ok(normalized)
    }

    pub(super) fn len(&self) -> usize {
        self.deleted.len().max(self.inserted.len())
    }

    pub(super) fn equivalent_end(&self, reference: &Reference) -> Result<u32> {
        let mut start = self.start;
        let mut sequence =
            if self.deleted.is_empty() { self.inserted.clone() } else { self.deleted.clone() };
        let deleted = u32::try_from(self.deleted.len())?;
        while let Some(next) = start.checked_add(deleted) {
            if reference.base_at(next) != sequence.first().copied() {
                break;
            }
            start = start.checked_add(1).ok_or_eyre("Indel equivalent placement overflow")?;
            sequence.rotate_left(1);
        }
        start.checked_add(deleted.max(1)).ok_or_eyre("Indel equivalent span overflow")
    }

    pub(super) fn deletes(&self, pos: u32) -> bool {
        pos.checked_sub(self.start)
            .is_some_and(|offset| u64::from(offset) < self.deleted.len() as u64)
    }
}

/// Inserted bases sort before the reference base at `pos`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SequencePosition {
    pub pos: u32,
    offset: SequenceOffset,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum SequenceOffset {
    Insertion(u32),
    Reference,
}

impl SequencePosition {
    pub fn reference(pos: u32) -> Self {
        Self { pos, offset: SequenceOffset::Reference }
    }

    pub fn insertion(pos: u32, offset: u32) -> Self {
        Self { pos, offset: SequenceOffset::Insertion(offset) }
    }

    pub fn is_reference(self) -> bool {
        self.offset == SequenceOffset::Reference
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SequenceObservation {
    pub position: SequencePosition,
    pub base: Base,
    pub quality: Option<u8>,
}

#[derive(Debug, Clone)]
pub struct SequenceEvidence {
    pub bases: Vec<SequenceObservation>,
    pub reference_skip: bool,
    pub indels: Vec<Indel>,
}

#[derive(Debug, Clone, Copy)]
struct TemplateBase {
    position: SequencePosition,
    base: Base,
    cpg: bool,
    cpg_key: Option<super::methylation::Cpg>,
}

pub(super) struct Template {
    bases: Vec<TemplateBase>,
    /// Each read-local window as the range of `bases` it holds; the bases are
    /// in position order, so a window is one contiguous run of them.
    local: BTreeMap<(u32, u32), std::ops::Range<usize>>,
}

impl Template {
    fn window(&self, start: u32, end: u32) -> std::ops::Range<usize> {
        let index = |pos: u32| self.bases.partition_point(|base| base.position.pos < pos);
        index(start)..index(end)
    }

    fn local_bases(&self, window: (u32, u32)) -> Option<&[TemplateBase]> {
        self.bases.get(self.local.get(&window)?.clone())
    }
}

/// The haplotype's alleles inside `start..end`, which is all a lookup in that
/// window can reach; the beam's haplotypes carry every site of the block.
fn windowed(haplotype: &Haplotype, start: u32, end: u32) -> Haplotype {
    Haplotype {
        alleles: haplotype
            .alleles
            .iter()
            .filter(|(pos, _)| (start..end).contains(pos))
            .copied()
            .collect(),
        indels: haplotype.indels.clone(),
    }
}

fn sequence(
    reference: &Reference,
    haplotype: &Haplotype,
    start: u32,
    end: u32,
) -> Result<Vec<TemplateBase>> {
    let haplotype = windowed(haplotype, start, end);
    let haplotype = &haplotype;
    let mut bases = Vec::with_capacity(usize::try_from(end.saturating_sub(start))?);
    'positions: for pos in start..end {
        for indel in &haplotype.indels {
            if pos == indel.start {
                for (offset, &base) in indel.inserted.iter().enumerate() {
                    bases.push(TemplateBase {
                        position: SequencePosition::insertion(pos, u32::try_from(offset)?),
                        base,
                        cpg: false,
                        cpg_key: None,
                    });
                }
            }
            if indel.deletes(pos) {
                continue 'positions;
            }
        }
        bases.push(TemplateBase {
            position: SequencePosition::reference(pos),
            base: haplotype
                .base_at(reference, pos)
                .ok_or_eyre("Realignment template outside reference")?,
            cpg: false,
            cpg_key: None,
        });
    }
    Ok(bases)
}

pub(super) fn cpgs(
    reference: &Reference,
    haplotype: &Haplotype,
    start: u32,
    end: u32,
) -> Result<BTreeSet<super::methylation::Cpg>> {
    let haplotype = windowed(haplotype, start, end);
    let mut cpgs = BTreeSet::new();
    let mut previous: Option<(SequencePosition, Base)> = None;
    let mut visit = |position, base| {
        if let Some((c, Base::C)) = previous
            && base == Base::G
        {
            cpgs.insert(super::methylation::Cpg { c, g: position });
        }
        previous = Some((position, base));
    };
    'positions: for pos in start..end {
        for indel in &haplotype.indels {
            if pos == indel.start {
                for (offset, &base) in indel.inserted.iter().enumerate() {
                    visit(SequencePosition::insertion(pos, u32::try_from(offset)?), base);
                }
            }
            if indel.deletes(pos) {
                continue 'positions;
            }
        }
        visit(
            SequencePosition::reference(pos),
            haplotype
                .base_at(reference, pos)
                .ok_or_eyre("Realignment template outside reference")?,
        );
    }
    Ok(cpgs)
}

type Profile = Vec<MolecularBase>;

/// A profile as compair scores it under every template: its bases and
/// qualities, with gap qualities from the reference's repeat context at the
/// positions they align to (as GATK takes them from the read's own context),
/// so the forward results under one template serve every other through the
/// memo. `key` is what determines `read`: the memo's key.
struct ProfileTerms {
    read: compair::Read,
    key: Box<[u8]>,
}

/// The gap qualities of every reference position of a window, from the
/// library's repeat profile; a position outside it has the plain qualities.
struct RepeatTracks {
    start: u32,
    gaps: Vec<[BaseQuality; 3]>,
    plain: [BaseQuality; 3],
}

impl RepeatTracks {
    fn new(
        reference: &Reference,
        window: std::ops::Range<u32>,
        profile: &super::ErrorProfile,
        gaps: GapQualities,
    ) -> Result<Self> {
        let start = window.start.max(reference.start);
        let end = window
            .end
            .max(start)
            .min(reference.start.saturating_add(u32::try_from(reference.bases.len())?));
        let bases = reference
            .bases
            .get(
                usize::try_from(start.saturating_sub(reference.start))?
                    ..usize::try_from(end.saturating_sub(reference.start))?,
            )
            .ok_or_eyre("Repeat track window outside the reference")?;
        Ok(Self {
            start,
            gaps: profile.tracks(bases, gaps.open, gaps.extend),
            plain: [gaps.open, gaps.open, gaps.extend],
        })
    }

    /// An inserted base takes the qualities of the reference base it follows.
    fn at(&self, position: SequencePosition) -> [BaseQuality; 3] {
        position
            .pos
            .checked_sub(self.start)
            .and_then(|offset| self.gaps.get(usize::try_from(offset).ok()?))
            .copied()
            .unwrap_or(self.plain)
    }
}

/// Where a joint profile sits: the band of every forward pass over it is
/// anchored on its first base and wide enough for its own indels.
#[derive(Debug, Clone, Copy)]
struct ProfileSpan {
    first: SequencePosition,
    /// Inserted bases plus deleted reference bases along the profile: how far
    /// its own alignment can leave the diagonal it starts on.
    indel_bases: usize,
}

struct MolecularEvidence {
    individual: Vec<Profile>,
    joint: Vec<Profile>,
    intervals: Vec<(u32, u32)>,
    spans: Vec<ProfileSpan>,
    rows: Vec<ProfileTerms>,
}

/// The gap model as compair's whole-Phred read qualities.
#[derive(Debug, Clone, Copy)]
struct GapQualities {
    open: BaseQuality,
    extend: BaseQuality,
}

impl GapQualities {
    fn new(params: &LikelihoodParams) -> Result<Self> {
        Ok(Self {
            open: phred(*params.gap_open, "gap open")?,
            extend: phred(*params.gap_extend, "gap extension")?,
        })
    }
}

/// `probability` as a Phred quality, which compair takes in whole steps; an
/// error rather than a rounded model where it is not one.
fn phred(probability: f64, name: &str) -> Result<BaseQuality> {
    let quality = (-10.0 * probability.log10()).round();
    ensure!(
        (0.0..=254.0).contains(&quality)
            && (10.0_f64.powf(-quality / 10.0) - probability).abs() <= 1e-9 * probability,
        "PairHMM {name} probability {probability} is not a whole Phred quality"
    );
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a whole number in 0..=254, checked above"
    )]
    Ok(BaseQuality::from_byte(quality as u8))
}

/// A read base no model can predict: high enough that compair's match
/// emission for an unknown base, `1 - eps`, is one in `f64`.
const UNINFORMATIVE: BaseQuality = BaseQuality::from_byte(254);

/// Beyond a pair's own indels, how far the band reaches either side of the
/// diagonal its first base starts on: room for an indel placed elsewhere in a
/// repeat, and for a first base that is not where its anchor says.
const BAND_SLACK: usize = 16;

/// The band for a pass whose first base aligns at template index `offset`
/// and whose path can leave that diagonal by up to `indel_bases`: at least
/// compair's default, and in its `16n - 2` widths, which fill every vector.
fn band(offset: usize, indel_bases: usize) -> Result<compair::Band> {
    let half_width =
        (compair::Band::DEFAULT_WIDTH as usize / 2).max(indel_bases.saturating_add(BAND_SLACK));
    let width = (2 * half_width + 2).div_ceil(16) * 16 - 2;
    let width = u32::try_from(width)?.min(compair::Band::MAX_WIDTH);
    Ok(compair::Band::new(width, i32::try_from(offset)?)?)
}

/// Inserted plus deleted bases along ascending sequence positions.
fn indel_bases(positions: impl IntoIterator<Item = SequencePosition>) -> usize {
    let mut previous: Option<SequencePosition> = None;
    let mut total = 0usize;
    for position in positions {
        if !position.is_reference() {
            total += 1;
        } else if let Some(previous) = previous.filter(|previous| previous.is_reference()) {
            total += position.pos.saturating_sub(previous.pos).saturating_sub(1) as usize;
        }
        if position.is_reference() {
            previous = Some(position);
        }
    }
    total
}

/// A template under one read orientation: its bases, with each column's
/// conversion rate carried as the haplotype's level there (see
/// [`TemplateEmission`]).
///
/// Most passes over a beam step find their forward results in the memo, so a
/// track is only its key: the compair haplotype is decoded from it for the
/// passes that have to run.
struct TemplateTrack {
    /// Each column's base, then its level's bits or `u64::MAX` for none: the
    /// memo's key, and the whole of the track.
    key: Box<[u64]>,
    /// Inserted plus deleted bases relative to the reference.
    indel_bases: usize,
}

impl TemplateTrack {
    fn len(&self) -> usize {
        self.key.len() / 2
    }

    fn haplotype(&self) -> Result<compair::Haplotype> {
        let mut bases = Vec::with_capacity(self.len());
        let mut levels = Vec::with_capacity(self.len());
        let (columns, rest) = self.key.as_chunks::<2>();
        ensure!(rest.is_empty(), "Malformed template key");
        for &[base, level] in columns {
            bases.push(Base::from(u8::try_from(base)?));
            levels.push(
                (level != u64::MAX).then(|| Probability::new(f64::from_bits(level))).transpose()?,
            );
        }
        Ok(compair::Haplotype::new(bases).with_betas(levels)?)
    }
}

/// The emission of every template column, read from the column itself.
///
/// A column with a level converts at that rate between the conversion pair
/// its base belongs to, `C`/`T` or `G`/`A`: an unconverted base and its
/// conversion product weigh `1 - rate` and `rate`, and a neutral column
/// weighs both a half whichever of the two the template holds. A template is
/// built for one read orientation, so only the cytosines that orientation can
/// convert carry a level, and the read's own strand is not needed.
struct TemplateEmission;

impl compair::Emission for TemplateEmission {
    fn site_weights(&self, site: compair::HapSite, _strand: Strand) -> compair::SiteWeights {
        match (site.beta, site.base) {
            (Some(rate), Base::C | Base::T) => {
                compair::SiteWeights { base: Base::C, converted: Some(Base::T), rate: *rate }
            }
            (Some(rate), Base::G | Base::A) => {
                compair::SiteWeights { base: Base::G, converted: Some(Base::A), rate: *rate }
            }
            _ => compair::SiteWeights::plain(site.base),
        }
    }

    fn epsilon(&self, observation: compair::Observation) -> f64 {
        observation.error_probability
    }
}

pub(super) struct Prepared {
    start: u32,
    end: u32,
    padding: u32,
    reference_bounds: std::ops::Range<u32>,
    fragments: Vec<Option<MolecularEvidence>>,
    localized: bool,
    /// The total indel length on the most edited candidate chromosome. A read can carry
    /// edits that its CIGAR does not show, so every band reaches this far,
    /// whatever the template it scores against holds.
    indel_reach: usize,
}

impl Prepared {
    pub(super) fn interval(&self) -> std::ops::Range<u32> {
        self.start..self.end
    }

    pub(super) fn new(
        reference: &Reference,
        sites: &[Site],
        indel: Option<&Indel>,
        fragments: &[&Fragment],
        params: &LikelihoodParams,
    ) -> Result<Self> {
        let mut start =
            sites.iter().map(|site| site.pos).min().ok_or_eyre("No realignment sites")?;
        let mut end = sites
            .iter()
            .map(|site| site.pos)
            .max()
            .ok_or_eyre("No realignment sites")?
            .checked_add(1)
            .ok_or_eyre("Site coordinate overflow")?;
        if let Some(indel) = indel {
            start = start.min(indel.start);
            end = end.max(
                indel
                    .start
                    .checked_add(u32::try_from(indel.deleted.len().max(1))?)
                    .ok_or_eyre("Indel coordinate overflow")?,
            );
        }
        start = start.saturating_sub(u32::from(params.realignment_flank)).max(reference.start);
        let reference_end = reference
            .start
            .checked_add(u32::try_from(reference.bases.len())?)
            .ok_or_eyre("Reference coordinate overflow")?;
        end = end.saturating_add(u32::from(params.realignment_flank)).min(reference_end);
        let gaps = GapQualities::new(params)?;
        let mut extent: Option<(u32, u32)> = None;
        for read in fragments.iter().flat_map(|fragment| &fragment.reads) {
            if let Some((first, last)) = read
                .sequence
                .as_ref()
                .and_then(|sequence| sequence.bases.first().zip(sequence.bases.last()))
            {
                let (first, last) = (first.position.pos, last.position.pos.saturating_add(1));
                extent = Some(extent.map_or((first, last), |(a, b)| (a.min(first), b.max(last))));
            }
        }
        let tracks = RepeatTracks::new(
            reference,
            extent.map_or(start..end, |(first, last)| first..last),
            &params.error_profile,
            gaps,
        )?;
        let prepared: Vec<_> = fragments
            .iter()
            .map(|fragment| {
                molecular_evidence(fragment, reference.start, reference_end, gaps, &tracks)
            })
            .collect::<Result<_>>()?;
        let flank =
            u32::from(params.realignment_flank) + u32::try_from(indel.map_or(0, Indel::len))?;
        for molecule in prepared.iter().flatten() {
            for &(first, last) in &molecule.intervals {
                start = start.min(first.saturating_sub(flank).max(reference.start));
                end = end.max(last.saturating_add(flank).min(reference_end));
            }
        }
        Ok(Self {
            start,
            end,
            padding: flank,
            reference_bounds: reference.start..reference_end,
            fragments: prepared,
            localized: false,
            indel_reach: indel.map_or(0, Indel::len),
        })
    }

    pub(super) fn reaching(mut self, haplotypes: &[Haplotype]) -> Result<Self> {
        let reach = haplotypes
            .iter()
            .map(|h| h.indels.iter().map(Indel::len).sum::<usize>())
            .fold(self.indel_reach, usize::max);
        let extra = u32::try_from(reach - self.indel_reach)?;
        self.start = self.start.saturating_sub(extra).max(self.reference_bounds.start);
        self.end = self.end.saturating_add(extra).min(self.reference_bounds.end);
        self.padding = self.padding.checked_add(extra).ok_or_eyre("Template padding overflow")?;
        self.indel_reach = reach;
        Ok(self)
    }

    pub(super) fn localize(mut self) -> Self {
        self.localized = true;
        self
    }

    #[cfg(any(feature = "experimental-seqair", test))]
    pub(super) fn signature(
        &self,
        haplotype: &Haplotype,
        params: &LikelihoodParams,
    ) -> Vec<(u32, u8)> {
        self.signature_for(haplotype, params, None)
    }

    pub(super) fn signature_at(
        &self,
        index: usize,
        haplotype: &Haplotype,
        params: &LikelihoodParams,
    ) -> Vec<(u32, u8)> {
        self.signature_for(haplotype, params, Some(index))
    }

    fn signature_for(
        &self,
        haplotype: &Haplotype,
        _params: &LikelihoodParams,
        index: Option<usize>,
    ) -> Vec<(u32, u8)> {
        let intervals: Vec<_> = self
            .fragments
            .iter()
            .enumerate()
            .filter(|(offset, _)| index.is_none_or(|index| *offset == index))
            .flat_map(|(_, fragment)| fragment)
            .flat_map(|molecule| &molecule.intervals)
            .copied()
            .collect();
        haplotype
            .alleles
            .iter()
            .filter(|(pos, _)| {
                if self.localized && !intervals.is_empty() {
                    intervals.iter().any(|&(start, end)| {
                        start.saturating_sub(self.padding + 1).max(self.start.saturating_sub(1))
                            <= *pos
                            && *pos <= end.saturating_add(self.padding).min(self.end)
                    })
                } else {
                    self.start.saturating_sub(1) <= *pos && *pos <= self.end
                }
            })
            .map(|&(pos, base)| (pos, base as u8))
            .collect()
    }

    pub(super) fn realigned_fragments(&self) -> usize {
        self.fragments.iter().filter(|fragment| fragment.is_some()).count()
    }

    pub(super) fn fallback_fragments(&self) -> usize {
        self.fragments.len() - self.realigned_fragments()
    }

    pub(super) fn template(
        &self,
        reference: &Reference,
        haplotype: &Haplotype,
        _loci: &[u32],
        _params: &LikelihoodParams,
    ) -> Result<Template> {
        let reference_end = reference
            .start
            .checked_add(u32::try_from(reference.bases.len())?)
            .ok_or_eyre("Reference coordinate overflow")?;
        let mut bases = sequence(
            reference,
            haplotype,
            self.start.saturating_sub(1).max(reference.start),
            self.end.saturating_add(1).min(reference_end),
        )?;
        let mut annotations = Vec::new();
        for (index, pair) in bases.windows(2).enumerate() {
            if let [c, g] = pair
                && c.base == Base::C
                && g.base == Base::G
            {
                annotations.push((index, super::methylation::Cpg { c: c.position, g: g.position }));
            }
        }
        for (index, key) in annotations {
            for offset in [0, 1] {
                let base = bases.get_mut(index + offset).ok_or_eyre("Invalid CpG index")?;
                base.cpg = true;
                base.cpg_key = Some(key);
            }
        }
        bases.retain(|base| (self.start..self.end).contains(&base.position.pos));
        let mut windows = BTreeSet::new();
        if self.localized {
            for molecule in self.fragments.iter().flatten() {
                for &(start, end) in &molecule.intervals {
                    windows.insert((
                        start.saturating_sub(self.padding).max(self.start),
                        end.saturating_add(self.padding).min(self.end),
                    ));
                }
            }
        }
        let mut template = Template { bases, local: BTreeMap::new() };
        template.local = windows
            .into_iter()
            .map(|(start, end)| ((start, end), template.window(start, end)))
            .collect();
        Ok(template)
    }

    #[cfg(test)]
    pub(super) fn fragment_likelihood(
        &self,
        index: usize,
        fragment: &Fragment,
        template: &Template,
        methylation: f64,
        params: &LikelihoodParams,
    ) -> Result<Option<f64>> {
        self.fragment_likelihood_with_model(
            index,
            fragment,
            template,
            params,
            &super::methylation::Scoring::Fixed(methylation),
        )
    }

    /// One fragment on its own. Production [`plan`](Self::plan)s every fragment
    /// and runs one [`Scorer`], which batches the forward passes.
    #[cfg(test)]
    fn fragment_likelihood_with_model(
        &self,
        index: usize,
        fragment: &Fragment,
        template: &Template,
        params: &LikelihoodParams,
        scoring: &super::methylation::Scoring<'_>,
    ) -> Result<Option<f64>> {
        let Some(plan) = self.plan(index, fragment, template, params, scoring)? else {
            return Ok(None);
        };
        let mut scorer = Scorer::new(params)?;
        let mut values = scorer.run(&[plan], false)?;
        values.pop().map(Some).ok_or_eyre("Scorer returned no value")
    }

    /// Everything the forward algorithm has to compute for this fragment under
    /// `template`, or `None` where the fragment falls back to aligned scoring.
    pub(super) fn plan<'a>(
        &'a self,
        index: usize,
        fragment: &Fragment,
        template: &Template,
        params: &LikelihoodParams,
        scoring: &super::methylation::Scoring<'_>,
    ) -> Result<Option<Plan<'a>>> {
        let Some(molecule) = self.fragments.get(index).ok_or_eyre("Missing molecular evidence")?
        else {
            return Ok(None);
        };
        let orientations = super::likelihood::orientations_from(
            fragment
                .reads
                .iter()
                .zip(&molecule.individual)
                .filter(|(_, profile)| !profile.is_empty())
                .map(|(read, _)| read.strand),
        );
        let mut plan = Plan { templates: Vec::new(), jobs: Vec::new(), hypotheses: Vec::new() };
        for &orientation in orientations {
            let mut tracks = Vec::with_capacity(molecule.joint.len());
            for (profile_index, profile) in molecule.joint.iter().enumerate() {
                if profile.is_empty() {
                    tracks.push(None);
                    continue;
                }
                let bases = if self.localized {
                    let &(start, end) = molecule
                        .intervals
                        .get(profile_index)
                        .ok_or_eyre("Missing molecular interval")?;
                    template
                        .local_bases((
                            start.saturating_sub(self.padding).max(self.start),
                            end.saturating_add(self.padding).min(self.end),
                        ))
                        .ok_or_eyre("Missing read-local haplotype template")?
                } else {
                    &template.bases
                };
                let track = template_track(bases, orientation, scoring, params)?;
                let span = molecule
                    .spans
                    .get(profile_index)
                    .ok_or_eyre("Missing molecular profile span")?;
                let offset = bases.partition_point(|base| base.position < span.first);
                let band = if self.localized {
                    band(
                        offset,
                        span.indel_bases
                            .saturating_add(track.indel_bases.max(self.indel_reach))
                            .max(self.padding as usize),
                    )?
                } else {
                    band(0, profile.len().max(bases.len()))?
                };
                plan.templates.push(track);
                tracks.push(Some((plan.templates.len() - 1, band)));
            }
            let mut hypothesis =
                Hypothesis { weight: -(orientations.len() as f64).ln(), jobs: Vec::new() };
            for (profile_index, (profile, terms)) in
                molecule.joint.iter().zip(&molecule.rows).enumerate()
            {
                let job = if profile.is_empty() {
                    Job::Known(0.0)
                } else {
                    let (template, band) = tracks
                        .get(profile_index)
                        .copied()
                        .flatten()
                        .ok_or_eyre("Missing template track")?;
                    Job::Forward { profile: terms, template, band }
                };
                hypothesis.jobs.push(plan.jobs.len());
                plan.jobs.push(job);
            }
            plan.hypotheses.push(hypothesis);
        }
        Ok(Some(plan))
    }
}

/// The forward passes one fragment needs under one template, and how their
/// scores combine over molecular strand orientations.
pub(super) struct Plan<'a> {
    templates: Vec<TemplateTrack>,
    jobs: Vec<Job<'a>>,
    hypotheses: Vec<Hypothesis>,
}

enum Job<'a> {
    Known(f64),
    Forward { profile: &'a ProfileTerms, template: usize, band: compair::Band },
}

/// One molecular strand orientation and the profiles it observes.
struct Hypothesis {
    weight: f64,
    jobs: Vec<usize>,
}

/// Runs plans' forward passes through compair in batches and remembers their
/// results.
///
/// A result is keyed by the full content of its read and template and by its
/// band, so a hit is an exact recomputation avoided, never an approximation.
/// Phasing's beam search rescores the same fragments against the same local
/// templates at every step, which is what the memo exists for.
pub(crate) struct Scorer {
    profiles: Interner<u8>,
    templates: Interner<u64>,
    values: FxHashMap<(u32, u32, i32, u32), f64>,
}

/// Largest footprint of interned keys before the memo starts over.
const MEMO_BYTES: usize = 256 << 20;

struct Progress {
    scores: Vec<Option<f64>>,
}

struct Interner<T> {
    ids: FxHashMap<Box<[T]>, u32>,
    bytes: usize,
}

impl<T> Default for Interner<T> {
    fn default() -> Self {
        Self { ids: FxHashMap::default(), bytes: 0 }
    }
}

impl<T: std::hash::Hash + Eq + Clone> Interner<T> {
    fn get(&self, key: &[T]) -> Option<u32> {
        self.ids.get(key).copied()
    }

    fn intern(&mut self, key: &[T]) -> Result<u32> {
        if let Some(&id) = self.ids.get(key) {
            return Ok(id);
        }
        let id = u32::try_from(self.ids.len())?;
        self.bytes += size_of_val(key);
        self.ids.insert(key.into(), id);
        Ok(id)
    }
}

impl Scorer {
    pub(crate) fn new(params: &LikelihoodParams) -> Result<Self> {
        params.validate()?;
        Ok(Self {
            profiles: Interner::default(),
            templates: Interner::default(),
            values: FxHashMap::default(),
        })
    }

    /// The likelihood of every plan, in order. With `remember`, keys not seen
    /// before join the memo; results are looked up either way.
    pub(super) fn run(&mut self, plans: &[Plan<'_>], remember: bool) -> Result<Vec<f64>> {
        if self.profiles.bytes + self.templates.bytes > MEMO_BYTES {
            self.profiles = Interner::default();
            self.templates = Interner::default();
            self.values.clear();
        }
        let mut progress: Vec<Progress> = plans
            .iter()
            .map(|plan| Progress {
                scores: plan
                    .jobs
                    .iter()
                    .map(|job| match job {
                        Job::Known(value) => Some(*value),
                        Job::Forward { .. } => None,
                    })
                    .collect(),
            })
            .collect();
        let wanted: Vec<_> = plans
            .iter()
            .map(|plan| {
                plan.jobs
                    .iter()
                    .enumerate()
                    .filter_map(|(index, job)| matches!(job, Job::Forward { .. }).then_some(index))
                    .collect()
            })
            .collect();
        self.compute(plans, &wanted, &mut progress, remember, &mut Interned::default())?;
        plans
            .iter()
            .zip(progress)
            .map(|(plan, progress)| {
                let mut result = f64::NEG_INFINITY;
                for hypothesis in &plan.hypotheses {
                    let mut value = hypothesis.weight;
                    for &job in &hypothesis.jobs {
                        value += progress
                            .scores
                            .get(job)
                            .copied()
                            .flatten()
                            .ok_or_eyre("Unscored molecular profile")?;
                    }
                    result = log_add(result, value);
                }
                ensure!(result.is_finite(), "PairHMM fragment likelihood is not finite");
                Ok(result)
            })
            .collect()
    }

    fn compute(
        &mut self,
        plans: &[Plan<'_>],
        wanted: &[Vec<usize>],
        progress: &mut [Progress],
        remember: bool,
        interned: &mut Interned,
    ) -> Result<()> {
        // The jobs the memo cannot answer: their profile, template and band.
        let mut pending = Vec::new();
        let mut keys = Vec::new();
        let mut slots = Vec::new();
        for (plan_index, (plan, wanted)) in plans.iter().zip(wanted).enumerate() {
            for &job in wanted {
                let Some(Job::Forward { profile, template, band }) = plan.jobs.get(job) else {
                    continue;
                };
                let track = plan.templates.get(*template).ok_or_eyre("Missing template track")?;
                let key = interned.key(
                    &mut self.profiles,
                    &mut self.templates,
                    (plan_index, *template),
                    profile,
                    track,
                    remember,
                )?;
                let key =
                    key.map(|(profile, template)| (profile, template, band.offset(), band.width()));
                match key.and_then(|key| self.values.get(&key).copied()) {
                    Some(value) => {
                        progress.score(plan_index, job, value)?;
                    }
                    None => {
                        pending.push((*profile, (plan_index, *template), *band));
                        keys.push(key);
                        slots.push((plan_index, job));
                    }
                }
            }
        }
        if pending.is_empty() {
            return Ok(());
        }
        // Profiles sharing a template use one decoded compair haplotype.
        let mut decoded = FxHashMap::<(usize, usize), usize>::default();
        let mut haplotypes = Vec::new();
        for &(_, (plan_index, template), _) in &pending {
            if let std::collections::hash_map::Entry::Vacant(entry) =
                decoded.entry((plan_index, template))
            {
                let track = plans
                    .get(plan_index)
                    .and_then(|plan| plan.templates.get(template))
                    .ok_or_eyre("Missing template track")?;
                entry.insert(haplotypes.len());
                haplotypes.push(track.haplotype()?);
            }
        }
        let pairs = pending
            .iter()
            .map(|(profile, slot, band)| {
                let haplotype = decoded
                    .get(slot)
                    .and_then(|&index| haplotypes.get(index))
                    .ok_or_eyre("Missing decoded template")?;
                Ok(compair::Pair::new(haplotype, &profile.read, *band))
            })
            .collect::<Result<Vec<_>>>()?;
        let values = align(&pairs)?;
        for ((key, &(plan_index, job)), &value) in keys.iter().zip(&slots).zip(&values) {
            if let Some(key) = key {
                self.values.insert(*key, value);
            }
            progress.score(plan_index, job, value)?;
        }
        Ok(())
    }
}

/// Pairs one rayon task scores before another thread can steal: eight full
/// groups of compair's pairs kernel.
const TASK_PAIRS: usize = 8 * compair::PAIRS;

thread_local! {
    static WORKSPACE: std::cell::RefCell<compair::Workspace> =
        std::cell::RefCell::new(compair::Workspace::new());
}

/// The natural-log forward likelihood of every pair, in order.
///
/// Pairs are scored grouped by band width, which a group of compair's pairs
/// kernel shares, and then by read and template, whose derived tracks it
/// caches; no grouping changes a score. Tasks run on the enclosing rayon pool
/// when there is one, so a worker phasing a dense segment borrows the idle
/// workers.
fn align(pairs: &[compair::Pair<'_>]) -> Result<Vec<f64>> {
    let mut order: Vec<usize> = (0..pairs.len()).collect();
    order.sort_by_key(|&index| {
        pairs.get(index).map(|pair| {
            (
                pair.band.width(),
                std::ptr::from_ref(pair.read).addr(),
                std::ptr::from_ref(pair.haplotype).addr(),
            )
        })
    });
    let sorted = order.iter().map(|&index| pairs.get(index).copied()).collect::<Option<Vec<_>>>();
    let sorted = sorted.ok_or_eyre("PairHMM pair outside its batch")?;
    let run = |chunk: &[compair::Pair<'_>]| -> Vec<compair::Log10Likelihood> {
        let mut out = Vec::with_capacity(chunk.len());
        let mut score = |workspace: &mut compair::Workspace| {
            if chunk.len() >= compair::PAIRS_BREAK_EVEN {
                workspace.align_pairs(chunk, &TemplateEmission, &mut out);
            } else {
                let mut one = Vec::with_capacity(1);
                for pair in chunk {
                    workspace.align_candidates(
                        std::slice::from_ref(pair.haplotype),
                        pair.read,
                        &TemplateEmission,
                        pair.band,
                        &mut one,
                    );
                    out.extend_from_slice(&one);
                }
            }
        };
        WORKSPACE.with(|workspace| match workspace.try_borrow_mut() {
            Ok(mut workspace) => score(&mut workspace),
            Err(_) => score(&mut compair::Workspace::new()),
        });
        out
    };
    let scored: Vec<compair::Log10Likelihood> =
        if rayon::current_thread_index().is_some() && sorted.len() > TASK_PAIRS {
            sorted.par_chunks(TASK_PAIRS).flat_map_iter(run).collect()
        } else {
            sorted.chunks(TASK_PAIRS).flat_map(run).collect()
        };
    ensure!(
        scored.len() == pairs.len(),
        "compair returned {} scores for {} pairs",
        scored.len(),
        pairs.len()
    );
    let mut values = vec![f64::NAN; pairs.len()];
    for (&index, score) in order.iter().zip(scored) {
        let value = score.get() * std::f64::consts::LN_10;
        ensure!(value.is_finite(), "PairHMM has no finite alignment path");
        *values.get_mut(index).ok_or_eyre("PairHMM result outside its batch")? = value;
    }
    Ok(values)
}

trait Scores {
    fn score(&mut self, plan: usize, job: usize, value: f64) -> Result<()>;
}

impl Scores for [Progress] {
    fn score(&mut self, plan: usize, job: usize, value: f64) -> Result<()> {
        *self
            .get_mut(plan)
            .and_then(|progress| progress.scores.get_mut(job))
            .ok_or_eyre("Scored job outside its plan")? = Some(value);
        Ok(())
    }
}

/// Memo ids are keyed by identity within one run; the borrowed plans keep
/// profile and template addresses stable until scoring finishes.
#[derive(Default)]
struct Interned {
    profiles: FxHashMap<usize, Option<u32>>,
    templates: FxHashMap<(usize, usize), Option<u32>>,
}

impl Interned {
    fn key(
        &mut self,
        profiles: &mut Interner<u8>,
        templates: &mut Interner<u64>,
        template_slot: (usize, usize),
        profile: &ProfileTerms,
        template: &TemplateTrack,
        remember: bool,
    ) -> Result<Option<(u32, u32)>> {
        let profile_slot = std::ptr::from_ref(profile).addr();
        let profile_id = match self.profiles.get(&profile_slot) {
            Some(&id) => id,
            None => {
                let id = if remember {
                    Some(profiles.intern(&profile.key)?)
                } else {
                    profiles.get(&profile.key)
                };
                self.profiles.insert(profile_slot, id);
                id
            }
        };
        let template_id = match self.templates.get(&template_slot) {
            Some(&id) => id,
            None => {
                let id = if remember {
                    Some(templates.intern(&template.key)?)
                } else {
                    templates.get(&template.key)
                };
                self.templates.insert(template_slot, id);
                id
            }
        };
        Ok(profile_id.zip(template_id))
    }
}

fn molecular_evidence(
    fragment: &Fragment,
    start: u32,
    end: u32,
    gaps: GapQualities,
    tracks: &RepeatTracks,
) -> Result<Option<MolecularEvidence>> {
    let mut reads = Vec::new();
    for read in &fragment.reads {
        let sequence =
            read.sequence.as_ref().ok_or_eyre("PairHMM requires CIGAR-ordered read sequences")?;
        ensure!(
            sequence.bases.iter().all(|base| base.quality != Some(255)),
            "Missing base quality must be None"
        );
        ensure!(
            sequence.bases.windows(2).all(|pair| match pair {
                [a, b] => a.position < b.position,
                _ => true,
            }),
            "Read sequence coordinates must be unique and ascending"
        );
        for observation in &read.observations {
            let index = sequence
                .bases
                .binary_search_by_key(&SequencePosition::reference(observation.pos), |base| {
                    base.position
                })
                .map_err(|_| {
                    color_eyre::eyre::eyre!("Aligned observation missing from read sequence")
                })?;
            let base = sequence.bases.get(index).ok_or_eyre("Missing read sequence observation")?;
            ensure!(
                base.base == observation.base && base.quality == observation.quality,
                "Aligned observation disagrees with read sequence"
            );
        }
        let mut bases: Vec<_> = sequence
            .bases
            .iter()
            .filter(|base| (start..end).contains(&base.position.pos))
            .copied()
            .collect();
        let first = bases.iter().position(|base| base.base.known_index().is_some());
        let last = bases.iter().rposition(|base| base.base.known_index().is_some());
        bases = match (first, last) {
            (Some(first), Some(last)) => {
                bases.get(first..=last).ok_or_eyre("Invalid masked sequence bounds")?.to_vec()
            }
            _ => Vec::new(),
        };
        if sequence.reference_skip && !bases.is_empty() {
            return Ok(None);
        }
        reads.push(bases);
    }
    let individual: Vec<Profile> = reads
        .iter()
        .enumerate()
        .map(|(read_index, read)| {
            read.iter()
                .map(|base| {
                    consensus([BaseMeasurement {
                        base: base.base,
                        quality: base.quality,
                        mapq: fragment.reads.get(read_index).and_then(|read| read.mapq),
                        read_index,
                    }])
                })
                .collect::<Result<Profile>>()
        })
        .collect::<Result<_>>()?;
    let span = |positions: &[SequencePosition]| -> Result<ProfileSpan> {
        Ok(ProfileSpan {
            first: *positions.first().ok_or_eyre("Empty molecular profile")?,
            indel_bases: indel_bases(positions.iter().copied()),
        })
    };
    let read_spans = |reads: &[Vec<SequenceObservation>]| -> Result<Vec<ProfileSpan>> {
        reads
            .iter()
            .map(|read| {
                if read.is_empty() {
                    Ok(ProfileSpan { first: SequencePosition::reference(0), indel_bases: 0 })
                } else {
                    span(&read.iter().map(|base| base.position).collect::<Vec<_>>())
                }
            })
            .collect()
    };
    let (joint, spans) = match reads.as_slice() {
        [a, b] if !a.is_empty() && !b.is_empty() => {
            let first = a
                .first()
                .ok_or_eyre("Empty read")?
                .position
                .pos
                .max(b.first().ok_or_eyre("Empty read")?.position.pos);
            let last = a
                .last()
                .ok_or_eyre("Empty read")?
                .position
                .pos
                .min(b.last().ok_or_eyre("Empty read")?.position.pos);
            if first <= last {
                let keys = |read: &[SequenceObservation]| {
                    read.iter()
                        .filter(|base| (first..=last).contains(&base.position.pos))
                        .map(|base| base.position)
                        .collect::<Vec<_>>()
                };
                // Conflicting CIGAR paths do not establish which mate bases
                // observe the same chemical event. Keep the aligned model
                // rather than inventing independent conversion draws.
                if keys(a) != keys(b) {
                    return Ok(None);
                }
                let mut merged: BTreeMap<SequencePosition, Vec<BaseMeasurement>> = BTreeMap::new();
                for (read_index, read) in reads.iter().enumerate() {
                    for base in read {
                        merged.entry(base.position).or_default().push(BaseMeasurement {
                            base: base.base,
                            quality: base.quality,
                            mapq: fragment.reads.get(read_index).and_then(|read| read.mapq),
                            read_index,
                        });
                    }
                }
                let positions: Vec<_> = merged.keys().copied().collect();
                (
                    vec![merged.into_values().map(consensus).collect::<Result<Profile>>()?],
                    vec![span(&positions)?],
                )
            } else {
                (individual.clone(), read_spans(&reads)?)
            }
        }
        [_] | [_, _] => (individual.clone(), read_spans(&reads)?),
        _ => color_eyre::eyre::bail!("Fragments require one or two reads"),
    };
    let intervals: Vec<_> =
        if joint.len() == 1 && reads.len() == 2 && reads.iter().all(|read| !read.is_empty()) {
            vec![(
                reads
                    .iter()
                    .filter_map(|read| read.first())
                    .map(|base| base.position.pos)
                    .min()
                    .ok_or_eyre("Missing molecular start")?,
                reads
                    .iter()
                    .filter_map(|read| read.last())
                    .map(|base| base.position.pos)
                    .max()
                    .ok_or_eyre("Missing molecular end")?
                    .checked_add(1)
                    .ok_or_eyre("Molecular interval overflow")?,
            )]
        } else {
            reads
                .iter()
                .map(|read| {
                    let start = read.first().map_or(0, |base| base.position.pos);
                    let end = read.last().map_or(0, |base| base.position.pos.saturating_add(1));
                    (start, end)
                })
                .collect()
        };
    let positions: Vec<Vec<SequencePosition>> = if joint.len() == 1 {
        vec![
            reads
                .iter()
                .flatten()
                .map(|base| base.position)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect(),
        ]
    } else {
        reads.iter().map(|read| read.iter().map(|base| base.position).collect()).collect()
    };
    let rows = joint
        .iter()
        .zip(&positions)
        .map(|(profile, positions)| {
            ensure!(positions.len() == profile.len(), "Profile coordinate count differs");
            profile_terms(
                profile,
                |index| positions.get(index).map_or(tracks.plain, |&position| tracks.at(position)),
                gaps,
                0,
            )
        })
        .collect::<Result<_>>()?;
    Ok(Some(MolecularEvidence { individual, joint, intervals, spans, rows }))
}

#[cfg(test)]
fn score(
    profile: &Profile,
    template: &Template,
    orientation: Strand,
    methylation: f64,
    params: &LikelihoodParams,
) -> Result<f64> {
    score_with_mapping(profile, template, orientation, methylation, params, 0)
}

#[cfg(test)]
fn score_with_mapping(
    profile: &Profile,
    template: &Template,
    orientation: Strand,
    methylation: f64,
    params: &LikelihoodParams,
    bad_reads: usize,
) -> Result<f64> {
    score_with_model(
        profile,
        template,
        orientation,
        &super::methylation::Scoring::Fixed(methylation),
        params,
        bad_reads,
    )
}

/// One profile under one template, through the same plan and scorer the
/// production matrix uses.
#[cfg(test)]
fn score_with_model(
    profile: &Profile,
    template: &Template,
    orientation: Strand,
    scoring: &super::methylation::Scoring<'_>,
    params: &LikelihoodParams,
    bad_reads: usize,
) -> Result<f64> {
    if profile.is_empty() {
        return Ok(0.0);
    }
    let gaps = GapQualities::new(params)?;
    let terms = profile_terms(profile, |_| [gaps.open, gaps.open, gaps.extend], gaps, bad_reads)?;
    // A band holding the whole matrix: these scores stand for the
    // unanchored recurrence.
    let band = band(0, profile.len().max(template.bases.len()))?;
    let plan = Plan {
        templates: vec![template_track(&template.bases, orientation, scoring, params)?],
        jobs: vec![Job::Forward { profile: &terms, template: 0, band }],
        hypotheses: vec![Hypothesis { weight: 0.0, jobs: vec![0] }],
    };
    Scorer::new(params)?.run(&[plan], false)?.pop().ok_or_eyre("Scorer returned no value")
}

/// `gaps_at` gives each profile base its insertion-open, deletion-open and
/// extension qualities; `gaps` are the plain ones, which the memo key omits.
fn profile_terms(
    profile: &[MolecularBase],
    gaps_at: impl Fn(usize) -> [BaseQuality; 3],
    gaps: GapQualities,
    bad_reads: usize,
) -> Result<ProfileTerms> {
    let bases: Vec<Base> = profile.iter().map(|observed| observed.base).collect();
    let qualities: Vec<_> = profile
        .iter()
        .map(|observed| {
            if observed.source.is_some_and(|index| bad_reads & (1 << index) != 0) {
                BaseQuality::from_byte(0)
            } else if observed.base.known_index().is_none() {
                UNINFORMATIVE
            } else {
                BaseQuality::from_byte(observed.quality.unwrap_or(0))
            }
        })
        .collect();
    let (mut insertion, mut deletion, mut extension) = (
        Vec::with_capacity(bases.len()),
        Vec::with_capacity(bases.len()),
        Vec::with_capacity(bases.len()),
    );
    let plain = [gaps.open, gaps.open, gaps.extend];
    let mut key: Vec<u8> = bases
        .iter()
        .zip(&qualities)
        .flat_map(|(base, quality)| [*base as u8, quality.as_byte()])
        .collect();
    for index in 0..bases.len() {
        let [i, d, e] = gaps_at(index);
        insertion.push(i);
        deletion.push(d);
        extension.push(e);
        if [i, d, e] != plain {
            key.extend(u16::try_from(index)?.to_le_bytes());
            key.extend([i.as_byte(), d.as_byte(), e.as_byte()]);
        }
    }
    // Compair's native insertion emission of one; indel costs are transitions.
    let read =
        compair::Read::new(bases, &qualities, &insertion, &deletion, &extension, Strand::OT)?;
    Ok(ProfileTerms { read, key: key.into() })
}

/// A template under one read orientation, as compair scores it.
///
/// Each column's emission is the one the mixed model has always used:
/// `1 - rate` on the unconverted base and `rate` on its conversion product
/// where this orientation can convert it, a half on each at a neutral column,
/// and the base itself elsewhere. An unknown template base matches every read
/// base, as in GATK.
fn template_track(
    template: &[TemplateBase],
    orientation: Strand,
    scoring: &super::methylation::Scoring<'_>,
    params: &LikelihoodParams,
) -> Result<TemplateTrack> {
    let mut key = Vec::with_capacity(2 * template.len());
    for base in template {
        let level = if base.base.known_index().is_none() || params.non_taps {
            None
        } else if let Some((a, b)) = scoring.neutral_bases(base.position, orientation, base.base)
            && a.known_index().is_some()
            && b.known_index().is_some()
        {
            Some(0.5)
        } else {
            [CpgSide::C, CpgSide::G]
                .into_iter()
                .any(|side| side.unmod_base() == base.base && side.strand() == orientation)
                .then(|| params.conversion(base.cpg.then_some(scoring.rate(base.cpg_key.as_ref()))))
        };
        let level = level.map(Probability::new).transpose()?;
        key.push(u64::from(base.base as u8));
        key.push(level.map_or(u64::MAX, |level| level.to_bits()));
    }
    Ok(TemplateTrack {
        key: key.into(),
        indel_bases: indel_bases(template.iter().map(|base| base.position)),
    })
}

#[cfg(test)]
fn forward(
    rows: usize,
    columns: usize,
    params: &LikelihoodParams,
    insertion: impl Fn(usize) -> f64,
    emission: impl Fn(usize, usize) -> f64,
) -> Result<f64> {
    forward_scaled_reference(
        rows,
        columns,
        params,
        |row| insertion(row).exp(),
        |row, column| emission(row, column).exp(),
    )
}

/// The recurrence cell by cell, as production computed it before the
/// lane-batched kernel: the oracle that kernel is held to bit for bit.
#[cfg(test)]
fn forward_scaled_reference(
    rows: usize,
    columns: usize,
    params: &LikelihoodParams,
    insertion: impl Fn(usize) -> f64,
    emission: impl Fn(usize, usize) -> f64,
) -> Result<f64> {
    if rows == 0 {
        return Ok(0.0);
    }
    ensure!(columns > 0, "Empty PairHMM haplotype");
    let width = columns.checked_add(1).ok_or_eyre("PairHMM width overflow")?;
    let mut previous = vec![[0.0; 3]; width];
    for cell in previous.iter_mut().take(columns) {
        if let Some(deletion) = cell.get_mut(2) {
            *deletion = 1.0 / columns as f64;
        }
    }
    let mut current = vec![[0.0; 3]; width];
    let (open, extend) = (*params.gap_open, *params.gap_extend);
    let transitions = [1.0 - 2.0 * open, 1.0 - extend, 1.0 - extend];
    let mut scale = 0.0;
    for row in 0..rows {
        current.fill([0.0; 3]);
        let insertion = insertion(row);
        for column in 1..width {
            let &[matched_above, inserted_above, _] =
                previous.get(column).ok_or_eyre("Missing HMM row cell")?;
            let &[matched_left, _, deleted_left] =
                current.get(column - 1).ok_or_eyre("Missing HMM column cell")?;
            let diagonal = previous.get(column - 1).ok_or_eyre("Missing diagonal HMM cell")?;
            let matched = emission(row, column - 1)
                * diagonal
                    .iter()
                    .zip(transitions)
                    .map(|(value, transition)| value * transition)
                    .sum::<f64>();
            let inserted = insertion * (matched_above * open + inserted_above * extend);
            let deleted = matched_left * open + deleted_left * extend;
            *current.get_mut(column).ok_or_eyre("Missing destination HMM cell")? =
                [matched, inserted, deleted];
        }
        let maximum = current.iter().flatten().copied().fold(0.0, f64::max);
        ensure!(maximum.is_finite() && maximum > 0.0, "PairHMM has no finite alignment path");
        for value in current.iter_mut().flatten() {
            *value /= maximum;
        }
        scale += maximum.ln();
        std::mem::swap(&mut previous, &mut current);
    }
    let mass: f64 = previous.iter().skip(1).flat_map(|cell| cell.iter().take(2)).sum();
    let value = scale + mass.ln();
    ensure!(value.is_finite(), "PairHMM forward likelihood is not finite");
    Ok(value)
}

#[cfg(test)]
fn forward_log(
    rows: usize,
    columns: usize,
    params: &LikelihoodParams,
    insertion: impl Fn(usize) -> f64,
    emission: impl Fn(usize, usize) -> f64,
) -> Result<f64> {
    if rows == 0 {
        return Ok(0.0);
    }
    ensure!(columns > 0, "Empty PairHMM haplotype");
    let width = columns.checked_add(1).ok_or_eyre("PairHMM width overflow")?;
    let mut previous = vec![[f64::NEG_INFINITY; 3]; width];
    for cell in previous.iter_mut().take(columns) {
        if let Some(deletion) = cell.get_mut(2) {
            *deletion = -(columns as f64).ln();
        }
    }
    let mut current = vec![[f64::NEG_INFINITY; 3]; width];
    let (open, extend) = (*params.gap_open, *params.gap_extend);
    let transitions = [(1.0 - 2.0 * open).ln(), (-extend).ln_1p(), (-extend).ln_1p()];
    let (open, extend) = (open.ln(), extend.ln());
    for i in 0..rows {
        current.fill([f64::NEG_INFINITY; 3]);
        for j in 1..width {
            let diagonal = previous.get(j - 1).ok_or_eyre("Missing diagonal HMM cell")?;
            let above = previous.get(j).ok_or_eyre("Missing HMM row cell")?;
            let left = current.get(j - 1).ok_or_eyre("Missing HMM column cell")?;
            let matched = emission(i, j - 1)
                + log_sum(
                    diagonal.iter().zip(transitions).map(|(&value, transition)| value + transition),
                );
            let inserted = insertion(i)
                + log_add(
                    above.first().copied().ok_or_eyre("Missing match state")? + open,
                    above.get(1).copied().ok_or_eyre("Missing insertion state")? + extend,
                );
            let deleted = log_add(
                left.first().copied().ok_or_eyre("Missing match state")? + open,
                left.get(2).copied().ok_or_eyre("Missing deletion state")? + extend,
            );
            *current.get_mut(j).ok_or_eyre("Missing destination HMM cell")? =
                [matched, inserted, deleted];
        }
        std::mem::swap(&mut previous, &mut current);
    }
    let value = log_sum(previous.iter().skip(1).flat_map(|cell| cell.iter().take(2)).copied());
    ensure!(value.is_finite(), "PairHMM forward likelihood is not finite");
    Ok(value)
}

#[cfg(test)]
mod tests;
