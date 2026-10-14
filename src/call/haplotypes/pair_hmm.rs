#[cfg(test)]
use super::likelihood::log_sum;
use super::likelihood::{BaseMeasurement, MolecularBase, consensus, log_add, sequencing_emission};
use super::{Fragment, Haplotype, LikelihoodParams, Reference, Site};
use crate::metrics::methylation::CpgSide;
use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, Strand};
use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet},
};

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
    path_masses: BTreeMap<usize, f64>,
    local: BTreeMap<(u32, u32), Box<Template>>,
}

fn sequence(
    reference: &Reference,
    haplotype: &Haplotype,
    start: u32,
    end: u32,
) -> Result<Vec<TemplateBase>> {
    let mut bases = Vec::new();
    for pos in start..end {
        if let Some(indel) = &haplotype.indel {
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
                continue;
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
    Ok(sequence(reference, haplotype, start, end)?
        .windows(2)
        .filter_map(|pair| match pair {
            [c, g] if c.base == Base::C && g.base == Base::G => {
                Some(super::methylation::Cpg { c: c.position, g: g.position })
            }
            _ => None,
        })
        .collect())
}

type Profile = Vec<MolecularBase>;

struct MolecularEvidence {
    individual: Vec<Profile>,
    joint: Vec<Profile>,
    intervals: Vec<(u32, u32)>,
}

pub(super) struct Prepared {
    start: u32,
    end: u32,
    fragments: Vec<Option<MolecularEvidence>>,
    lengths: BTreeSet<usize>,
    localized: bool,
    normalizers: RefCell<BTreeMap<(usize, usize, u64, u64), f64>>,
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
        let mut prepared = Vec::with_capacity(fragments.len());
        let mut lengths = BTreeSet::new();
        for fragment in fragments {
            let molecule = molecular_evidence(fragment, start, end)?;
            if let Some(molecule) = &molecule {
                lengths.extend(molecule.joint.iter().map(Vec::len));
            }
            prepared.push(molecule);
        }
        Ok(Self {
            start,
            end,
            fragments: prepared,
            lengths,
            localized: false,
            normalizers: RefCell::new(BTreeMap::new()),
        })
    }

    #[cfg(any(feature = "experimental-seqair", test))]
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
        let intervals: Vec<_> = self
            .fragments
            .iter()
            .flatten()
            .flat_map(|molecule| &molecule.intervals)
            .copied()
            .collect();
        haplotype
            .alleles
            .iter()
            .filter(|(pos, _)| {
                if self.localized && !intervals.is_empty() {
                    intervals.iter().any(|&(start, end)| {
                        start
                            .saturating_sub(u32::from(params.realignment_flank) + 1)
                            .max(self.start.saturating_sub(1))
                            <= *pos
                            && *pos
                                <= end
                                    .saturating_add(u32::from(params.realignment_flank))
                                    .min(self.end)
                    })
                } else {
                    self.start.saturating_sub(1) <= *pos && *pos <= self.end
                }
            })
            .map(|&(pos, base)| (pos, base as u8))
            .collect()
    }

    fn normalization(&self, rows: usize, columns: usize, params: &LikelihoodParams) -> Result<f64> {
        let key = (rows, columns, params.gap_open.to_bits(), params.gap_extend.to_bits());
        if let Some(&value) = self.normalizers.borrow().get(&key) {
            return Ok(value);
        }
        let value = forward(rows, columns, params, |_| 0.0, |_, _| 0.0)?;
        self.normalizers.borrow_mut().insert(key, value);
        Ok(value)
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
        params: &LikelihoodParams,
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
        let mut local = BTreeMap::new();
        if self.localized {
            for molecule in self.fragments.iter().flatten() {
                for (profile, &(start, end)) in molecule.joint.iter().zip(&molecule.intervals) {
                    let start =
                        start.saturating_sub(u32::from(params.realignment_flank)).max(self.start);
                    let end = end.saturating_add(u32::from(params.realignment_flank)).min(self.end);
                    let subset = bases
                        .iter()
                        .filter(|base| (start..end).contains(&base.position.pos))
                        .copied()
                        .collect::<Vec<_>>();
                    let child = local.entry((start, end)).or_insert_with(|| {
                        Box::new(Template {
                            bases: subset,
                            path_masses: BTreeMap::new(),
                            local: BTreeMap::new(),
                        })
                    });
                    if !profile.is_empty() && !child.path_masses.contains_key(&profile.len()) {
                        child.path_masses.insert(
                            profile.len(),
                            self.normalization(profile.len(), child.bases.len(), params)?,
                        );
                    }
                }
            }
        }
        let path_masses = if self.localized {
            BTreeMap::new()
        } else {
            self.lengths
                .iter()
                .map(|&length| Ok((length, self.normalization(length, bases.len(), params)?)))
                .collect::<Result<BTreeMap<_, _>>>()?
        };
        Ok(Template { bases, path_masses, local })
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

    pub(super) fn fragment_likelihood_with_model(
        &self,
        index: usize,
        fragment: &Fragment,
        template: &Template,
        params: &LikelihoodParams,
        scoring: &super::methylation::Scoring<'_>,
    ) -> Result<Option<f64>> {
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
        let mut result = f64::NEG_INFINITY;
        for &orientation in orientations {
            let mut profile_scores = vec![[None::<f64>; 4]; molecule.joint.len()];
            for mask in 0..(1usize << fragment.reads.len()) {
                let owns = |read_index| {
                    molecule.joint.iter().flatten().any(|base| base.source == Some(read_index))
                };
                if fragment
                    .reads
                    .iter()
                    .enumerate()
                    .any(|(read_index, _)| !owns(read_index) && mask & (1 << read_index) != 0)
                {
                    continue;
                }
                let mut weight = -(orientations.len() as f64).ln();
                for (read_index, read) in fragment.reads.iter().enumerate() {
                    if !owns(read_index) {
                        continue;
                    }
                    let bad = mask & (1 << read_index) != 0;
                    let error = read.mapq.map_or(1.0, |q| 10.0_f64.powf(-f64::from(q) / 10.0));
                    weight += if bad { error.ln() } else { (-error).ln_1p() };
                }
                if !weight.is_finite() {
                    continue;
                }
                let likelihood = molecule
                    .joint
                    .iter()
                    .enumerate()
                    .map(|(index, profile)| {
                        let sources = profile
                            .iter()
                            .filter_map(|base| base.source)
                            .fold(0usize, |sources, source| sources | (1 << source));
                        let local_mask = mask & sources;
                        let cached = profile_scores
                            .get_mut(index)
                            .and_then(|scores| scores.get_mut(local_mask))
                            .ok_or_eyre("Invalid profile mapping hypothesis")?;
                        if let Some(score) = cached {
                            return Ok(*score);
                        }
                        let local_template;
                        let template = if self.localized && !profile.is_empty() {
                            let &(start, end) = molecule
                                .intervals
                                .get(index)
                                .ok_or_eyre("Missing molecular interval")?;
                            local_template = (
                                start
                                    .saturating_sub(u32::from(params.realignment_flank))
                                    .max(self.start),
                                end.saturating_add(u32::from(params.realignment_flank))
                                    .min(self.end),
                            );
                            template
                                .local
                                .get(&local_template)
                                .ok_or_eyre("Missing read-local haplotype template")?
                                .as_ref()
                        } else {
                            template
                        };
                        let score = score_with_model(
                            profile,
                            template,
                            orientation,
                            scoring,
                            params,
                            local_mask,
                        )?;
                        *cached = Some(score);
                        Ok(score)
                    })
                    .collect::<Result<Vec<_>>>()?
                    .iter()
                    .sum::<f64>();
                result = log_add(result, weight + likelihood);
            }
        }
        ensure!(result.is_finite(), "PairHMM fragment likelihood is not finite");
        Ok(Some(result))
    }
}

fn molecular_evidence(
    fragment: &Fragment,
    start: u32,
    end: u32,
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
    let joint = match reads.as_slice() {
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
                vec![merged.into_values().map(consensus).collect::<Result<Profile>>()?]
            } else {
                individual.clone()
            }
        }
        [_] | [_, _] => individual.clone(),
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
    Ok(Some(MolecularEvidence { individual, joint, intervals }))
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
    if profile.iter().all(|base| base.source.is_none_or(|index| bad_reads & (1 << index) != 0)) {
        // Every path emits the same uniform observation probability. Path-mass
        // normalization cancels the entire alignment sum, including indels.
        return Ok(profile.iter().filter(|base| base.base.known_index().is_some()).count() as f64
            * 0.25_f64.ln());
    }
    let insertions: Vec<_> = profile
        .iter()
        .map(|observed| if observed.base.known_index().is_some() { 0.25 } else { 1.0 })
        .collect();
    let emissions: Vec<[f64; 4]> = profile
        .iter()
        .map(|observed| {
            Base::KNOWN.map(|base| {
                if observed.source.is_some_and(|index| bad_reads & (1 << index) != 0) {
                    0.25
                } else {
                    sequencing_emission(observed.base, base, observed.quality)
                }
            })
        })
        .collect();
    let conversions: Vec<_> = template
        .bases
        .iter()
        .map(|base| {
            let side = [CpgSide::C, CpgSide::G].into_iter().find(|side| {
                side.unmod_base() == base.base && side.strand() == orientation && !params.non_taps
            });
            let rate = params.conversion(base.cpg.then_some(scoring.rate(base.cpg_key.as_ref())));
            (
                base.base.known_index(),
                side.and_then(|side| side.mod_base().known_index()),
                1.0 - rate,
                rate,
                if params.non_taps {
                    None
                } else {
                    scoring.neutral_bases(base.position, orientation, base.base)
                },
            )
        })
        .collect();
    let insertion = |i: usize| insertions.get(i).copied().unwrap_or(0.0);
    let emission = |i: usize, j: usize| {
        let Some(&(ordinary, modified, unconverted_rate, converted_rate, masked)) =
            conversions.get(j)
        else {
            return 0.0;
        };
        let Some(ordinary) = ordinary else { return insertion(i) };
        let Some(values) = emissions.get(i) else { return 0.0 };
        if let Some((a, b)) = masked {
            let Some(a) = a.known_index().and_then(|index| values.get(index)) else {
                return 0.0;
            };
            let Some(b) = b.known_index().and_then(|index| values.get(index)) else {
                return 0.0;
            };
            return (*a + *b) * 0.5;
        }
        let ordinary = values.get(ordinary).copied().unwrap_or(0.0);
        let Some(modified) = modified.and_then(|index| values.get(index).copied()) else {
            return ordinary;
        };
        unconverted_rate * ordinary + converted_rate * modified
    };
    let value = forward_scaled(profile.len(), template.bases.len(), params, insertion, emission)?;
    // Condition on the observed profile length and containment in this local
    // window. Uniform bases then carry no evidence for a longer haplotype.
    let mass = *template
        .path_masses
        .get(&profile.len())
        .ok_or_eyre("Missing alignment path normalization")?;
    Ok(value - mass)
}

fn forward(
    rows: usize,
    columns: usize,
    params: &LikelihoodParams,
    insertion: impl Fn(usize) -> f64,
    emission: impl Fn(usize, usize) -> f64,
) -> Result<f64> {
    forward_scaled(
        rows,
        columns,
        params,
        |row| insertion(row).exp(),
        |row, column| emission(row, column).exp(),
    )
}

fn forward_scaled(
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
        // A common row scale leaves every relative path weight intact and prevents underflow on long reads.
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
