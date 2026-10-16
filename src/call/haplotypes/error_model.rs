use color_eyre::eyre::{OptionExt as _, Result, ensure};
use seqair_types::{Base, BaseQuality};

/// Library error parameters are part of the forest's feature definition.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ErrorProfile {
    pub version: u32,
    pub platform: String,
    pub library: String,
    pub maximum_evidence_phred: u8,
    pub repeats: Vec<RepeatError>,
    /// When present, whole-unit copy-number errors of repeat tracts are a
    /// molecule-level mixture over the chromosome's stuttered copies rather
    /// than per-position gap qualities, and `repeats` is not applied: the
    /// forward sum would otherwise count one stutter once per placement, and
    /// one gap quality cannot serve both the presence and the dosage of an
    /// indel in a library whose stutter is this common.
    #[serde(default)]
    pub stutter: Vec<StutterRate>,
}

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StutterRate {
    pub period: u8,
    pub minimum_copies: u16,
    /// Per-molecule probabilities of reading the tract with two units fewer,
    /// one fewer, one more and two more.
    pub steps: [f64; 4],
}

/// The copy-number changes `StutterRate::steps` are given for, in order.
pub(super) const STUTTER_STEPS: [i64; 4] = [-2, -1, 1, 2];

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RepeatError {
    pub period: u8,
    pub minimum_copies: u16,
    pub insertion_open_phred: u8,
    pub deletion_open_phred: u8,
    pub extension_phred: u8,
}

impl Default for ErrorProfile {
    fn default() -> Self {
        // An explicit conservative starting profile, replaced by library calibration.
        // The repeat-dependent opening curve follows GATK 4.5's conservative PCR model.
        let repeats = (1..=6)
            .flat_map(|period| {
                (2..=20).map(move |copies| {
                    let quality = (40.0 - (f64::from(copies) / (3.0 * std::f64::consts::PI)).exp()
                        + 1.0)
                        .round()
                        .clamp(10.0, 40.0);
                    #[expect(
                        clippy::cast_possible_truncation,
                        clippy::cast_sign_loss,
                        reason = "Rounded Phred in 10..=40"
                    )]
                    let quality = quality as u8;
                    RepeatError {
                        period,
                        minimum_copies: copies,
                        insertion_open_phred: quality,
                        deletion_open_phred: quality,
                        extension_phred: 10,
                    }
                })
            })
            .collect();
        Self {
            version: 1,
            platform: "illumina".into(),
            library: "conservative".into(),
            maximum_evidence_phred: 45,
            repeats,
            stutter: Vec::new(),
        }
    }
}

impl ErrorProfile {
    pub(super) fn models_stutter(&self) -> bool {
        !self.stutter.is_empty()
    }

    /// The stutter steps of a tract of `period` with `copies` copies: the
    /// entry for the most copies at or under it, so long tracts keep the
    /// longest calibrated context.
    pub(super) fn stutter_at(&self, period: usize, copies: usize) -> Option<[f64; 4]> {
        self.stutter
            .iter()
            .filter(|entry| {
                usize::from(entry.period) == period && usize::from(entry.minimum_copies) <= copies
            })
            .max_by_key(|entry| entry.minimum_copies)
            .map(|entry| entry.steps)
    }

    pub fn validate(&self) -> Result<()> {
        let mut contexts = std::collections::BTreeSet::new();
        for entry in &self.stutter {
            ensure!(
                (1..=6).contains(&entry.period) && entry.minimum_copies >= 2,
                "Invalid stutter repeat context"
            );
            ensure!(
                contexts.insert((entry.period, entry.minimum_copies)),
                "Duplicate stutter context"
            );
            ensure!(
                entry.steps.iter().all(|step| step.is_finite() && (0.0..=1.0).contains(step))
                    && entry.steps.iter().sum::<f64>() <= 1.0,
                "Stutter steps must be probabilities summing to at most one"
            );
        }
        ensure!(self.version == 1, "Unsupported haplotype error profile version");
        ensure!(
            !self.platform.is_empty() && !self.library.is_empty(),
            "Error profile needs platform and library names"
        );
        ensure!(
            (1..=99).contains(&self.maximum_evidence_phred),
            "Molecular evidence cap must be Q1..Q99"
        );
        let mut keys = std::collections::BTreeSet::new();
        for entry in &self.repeats {
            ensure!(
                (1..=6).contains(&entry.period) && entry.minimum_copies >= 2,
                "Invalid repeat context"
            );
            ensure!(
                keys.insert((entry.period, entry.minimum_copies)),
                "Duplicate repeat error context"
            );
            ensure!(
                (6..=60).contains(&entry.insertion_open_phred)
                    && (6..=60).contains(&entry.deletion_open_phred)
                    && (1..=60).contains(&entry.extension_phred),
                "Invalid repeat gap qualities"
            );
        }
        Ok(())
    }

    pub(super) fn tracks(
        &self,
        bases: &[Base],
        open: BaseQuality,
        extend: BaseQuality,
    ) -> Vec<[BaseQuality; 3]> {
        let mut tracks = vec![[open, open, extend]; bases.len()];
        if self.models_stutter() {
            return tracks;
        }
        let mut spans = vec![0usize; bases.len()];
        for period in 1usize..=6 {
            let mut start = period;
            while start < bases.len() {
                let matches = |index: usize| {
                    bases
                        .get(index)
                        .zip(bases.get(index - period))
                        .is_some_and(|(a, b)| a.known_index().is_some() && a == b)
                };
                if !matches(start) {
                    start += 1;
                    continue;
                }
                let mut end = start + 1;
                while end < bases.len() && matches(end) {
                    end += 1;
                }
                let motif = bases.get(start - period..start).unwrap_or_default();
                if (1..period).any(|shorter| {
                    period.is_multiple_of(shorter)
                        && motif
                            .iter()
                            .enumerate()
                            .all(|(index, base)| motif.get(index % shorter) == Some(base))
                }) {
                    start = end;
                    continue;
                }
                let copies = (end - start + period) / period;
                if let Some(entry) = self
                    .repeats
                    .iter()
                    .filter(|entry| {
                        usize::from(entry.period) == period
                            && usize::from(entry.minimum_copies) <= copies
                    })
                    .max_by_key(|entry| entry.minimum_copies)
                {
                    for (track, span) in
                        tracks.iter_mut().zip(&mut spans).take(end).skip(start - period)
                    {
                        let length = end - start + period;
                        if length > *span {
                            *span = length;
                            *track = [
                                entry.insertion_open_phred,
                                entry.deletion_open_phred,
                                entry.extension_phred,
                            ]
                            .map(BaseQuality::from_byte);
                        }
                    }
                }
                start = end;
            }
        }
        tracks
    }
}

/// A repeat tract one of the compared indels changes the copy number of.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Tract {
    start: u32,
    unit: Vec<Base>,
    /// Copies on the reference.
    copies: usize,
}

/// The smallest period `sequence` repeats with, its own length if none.
fn primitive_period(sequence: &[Base]) -> usize {
    (1..sequence.len())
        .find(|&period| {
            sequence.len().is_multiple_of(period)
                && sequence
                    .iter()
                    .enumerate()
                    .all(|(index, base)| sequence.get(index % period) == Some(base))
        })
        .unwrap_or(sequence.len())
}

fn tract_of(reference: &super::Reference, indel: &super::Indel) -> Result<Option<Tract>> {
    let indel = indel.normalized(reference)?;
    let sequence = if indel.deleted.is_empty() { &indel.inserted } else { &indel.deleted };
    let period = primitive_period(sequence);
    if period == 0 {
        return Ok(None);
    }
    let end = indel.equivalent_end(reference)?;
    let bases = usize::try_from(end.saturating_sub(indel.start))?
        .saturating_sub(usize::from(indel.deleted.is_empty()));
    let copies = bases / period;
    Ok((copies >= 2).then(|| Tract {
        start: indel.start,
        unit: sequence.get(..period).unwrap_or_default().to_vec(),
        copies,
    }))
}

/// How a haplotype's own indels change `tract`'s copy number, in units, or
/// `None` where one of them is not a whole number of units: that chromosome
/// has no simple stuttered copies.
fn units_at(
    reference: &super::Reference,
    haplotype: &super::Haplotype,
    tract: &Tract,
) -> Result<Option<i64>> {
    let end = tract
        .start
        .checked_add(u32::try_from(tract.copies * tract.unit.len())?)
        .ok_or_eyre("Tract overflow")?;
    let mut units = 0i64;
    for edit in &haplotype.indels {
        let edit = edit.normalized(reference)?;
        let edit_end = edit
            .start
            .checked_add(u32::try_from(edit.deleted.len().max(1))?)
            .ok_or_eyre("Indel overflow")?;
        if edit_end <= tract.start || edit.start >= end {
            continue;
        }
        let sequence = if edit.deleted.is_empty() { &edit.inserted } else { &edit.deleted };
        let whole = edit.start == tract.start
            && sequence.len().is_multiple_of(tract.unit.len())
            && sequence.chunks(tract.unit.len()).all(|chunk| chunk == tract.unit.as_slice());
        if !whole {
            return Ok(None);
        }
        let count = i64::try_from(sequence.len() / tract.unit.len())?;
        units += if edit.deleted.is_empty() { count } else { -count };
    }
    Ok(Some(units))
}

/// `haplotype` read with `tract`'s copy number changed by `step` units.
fn stuttered(
    reference: &super::Reference,
    haplotype: &super::Haplotype,
    tract: &Tract,
    units: i64,
    step: i64,
) -> Result<Option<super::Haplotype>> {
    let net = units + step;
    let copies = i64::try_from(tract.copies)? + net;
    if copies < 0 {
        return Ok(None);
    }
    let end = tract
        .start
        .checked_add(u32::try_from(tract.copies * tract.unit.len())?)
        .ok_or_eyre("Tract overflow")?;
    let mut result = haplotype.clone();
    result.indels.retain(|edit| {
        let edit_end =
            edit.start.saturating_add(u32::try_from(edit.deleted.len().max(1)).unwrap_or(u32::MAX));
        edit_end <= tract.start || edit.start >= end
    });
    if net != 0 {
        let length = usize::try_from(net.unsigned_abs())? * tract.unit.len();
        let edit = if net < 0 {
            let deleted = (0..length)
                .map(|offset| {
                    reference.base_at(tract.start.checked_add(u32::try_from(offset).ok()?)?)
                })
                .collect::<Option<Vec<_>>>()
                .ok_or_eyre("Stutter deletion leaves the reference")?;
            super::Indel { start: tract.start, deleted, inserted: Vec::new() }
        } else {
            super::Indel {
                start: tract.start,
                deleted: Vec::new(),
                inserted: tract.unit.iter().copied().cycle().take(length).collect(),
            }
        };
        result.alleles.retain(|&(pos, _)| !edit.deletes(pos));
        result.indels.push(edit);
        result.indels.sort();
    }
    Ok(Some(result))
}

/// The chromosomes a likelihood pass scores when the profile models stutter:
/// every compared haplotype, then every distinct stuttered copy of one, and
/// for each compared haplotype the rows its molecule likelihood mixes with
/// their weights.
pub(super) struct StutterMixture {
    pub extended: Vec<super::Haplotype>,
    pub mixtures: Vec<Vec<(usize, f64)>>,
}

pub(super) fn stutter_mixture(
    reference: &super::Reference,
    haplotypes: &[super::Haplotype],
    indel: Option<&super::Indel>,
    profile: &ErrorProfile,
) -> Result<Option<StutterMixture>> {
    if !profile.models_stutter() {
        return Ok(None);
    }
    let mut tracts: Vec<Tract> = Vec::new();
    for edit in haplotypes.iter().flat_map(|haplotype| &haplotype.indels).chain(indel) {
        if let Some(tract) = tract_of(reference, edit)?
            && !tracts.contains(&tract)
        {
            tracts.push(tract);
        }
    }
    if tracts.is_empty() {
        return Ok(None);
    }
    let same =
        |a: &super::Haplotype, b: &super::Haplotype| a.alleles == b.alleles && a.indels == b.indels;
    let mut extended: Vec<super::Haplotype> = haplotypes.to_vec();
    let mut mixtures = Vec::with_capacity(haplotypes.len());
    let mut any = false;
    for (index, haplotype) in haplotypes.iter().enumerate() {
        let mut mixture = Vec::new();
        let mut remaining = 1.0;
        for tract in &tracts {
            let Some(units) = units_at(reference, haplotype, tract)? else { continue };
            let copies = usize::try_from(i64::try_from(tract.copies)? + units).unwrap_or(0);
            let Some(steps) = profile.stutter_at(tract.unit.len(), copies) else { continue };
            for (&step, &rate) in STUTTER_STEPS.iter().zip(&steps) {
                if rate <= 0.0 {
                    continue;
                }
                let Some(variant) = stuttered(reference, haplotype, tract, units, step)? else {
                    continue;
                };
                let row = match extended.iter().position(|other| same(other, &variant)) {
                    Some(row) => row,
                    None => {
                        extended.push(variant);
                        extended.len() - 1
                    }
                };
                mixture.push((row, rate));
                remaining -= rate;
                any = true;
            }
        }
        ensure!(remaining > 0.0, "Stutter rates of a chromosome exceed one");
        mixture.insert(0, (index, remaining));
        mixtures.push(mixture);
    }
    Ok(any.then_some(StutterMixture { extended, mixtures }))
}

pub(super) fn cap_molecules(
    matrix: &mut [Vec<f64>],
    fragments: &[&super::Fragment],
    reference: &super::Reference,
    params: &super::LikelihoodParams,
) -> Result<()> {
    for (index, fragment) in fragments.iter().enumerate() {
        let best = matrix
            .iter()
            .filter_map(|row| row.get(index))
            .copied()
            .reduce(f64::max)
            .ok_or_else(|| color_eyre::eyre::eyre!("Empty molecular likelihood matrix"))?;
        let active = fragment.reads.iter().filter(|read| {
            read.observations.iter().any(|obs| reference.base_at(obs.pos).is_some())
        });
        let mut quality = params.error_profile.maximum_evidence_phred;
        let mut informative = false;
        for read in active {
            quality = quality.min(read.mapq.unwrap_or(0));
            informative |= read.observations.iter().any(|obs| {
                reference.base_at(obs.pos).is_some()
                    && obs.base.known_index().is_some()
                    && obs.quality.is_some_and(|q| q >= 2)
            });
        }
        if !informative {
            quality = 0;
        }
        let floor = best - f64::from(quality) * std::f64::consts::LN_10 / 10.0;
        for row in matrix.iter_mut() {
            let value = row.get_mut(index).ok_or_else(|| {
                color_eyre::eyre::eyre!("Molecular likelihood dimensions disagree")
            })?;
            ensure!(value.is_finite(), "Non-finite molecular likelihood");
            *value = value.max(floor);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primitive_dna_repeats_select_separate_insertion_and_deletion_errors() -> Result<()> {
        let profile = ErrorProfile {
            repeats: vec![
                RepeatError {
                    period: 1,
                    minimum_copies: 6,
                    insertion_open_phred: 20,
                    deletion_open_phred: 15,
                    extension_phred: 7,
                },
                RepeatError {
                    period: 2,
                    minimum_copies: 3,
                    insertion_open_phred: 6,
                    deletion_open_phred: 6,
                    extension_phred: 1,
                },
            ],
            ..Default::default()
        };
        profile.validate()?;
        let dna: Vec<_> = b"CAAAAAAG".iter().copied().map(Base::from).collect();
        let tracks = profile.tracks(&dna, BaseQuality::from_byte(45), BaseQuality::from_byte(10));
        assert_eq!(tracks.get(3).map(|q| q.map(BaseQuality::as_byte)), Some([20, 15, 7]));
        assert_eq!(tracks.first().map(|q| q.map(BaseQuality::as_byte)), Some([45, 45, 10]));
        let encoded = serde_json::to_vec(&profile)?;
        assert_eq!(serde_json::from_slice::<ErrorProfile>(&encoded)?, profile);
        Ok(())
    }

    #[test]
    fn stuttered_copies_of_a_homopolymer_share_rows_and_rates_follow_each_chromosome() -> Result<()>
    {
        use super::super::{Haplotype, Indel, Reference};
        let mut bases = vec![Base::C; 20];
        bases.extend(vec![Base::A; 14]);
        bases.extend(vec![Base::G; 20]);
        let reference = Reference { start: 1000, bases };
        let deletion = Indel { start: 1020, deleted: vec![Base::A], inserted: Vec::new() };
        let haplotypes = [
            Haplotype { alleles: Vec::new(), indels: Vec::new() },
            Haplotype { alleles: Vec::new(), indels: vec![deletion.clone()] },
        ];
        let mut profile = ErrorProfile::default();
        assert!(stutter_mixture(&reference, &haplotypes, None, &profile)?.is_none());
        profile.stutter = vec![
            StutterRate { period: 1, minimum_copies: 2, steps: [0.0, 0.01, 0.005, 0.0] },
            StutterRate { period: 1, minimum_copies: 14, steps: [0.02, 0.2, 0.1, 0.01] },
        ];
        profile.validate()?;
        let mixture =
            stutter_mixture(&reference, &haplotypes, None, &profile)?.ok_or_eyre("mixture")?;
        let find = |alleles: &[(u32, Base)], deleted: usize, inserted: usize| {
            mixture.extended.iter().position(|h| {
                h.alleles == alleles
                    && h.indels.iter().map(|e| e.deleted.len()).sum::<usize>() == deleted
                    && h.indels.iter().map(|e| e.inserted.len()).sum::<usize>() == inserted
            })
        };
        // The reference has fourteen copies: its one-unit loss is the ALT
        // chromosome itself; ALT has thirteen, so its steps are the shorter
        // context's, and its one-unit gain is the reference.
        let reference_row = mixture.mixtures.first().ok_or_eyre("reference mixture")?;
        assert!(reference_row.contains(&(1usize, 0.2)));
        assert!(reference_row.contains(&(find(&[], 2, 0).ok_or_eyre("two fewer")?, 0.02)));
        assert!(reference_row.contains(&(find(&[], 0, 1).ok_or_eyre("one more")?, 0.1)));
        assert!((reference_row.iter().map(|(_, w)| w).sum::<f64>() - 1.0).abs() < 1e-12);
        let alt_row = mixture.mixtures.get(1).ok_or_eyre("alt mixture")?;
        assert!(alt_row.contains(&(0usize, 0.005)));
        assert!(alt_row.contains(&(find(&[], 2, 0).ok_or_eyre("loss")?, 0.01)));
        assert_eq!(alt_row.first(), Some(&(1usize, 1.0 - 0.015)));
        let distinct: std::collections::BTreeSet<_> = mixture
            .extended
            .iter()
            .map(|h| {
                (
                    h.alleles.iter().map(|&(pos, base)| (pos, base as u8)).collect::<Vec<_>>(),
                    h.indels.clone(),
                )
            })
            .collect();
        assert_eq!(distinct.len(), mixture.extended.len());
        // REF, ALT, and the reference with two fewer, one more and two more
        // copies; ALT's own steps land on REF and on the two-fewer row.
        assert_eq!(mixture.extended.len(), 5);
        Ok(())
    }

    #[test]
    fn stutter_tracts_need_whole_repeat_units() -> Result<()> {
        use super::super::{Indel, Reference};
        let bases: Vec<Base> = b"CCCCATATATATATGGGG".iter().copied().map(Base::from).collect();
        let reference = Reference { start: 0, bases };
        let dinucleotide =
            Indel { start: 4, deleted: vec![Base::A, Base::T], inserted: Vec::new() };
        let tract = tract_of(&reference, &dinucleotide)?.ok_or_eyre("tract")?;
        assert_eq!((tract.start, tract.unit.clone(), tract.copies), (4, vec![Base::A, Base::T], 5));
        let partial = Indel { start: 4, deleted: vec![Base::A], inserted: Vec::new() };
        assert!(tract_of(&reference, &partial)?.is_none());
        assert_eq!(primitive_period(&[Base::A, Base::T, Base::A, Base::T]), 2);
        assert_eq!(primitive_period(&[Base::A, Base::T, Base::A]), 3);
        Ok(())
    }

    #[test]
    fn molecular_caps_are_symmetric_and_mapping_quality_can_neutralize_evidence() -> Result<()> {
        use super::super::{Fragment, Observation, ReadEvidence, Reference};
        let reference = Reference { start: 1000, bases: vec![Base::A] };
        let make = |mapq| Fragment {
            reads: vec![ReadEvidence {
                mapq,
                strand: seqair_types::Strand::OT,
                observations: vec![Observation { pos: 1000, base: Base::A, quality: Some(40) }],
                sequence: None,
                annotations: None,
            }],
        };
        let fragments = [make(Some(60)), make(Some(60)), make(None)];
        let mut matrix = vec![vec![-1.0, -100.0, -1.0], vec![-100.0, -1.0, -100.0]];
        let params = super::super::features::parameters()?;
        cap_molecules(&mut matrix, &fragments.iter().collect::<Vec<_>>(), &reference, &params)?;
        let delta =
            f64::from(params.error_profile.maximum_evidence_phred) * std::f64::consts::LN_10 / 10.0;
        assert!(
            (matrix.first().and_then(|row| row.get(1)).copied().expect("score") + 1.0 + delta)
                .abs()
                < 1e-10
        );
        assert_eq!(
            matrix.first().and_then(|row| row.get(1)),
            matrix.get(1).and_then(|row| row.first())
        );
        assert_eq!(
            matrix.first().and_then(|row| row.get(2)),
            matrix.get(1).and_then(|row| row.get(2))
        );
        Ok(())
    }
}
