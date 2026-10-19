use color_eyre::eyre::{Result, ensure};
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
}

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
        }
    }
}

impl ErrorProfile {
    pub fn validate(&self) -> Result<()> {
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

    /// Per base of `bases`, the insertion-open, deletion-open and extension
    /// qualities of the longest calibrated repeat tract it lies in, or the
    /// plain ones outside any.
    pub(crate) fn tracks(
        &self,
        bases: &[Base],
        open: BaseQuality,
        extend: BaseQuality,
    ) -> Vec<[BaseQuality; 3]> {
        let mut tracks = vec![[open, open, extend]; bases.len()];
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
    fn molecular_caps_are_symmetric_and_mapping_quality_can_neutralize_evidence() -> Result<()> {
        use super::super::{Fragment, Observation, ReadEvidence, Reference};
        let reference = Reference { start: 1000, bases: vec![Base::A] };
        let make = |mapq| Fragment {
            reads: vec![ReadEvidence {
                mapq,
                strand: seqair_types::Strand::OT,
                observations: vec![Observation { pos: 1000, base: Base::A, quality: Some(40) }],
                sequence: None,
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
