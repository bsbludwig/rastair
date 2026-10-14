use super::{
    Fragment, Indel, Observation, ReadEvidence, Reference, SequenceEvidence, SequenceObservation,
    SequencePosition,
};
use crate::call::{
    pileup::{
        has_terminal_repeat,
        indels::{IndelAllele, IndelObservation},
    },
    process::PileupMappingParams,
};
#[cfg(not(feature = "experimental-seqair"))]
use color_eyre::eyre::WrapErr as _;
use color_eyre::eyre::{OptionExt as _, Result, bail, ensure};
use seqair::bam::cigar::{CigarOp, CigarOpType};
#[cfg(feature = "experimental-seqair")]
use seqair::bam::record_store::RecordStore;
use seqair_types::{Base, SmallVec, Strand};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone)]
pub struct NormalizedObservation {
    pub raw_start: u32,
    pub indel: Indel,
    pub observation: IndelObservation,
    pub boundary_aligned: bool,
}

#[derive(Debug, Clone)]
pub struct ReadAnnotations {
    pub indels: Vec<NormalizedObservation>,
    pub soft_clipped: bool,
    pub noisy: bool,
}

#[derive(Default)]
pub(crate) struct Collector {
    named: BTreeMap<(Vec<u8>, Vec<u8>), Vec<ReadEvidence>>,
    seen: BTreeSet<(Vec<u8>, Vec<u8>, u64, u16)>,
    pub nameless: bool,
}

struct Alignment<'a> {
    name: &'a [u8],
    group: Vec<u8>,
    start: u32,
    flags: u16,
    mapq: u8,
    strand: Strand,
    bases: &'a [Base],
    qualities: &'a [u8],
    cigar: &'a [CigarOp],
}

impl Collector {
    #[cfg(feature = "experimental-seqair")]
    pub fn capture(
        &mut self,
        store: &RecordStore<crate::sequence::RastairReadExtras>,
        reference: &Reference,
        params: &PileupMappingParams,
    ) -> Result<()> {
        for record in store.records() {
            let name = record.qname(store)?;
            let group = match record.aux(store)?.get::<&[u8]>(*b"RG") {
                Ok(group) => group.to_vec(),
                Err(seqair::bam::aux::GetAuxError::TagNotFound { .. }) => Vec::new(),
                Err(error) => return Err(error.into()),
            };
            let qualities: Vec<_> =
                record.qual(store)?.iter().map(|q| q.get().unwrap_or(255)).collect();
            self.add(
                Alignment {
                    name,
                    group,
                    start: u32::try_from(record.pos.as_u64())?,
                    flags: record.flags.raw(),
                    mapq: record.mapq,
                    strand: record.extra(store)?.strand,
                    bases: record.seq(store)?,
                    qualities: &qualities,
                    cigar: record.cigar(store)?,
                },
                reference,
                params,
            )?;
        }
        Ok(())
    }

    fn add(
        &mut self,
        mut alignment: Alignment<'_>,
        reference: &Reference,
        params: &PileupMappingParams,
    ) -> Result<()> {
        if alignment.mapq < params.quality.min_mapq || alignment.mapq == 255 {
            return Ok(());
        }
        if alignment.name.is_empty() || alignment.name == b"*" {
            self.nameless = true;
            return Ok(());
        }
        if !self.seen.insert((
            alignment.group.clone(),
            alignment.name.to_vec(),
            u64::from(alignment.start),
            alignment.flags,
        )) {
            return Ok(());
        }
        if params.unpaired {
            alignment.strand = if alignment.flags & 16 != 0 { Strand::OB } else { Strand::OT };
        }
        if let Some(read) = evidence(&alignment, reference, params)? {
            let mut key = alignment.name.to_vec();
            if params.unpaired || alignment.flags & 1 == 0 {
                key.extend_from_slice(&alignment.start.to_le_bytes());
                key.extend_from_slice(&alignment.flags.to_le_bytes());
            }
            self.named.entry((alignment.group, key)).or_default().push(read);
        }
        Ok(())
    }

    pub fn finish(self) -> Result<Vec<Fragment>> {
        warn_nameless(self.nameless);
        let mut fragments: Vec<Fragment> = self
            .named
            .into_values()
            .map(|reads| {
                ensure!(
                    reads.len() <= 2,
                    "More than two primary alignments share a read name and read group in a segment"
                );
                Ok(Fragment { reads })
            })
            .collect::<Result<_>>()?;
        fragments.sort_by_key(|fragment| {
            fragment
                .reads
                .iter()
                .filter_map(|read| read.sequence.as_ref())
                .filter_map(|sequence| sequence.bases.first())
                .map(|base| base.position.pos)
                .min()
        });
        Ok(fragments)
    }
}

pub(crate) fn warn_nameless(nameless: bool) {
    static WARNED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if nameless && !WARNED.swap(true, std::sync::atomic::Ordering::Relaxed) {
        tracing::warn!("Reads without names cannot contribute molecular haplotype evidence");
    }
}

#[cfg(not(feature = "experimental-seqair"))]
pub(crate) fn fetch_hts(
    params: &crate::sequence::ReaderParams,
    segment: &crate::sequence::Segment,
    mapping: &PileupMappingParams,
) -> Result<Vec<Fragment>> {
    use rust_htslib::bam::{self, Read as _, record::Aux};
    let mut reader = bam::IndexedReader::from_path(params.bam_file.path())?;
    reader.set_reference(params.fasta_file.path())?;
    reader.fetch((segment.contig.as_str(), segment.start, segment.end.saturating_add(1)))?;
    let reference = Reference {
        start: u32::try_from(segment.start)?,
        bases: segment.sequence.iter().copied().map(Base::from).collect(),
    };
    let mut collector = Collector::default();
    for record in reader.records() {
        let record = record.wrap_err("Cannot read BAM for haplotype features")?;
        if !mapping.read_flags.filter_flags(record.flags(), mapping.unpaired) {
            continue;
        }
        // SAFETY: The record owns the BAM allocation for the lifetime of this view.
        let view = unsafe { bam::RecordView::from_raw(record.inner() as *const _) };
        if !mapping.require_tags.allows(&view) {
            continue;
        }
        let group = match record.aux(b"RG") {
            Ok(Aux::String(group)) => group.as_bytes().to_vec(),
            Err(rust_htslib::bam::BamError::AuxTagNotFound) => Vec::new(),
            Ok(_) => bail!("BAM read-group tag is not a string"),
            Err(error) => return Err(error.into()),
        };
        let bases: Vec<_> = record.seq().as_bytes().into_iter().map(Base::from).collect();
        let cigar: Vec<_> = record.raw_cigar().iter().copied().map(CigarOp::from_bam_u32).collect();
        let strand = if mapping.guess_read_orientation {
            crate::call::pileup::hts_utils::infer_strand_from_mismatch_motifs(&record, segment)
        } else {
            Strand::from(seqair_types::BamFlags::from(record.flags()))
        };
        collector.add(
            Alignment {
                name: record.qname(),
                group,
                start: u32::try_from(record.pos())?,
                flags: record.flags(),
                mapq: record.mapq(),
                strand,
                bases: &bases,
                qualities: record.qual(),
                cigar: &cigar,
            },
            &reference,
            mapping,
        )?;
    }
    collector.finish()
}

fn evidence(
    alignment: &Alignment<'_>,
    reference: &Reference,
    params: &PileupMappingParams,
) -> Result<Option<ReadEvidence>> {
    let Alignment { bases, qualities, cigar, strand, .. } = alignment;
    let strand = *strand;
    let read_length = u32::try_from(bases.len())?;
    ensure!(bases.len() == qualities.len(), "Read bases and qualities disagree");
    let reverse = alignment.flags & 16 != 0;
    let Some(mask) = params.read_masking.keep_range(strand, reverse, read_length) else {
        return Ok(None);
    };
    let soft_clipped = cigar.iter().any(|op| op.op_type() == CigarOpType::SoftClip);
    let has_repeat = has_terminal_repeat(bases.len(), |i| bases.get(i));
    let matching_bases = cigar
        .iter()
        .filter(|op| {
            matches!(
                op.op_type(),
                CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch
            )
        })
        .map(|op| op.len())
        .sum();
    let num_indels = cigar
        .iter()
        .filter(|op| matches!(op.op_type(), CigarOpType::Insertion | CigarOpType::Deletion))
        .map(|op| op.len())
        .sum();
    let mut annotations =
        ReadAnnotations { indels: Vec::new(), soft_clipped, noisy: has_repeat || soft_clipped };
    let mut observations = Vec::new();
    let mut sequence =
        SequenceEvidence { bases: Vec::new(), reference_skip: false, indels: Vec::new() };
    let mut pos = alignment.start;
    let mut query = 0u32;
    let mut previous_match = 0;
    let mut mismatches = 0;
    for op in *cigar {
        let len = op.len();
        let kind = op.op_type();
        match kind {
            CigarOpType::Match
            | CigarOpType::SeqMatch
            | CigarOpType::SeqMismatch
            | CigarOpType::Insertion => {
                for offset in 0..len {
                    let genomic = if kind == CigarOpType::Insertion {
                        pos
                    } else {
                        pos.checked_add(offset).ok_or_eyre("Read coordinate overflow")?
                    };
                    if reference.base_at(genomic).is_none() {
                        continue;
                    }
                    let qpos = query.checked_add(offset).ok_or_eyre("Query position overflow")?;
                    let index = usize::try_from(qpos)?;
                    let byte = *qualities.get(index).ok_or_eyre("CIGAR exceeds read qualities")?;
                    let quality = (byte != 255).then_some(byte);
                    let observed = *bases.get(index).ok_or_eyre("CIGAR exceeds read bases")?;
                    if kind != CigarOpType::Insertion
                        && let Some(reference) = reference.base_at(genomic)
                    {
                        let conversion =
                            (reference == Base::C && observed == Base::T && strand != Strand::OB)
                                || (reference == Base::G
                                    && observed == Base::A
                                    && strand != Strand::OT);
                        mismatches += u32::from(
                            reference != Base::Unknown
                                && observed != Base::Unknown
                                && reference != observed
                                && !conversion,
                        );
                    }
                    let accepted = mask.contains(&qpos)
                        && quality.is_some_and(|q| q >= params.quality.min_baseq);
                    let base = if accepted { observed } else { Base::Unknown };
                    let position = if kind == CigarOpType::Insertion {
                        SequencePosition::insertion(genomic, offset)
                    } else {
                        SequencePosition::reference(genomic)
                    };
                    sequence.bases.push(SequenceObservation { position, base, quality });
                    if position.is_reference() && base.known_index().is_some() {
                        observations.push(Observation { pos: genomic, base, quality });
                    }
                }
            }
            CigarOpType::RefSkip => sequence.reference_skip = true,
            CigarOpType::Unknown(code) => bail!("Unsupported BAM CIGAR opcode {code}"),
            _ => {}
        }
        if matches!(kind, CigarOpType::Insertion | CigarOpType::Deletion)
            && reference.base_at(pos).is_some()
        {
            let proposal = if kind == CigarOpType::Insertion {
                let inserted: Vec<_> = (query
                    ..query.checked_add(len).ok_or_eyre("Insertion query overflow")?)
                    .map(|q| {
                        bases
                            .get(q as usize)
                            .copied()
                            .filter(|base| mask.contains(&q) && base.known_index().is_some())
                    })
                    .collect::<Option<_>>()
                    .unwrap_or_default();
                (inserted.len() == len as usize).then_some(Indel {
                    start: pos,
                    inserted,
                    deleted: Vec::new(),
                })
            } else if query.checked_sub(1).is_some_and(|q| mask.contains(&q))
                && mask.contains(&query)
            {
                let deleted = (0..len)
                    .map(|offset| {
                        reference
                            .base_at(pos.checked_add(offset)?)
                            .filter(|base| base.known_index().is_some())
                    })
                    .collect::<Option<Vec<_>>>();
                deleted.map(|deleted| Indel { start: pos, deleted, inserted: Vec::new() })
            } else {
                None
            };
            if let Some(raw) = proposal {
                let indel = raw.normalized(reference)?;
                if len <= 20
                    && (kind == CigarOpType::Deletion
                        || (query..query.saturating_add(len)).all(|q| {
                            qualities
                                .get(q as usize)
                                .is_some_and(|&v| v != 255 && v >= params.quality.min_baseq)
                        }))
                {
                    sequence.indels.push(indel.clone());
                }
                let shift = raw.start - indel.start;
                let anchor_query = (shift <= previous_match && query > shift)
                    .then(|| query - shift - 1)
                    .filter(|q| mask.contains(q));
                let boundary_aligned = anchor_query.is_some();
                let raw_anchor = query.checked_sub(1);
                if raw_anchor.is_some_and(|q| {
                    mask.contains(&q)
                        && qualities
                            .get(q as usize)
                            .is_some_and(|&v| v != 255 && v >= params.quality.min_baseq)
                        && q as usize >= params.indel_end_of_read_cutoff
                        && (q as usize)
                            < bases.len().saturating_sub(params.indel_end_of_read_cutoff)
                }) {
                    let quality = |index: u32| {
                        qualities.get(index as usize).copied().filter(|&q| q != 255).unwrap_or(0)
                    };
                    let normalized_query = anchor_query.map(|q| q + 1);
                    let inserted_qualities: SmallVec<u8, 4> = if indel.inserted.is_empty() {
                        SmallVec::new()
                    } else {
                        normalized_query
                            .into_iter()
                            .flat_map(|start| (start..start.saturating_add(len)).map(quality))
                            .collect()
                    };
                    annotations.indels.push(NormalizedObservation {
                        raw_start: raw.start,
                        indel: indel.clone(),
                        boundary_aligned,
                        observation: IndelObservation {
                            allele: if indel.inserted.is_empty() {
                                IndelAllele::Deletion(indel.deleted.iter().copied().collect())
                            } else {
                                IndelAllele::Insertion(indel.inserted.iter().copied().collect())
                            },
                            strand,
                            reverse,
                            pos_in_read: anchor_query.unwrap_or(0),
                            read_length,
                            mapq: alignment.mapq,
                            base_qual: anchor_query.map(quality).unwrap_or(0),
                            matching_bases,
                            num_indels_in_read: num_indels,
                            insertion_base_quals: inserted_qualities,
                            post_del_base_qual: if indel.deleted.is_empty() {
                                0
                            } else {
                                normalized_query.map(quality).unwrap_or(0)
                            },
                            has_repeat,
                            noisy: annotations.noisy,
                        },
                    });
                }
            }
        }
        if kind.consumes_ref() {
            pos = pos.checked_add(len).ok_or_eyre("Read reference coordinate overflow")?;
        }
        if kind.consumes_query() {
            query = query.checked_add(len).ok_or_eyre("Read query coordinate overflow")?;
        }
        previous_match = if matches!(
            kind,
            CigarOpType::Match | CigarOpType::SeqMatch | CigarOpType::SeqMismatch
        ) {
            len
        } else {
            0
        };
    }
    ensure!(query == read_length, "CIGAR query length does not match sequence length");
    if mismatches > params.indel_max_mismatches {
        annotations.indels.clear();
    }
    Ok((!observations.is_empty()).then_some(ReadEvidence {
        mapq: Some(alignment.mapq),
        strand,
        observations,
        sequence: Some(sequence),
        annotations: Some(Box::new(annotations)),
    }))
}

#[cfg(test)]
pub(super) fn fixture_read(
    reference: &Reference,
    bases: &[u8],
    cigar: &[(CigarOpType, u32)],
) -> Result<ReadEvidence> {
    let bases: Vec<_> = bases.iter().copied().map(Base::from).collect();
    let qualities = vec![40; bases.len()];
    let cigar: Vec<_> = cigar.iter().map(|&(kind, len)| CigarOp::new(kind, len)).collect();
    evidence(
        &Alignment {
            name: b"fixture",
            group: Vec::new(),
            start: reference.start,
            flags: 65,
            mapq: 60,
            strand: Strand::OT,
            bases: &bases,
            qualities: &qualities,
            cigar: &cigar,
        },
        reference,
        &PileupMappingParams::default(),
    )?
    .ok_or_eyre("Empty fixture read")
}

#[cfg(test)]
mod tests {
    use super::*;
    use CigarOpType::{Deletion, Insertion, Match};

    fn reference() -> Reference {
        Reference { start: 1000, bases: b"GCAAAACTG".iter().copied().map(Base::from).collect() }
    }

    #[test]
    fn normalized_boundaries_use_shifted_query_qualities() -> Result<()> {
        let reference = reference();
        for (kind, sequence, remainder, post_quality) in [
            (Insertion, b"GCAAAAACTG".as_slice(), 4, 24),
            (Deletion, b"GCAAACTG".as_slice(), 3, 24),
        ] {
            let bases: Vec<_> = sequence.iter().copied().map(Base::from).collect();
            let qualities =
                (0..bases.len()).map(|i| Ok(22 + u8::try_from(i)?)).collect::<Result<Vec<_>>>()?;
            let cigar =
                [CigarOp::new(Match, 5), CigarOp::new(kind, 1), CigarOp::new(Match, remainder)];
            let alignment = Alignment {
                name: b"x",
                group: Vec::new(),
                start: 1000,
                flags: 65,
                mapq: 60,
                strand: Strand::OT,
                bases: &bases,
                qualities: &qualities,
                cigar: &cigar,
            };
            let read = evidence(&alignment, &reference, &PileupMappingParams::default())?
                .ok_or_eyre("read")?;
            let annotation =
                read.annotations.as_ref().and_then(|a| a.indels.first()).ok_or_eyre("indel")?;
            assert_eq!(annotation.raw_start, 1005);
            assert_eq!(annotation.indel.start, 1002);
            assert!(annotation.boundary_aligned);
            assert_eq!(annotation.observation.pos_in_read, 1);
            assert_eq!(annotation.observation.base_qual, 23);
            if kind == Insertion {
                assert_eq!(annotation.observation.insertion_base_quals.as_slice(), &[post_quality]);
            } else {
                assert_eq!(annotation.observation.post_del_base_qual, post_quality);
            }
            let sequence = read.sequence.as_ref().ok_or_eyre("sequence")?;
            assert_eq!(
                sequence.bases.iter().any(|b| b.position.pos == 1005),
                kind == Insertion,
                "Original CIGAR geometry is retained"
            );
        }
        Ok(())
    }

    #[test]
    fn insertion_outside_the_captured_reference_is_ignored() -> Result<()> {
        let read = fixture_read(&reference(), b"GCAAAACTGT", &[(Match, 9), (Insertion, 1)])?;
        assert!(read.annotations.as_ref().ok_or_eyre("annotations")?.indels.is_empty());
        assert_eq!(read.observations.len(), 9);
        Ok(())
    }

    #[test]
    fn reference_equivalent_shift_before_a_read_has_missing_boundaries() -> Result<()> {
        let reference = reference();
        let bases: Vec<_> = b"AAACTG".iter().copied().map(Base::from).collect();
        let cigar = [CigarOp::new(Match, 1), CigarOp::new(Deletion, 1), CigarOp::new(Match, 5)];
        let read = evidence(
            &Alignment {
                name: b"x",
                group: Vec::new(),
                start: 1002,
                flags: 65,
                mapq: 60,
                strand: Strand::OT,
                bases: &bases,
                qualities: &[40; 6],
                cigar: &cigar,
            },
            &reference,
            &PileupMappingParams::default(),
        )?
        .ok_or_eyre("read")?;
        let annotation =
            read.annotations.as_ref().and_then(|a| a.indels.first()).ok_or_eyre("indel")?;
        assert!(!annotation.boundary_aligned);
        assert_eq!(annotation.observation.base_qual, 0);
        assert_eq!(annotation.indel.start, 1002);
        Ok(())
    }

    #[test]
    fn read_group_separates_molecules_and_repeated_fetches_do_not_duplicate_mates() -> Result<()> {
        let reference = reference();
        let bases = reference.bases.clone();
        let cigar = [CigarOp::new(Match, 9)];
        let mut collector = Collector::default();
        for group in [b"a", b"b"] {
            for flags in [65, 129, 65, 129] {
                collector.add(
                    Alignment {
                        name: b"same-name",
                        group: group.to_vec(),
                        start: 1000,
                        flags,
                        mapq: 60,
                        strand: Strand::OT,
                        bases: &bases,
                        qualities: &[40; 9],
                        cigar: &cigar,
                    },
                    &reference,
                    &PileupMappingParams::default(),
                )?;
            }
        }
        let fragments = collector.finish()?;
        assert_eq!(fragments.len(), 2);
        assert!(fragments.iter().all(|f| f.reads.len() == 2));
        Ok(())
    }
}
