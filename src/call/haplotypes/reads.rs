use super::{
    Fragment, Indel, Observation, ReadEvidence, Reference, SequenceEvidence, SequenceObservation,
    SequencePosition,
};
use crate::{call::process::PileupMappingParams, sequence::RastairReadExtras};
use color_eyre::eyre::{OptionExt as _, Result, bail, ensure};
use seqair::bam::{
    cigar::CigarOpType,
    record_store::{RecordStore, SlimRecord},
};
use seqair_types::{Base, Strand};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Default)]
pub(crate) struct Collector {
    named: BTreeMap<(Vec<u8>, Vec<u8>), Vec<ReadEvidence>>,
    seen: BTreeSet<(Vec<u8>, Vec<u8>, u64, u16)>,
    pub nameless: bool,
}

impl Collector {
    pub fn capture(
        &mut self,
        store: &RecordStore<RastairReadExtras>,
        reference: &Reference,
        params: &PileupMappingParams,
    ) -> Result<()> {
        for record in store.records() {
            if record.mapq < params.quality.min_mapq || record.mapq == 255 {
                continue;
            }
            let name = record.qname(store)?;
            if name.is_empty() || name == b"*" {
                self.nameless = true;
                continue;
            }
            let group = match record.aux(store)?.get::<&[u8]>(*b"RG") {
                Ok(group) => group.to_vec(),
                Err(seqair::bam::aux::GetAuxError::TagNotFound { .. }) => Vec::new(),
                Err(error) => return Err(error.into()),
            };
            // Byte-budget subdivisions fetch the same boundary alignment again.
            if !self.seen.insert((
                group.clone(),
                name.to_vec(),
                record.pos.as_u64(),
                record.flags.raw(),
            )) {
                continue;
            }
            if let Some(read) = evidence(record, store, reference, params)? {
                let mut key = name.to_vec();
                if params.unpaired || !record.flags.is_paired() {
                    key.extend_from_slice(&record.pos.as_u64().to_le_bytes());
                    key.extend_from_slice(&record.flags.raw().to_le_bytes());
                }
                self.named.entry((group, key)).or_default().push(read);
            }
        }
        Ok(())
    }

    pub fn finish(self) -> Result<Vec<Fragment>> {
        self.named.into_values().map(|reads| {
            ensure!(reads.len() <= 2, "More than two primary alignments share a read name and read group in a phase segment");
            Ok(Fragment { reads })
        }).collect()
    }
}

fn evidence(
    record: &SlimRecord,
    store: &RecordStore<RastairReadExtras>,
    reference: &Reference,
    params: &PileupMappingParams,
) -> Result<Option<ReadEvidence>> {
    let extra = record.extra(store)?;
    let strand = if params.unpaired {
        if record.flags.is_reverse() { Strand::OB } else { Strand::OT }
    } else {
        extra.strand
    };
    let Some(mask) =
        params.read_masking.keep_range(strand, record.flags.is_reverse(), record.seq_len)
    else {
        return Ok(None);
    };
    let bases = record.seq(store)?;
    let qualities = record.qual(store)?;
    let mut observations = Vec::new();
    let mut sequence =
        SequenceEvidence { bases: Vec::new(), reference_skip: false, indels: Vec::new() };
    let mut pos = record.pos.as_u64();
    let mut query = 0u32;
    for op in record.cigar(store)? {
        let len = op.len();
        match op.op_type() {
            CigarOpType::Match
            | CigarOpType::SeqMatch
            | CigarOpType::SeqMismatch
            | CigarOpType::Insertion => {
                for offset in 0..len {
                    let genomic = if op.op_type() == CigarOpType::Insertion {
                        pos
                    } else {
                        pos + u64::from(offset)
                    };
                    let genomic = u32::try_from(genomic)?;
                    if reference.base_at(genomic).is_none() {
                        continue;
                    }
                    let qpos = query.checked_add(offset).ok_or_eyre("Query position overflow")?;
                    let index = usize::try_from(qpos)?;
                    let quality =
                        qualities.get(index).ok_or_eyre("CIGAR exceeds quality sequence")?.get();
                    let base = *bases.get(index).ok_or_eyre("CIGAR exceeds base sequence")?;
                    let accepted = mask.contains(&qpos)
                        && quality.is_some_and(|q| q >= params.quality.min_baseq);
                    let base = if accepted { base } else { Base::Unknown };
                    let position = if op.op_type() == CigarOpType::Insertion {
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
        if len <= 20
            && (op.op_type() == CigarOpType::Insertion || op.op_type() == CigarOpType::Deletion)
        {
            let start = u32::try_from(pos)?;
            let proposal = if op.op_type() == CigarOpType::Insertion {
                let inserted: Vec<_> = sequence
                    .bases
                    .iter()
                    .filter(|base| base.position.pos == start && !base.position.is_reference())
                    .map(|base| base.base)
                    .collect();
                (inserted.len() == len as usize
                    && inserted.iter().all(|base| base.known_index().is_some()))
                .then_some(Indel { start, inserted, deleted: Vec::new() })
            } else if query.checked_sub(1).is_some_and(|query| mask.contains(&query))
                && mask.contains(&query)
            {
                let deleted: Option<Vec<_>> = (0..len)
                    .map(|offset| {
                        reference
                            .base_at(start.checked_add(offset)?)
                            .filter(|base| base.known_index().is_some())
                    })
                    .collect();
                deleted.map(|deleted| Indel { start, deleted, inserted: Vec::new() })
            } else {
                None
            };
            if let Some(indel) = proposal {
                sequence.indels.push(indel.normalized(reference)?);
            }
        }
        if op.op_type().consumes_ref() {
            pos += u64::from(len);
        }
        if op.op_type().consumes_query() {
            query = query.checked_add(len).ok_or_eyre("Query position overflow")?;
        }
    }
    ensure!(query == record.seq_len, "CIGAR query length does not match sequence length");
    Ok((!observations.is_empty()).then_some(ReadEvidence {
        mapq: Some(record.mapq),
        strand,
        observations,
        sequence: Some(sequence),
    }))
}
