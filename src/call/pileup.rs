#[cfg(test)]
use crate::{
    call::process::PileupMappingParams, metrics::PileupMetrics, sequence::RastairReadExtras,
};
#[cfg(test)]
use crate::{sequence::ChunkRegion, utils::Base, vcf::SequenceContext};
#[cfg(test)]
use color_eyre::eyre::{ContextCompat as _, Result, ensure};
#[cfg(test)]
use seqair::bam::record_store::RecordStore;
#[cfg(test)]
use seqair_types::{SmallVec, SmolStr};

pub mod indels;
#[cfg(test)]
mod read;
#[cfg(test)]
pub use read::*;
pub(crate) mod from_seqair;

/// A hand-built pileup column: the reads at one position, for tests to turn
/// into [`PileupMetrics`](crate::metrics::PileupMetrics) without a BAM.
#[cfg(test)]
#[derive(Debug, Clone)]
pub struct Pileup {
    /// Region of the chunk this pileup belongs to
    pub region: std::sync::Arc<ChunkRegion>,
    /// Sequence context around the position in the reference
    pub context: SequenceContext,
    /// Position in the sequence, 0-based
    pub pos: u32,
    // if this becomes slow, consider boxing or Rc-ing this
    pub reads: SimpleReads,
    /// Reference base at this position
    pub reference_base: Base,
}

#[cfg(test)]
impl Pileup {
    /// Chromosome name of the segment
    pub fn contig(&self) -> SmolStr {
        self.region.contig.clone()
    }

    /// Position in the segment sequence, 0-based
    pub fn idx(&self) -> usize {
        self.region.pos_to_idx(self.pos).expect("valid position")
    }

    /// Reference base right before the variant position
    pub fn ref_before(&self) -> Option<Base> {
        self.context.before_1
    }

    /// Reference base right after the variant position
    pub fn ref_after(&self) -> Option<Base> {
        self.context.after_1
    }

    pub fn alleles(&self) -> SmallVec<Base, 4> {
        let mut res = SmallVec::new();
        res.push(self.reference_base);
        self.reads.iter().map(|b| b.base).fold(res, |mut acc, base| {
            if !acc.contains(&base) {
                acc.push(base);
            }
            acc
        })
    }
}

#[cfg(test)]
impl Pileup {
    /// The column as production builds it: each read becomes a record in a seqair
    /// store, and the engine's column at `pos` goes through
    /// [`ColumnDraft::accumulate`](from_seqair::ColumnDraft::accumulate) and
    /// `finish`; a column no read reaches is
    /// [`ColumnDraft::uncovered_cpg`](from_seqair::ColumnDraft::uncovered_cpg).
    ///
    /// Each read is a fragment of its own, at its `position` in a read of its
    /// length, so overlap dedup never applies; a missing `before_base` or
    /// `after_base` is an `N` there, which no adjacency count takes. The column
    /// takes every read the fixture lists: base quality and read-end masking are
    /// not applied.
    pub fn metrics(&self, segment: &crate::sequence::Segment) -> Result<PileupMetrics> {
        use from_seqair::{ColumnDraft, ColumnScratch};
        use seqair::bam::pileup::PileupEngine;
        use seqair_types::Pos0;

        let mut params = PileupMappingParams::default();
        params.variant_calling.quality.min_baseq = 0;
        let segment = std::rc::Rc::new(segment.clone());
        let store = self.store()?;
        let pos = Pos0::new(self.pos).wrap_err("pileup position out of range")?;
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, (pos..=pos).into());
        let draft = match engine.pileups() {
            Some(column) => ColumnDraft::accumulate(
                &column,
                std::rc::Rc::clone(&segment),
                &params,
                &mut ColumnScratch::default(),
            )?,
            None => ColumnDraft::uncovered_cpg(u64::from(self.pos), &segment)?
                .wrap_err("a column without reads has metrics only in a CpG")?,
        };
        let metrics = draft.finish()?;
        ensure!(
            metrics.reference_base == self.reference_base && metrics.context == self.context,
            "the pileup's reference ({} in {:?}) is not the segment's at {} ({} in {:?})",
            self.reference_base,
            self.context,
            self.pos,
            metrics.reference_base,
            metrics.context,
        );
        Ok(metrics)
    }

    fn store(&self) -> Result<RecordStore<RastairReadExtras>> {
        use seqair::bam::{
            BamHeader, MateFields,
            cigar::{CigarOp, CigarOpType},
        };
        use seqair::reader::ResolveTid;
        use seqair_types::{BamFlags, Pos0};

        let tid = 0u32.resolve_tid(&BamHeader::synthetic(1, u64::from(u32::MAX)))?;
        let mut extras = FixtureExtras { strand: seqair_types::Strand::Unknown };
        let mut store = RecordStore::<RastairReadExtras>::new();
        for (i, read) in self.reads.iter().enumerate() {
            let PositionInRead { pos: qpos, read_length } = read.position;
            let offset = |i: Option<u32>| usize::try_from(i?).ok();
            let mut bases = vec![Base::Unknown; usize::try_from(read_length)?];
            let mut place = |at: Option<usize>, base: Base| -> Result<()> {
                let slot = at.and_then(|i| bases.get_mut(i)).wrap_err("base outside the read")?;
                *slot = base;
                Ok(())
            };
            place(offset(Some(qpos)), read.base)?;
            if let Some(before) = read.before_base {
                place(offset(qpos.checked_sub(1)), before)?;
            }
            if let Some(after) = read.after_base {
                place(offset(qpos.checked_add(1)), after)?;
            }
            let start = self.pos.checked_sub(qpos).wrap_err("read starts before the contig")?;
            let end = (start + read_length).checked_sub(1).wrap_err("empty read")?;
            let mut flags = if read.second { 0x81 } else { 0x41 };
            if read.reverse {
                flags |= 0x10;
            }
            extras.strand = read.strand;
            store.push_fields(
                Pos0::new(start).wrap_err("read start out of range")?,
                Pos0::new(end).wrap_err("read end out of range")?,
                BamFlags::from(flags),
                read.mapq,
                read.matching_bases,
                read.indels,
                format!("read{i}").as_bytes(),
                &[CigarOp::new(CigarOpType::Match, read_length)],
                &bases,
                &vec![read.qual; bases.len()],
                &[],
                Some(tid),
                MateFields { tid: None, pos: None, template_len: None },
                &mut extras,
            )?;
        }
        Ok(store)
    }
}

/// The strand of the read being pushed; nothing masked.
#[cfg(test)]
#[derive(Clone)]
struct FixtureExtras {
    strand: seqair_types::Strand,
}

#[cfg(test)]
impl seqair::bam::record_store::CustomizeRecordStore for FixtureExtras {
    type Extra = RastairReadExtras;

    fn compute(
        &mut self,
        rec: &seqair::bam::record_store::SlimRecord,
        _store: &RecordStore<RastairReadExtras>,
    ) -> RastairReadExtras {
        RastairReadExtras { strand: self.strand, mask: 0..rec.seq_len }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        sequence::{ChunkRegion, Region, Segment},
        utils::default,
    };
    use insta::assert_debug_snapshot;
    use seqair_types::Strand;

    #[test]
    fn test_alleles_in_order() {
        let bases = SimpleReads(
            vec![
                SimpleRead { base: Base::A, strand: Strand::OT, ..default() },
                SimpleRead {
                    base: Base::C,
                    strand: Strand::OB,
                    reverse: true,
                    second: false,
                    ..default()
                },
                SimpleRead { base: Base::A, strand: Strand::OT, ..default() },
            ]
            .into(),
        );

        let segment = Segment {
            range: std::sync::Arc::new(ChunkRegion {
                region: Region { contig: "chr19".into(), start: 1000, end: 1100 },
                last_position: 2000,
                overlap_start: 0,
                overlap_end: 0,
            }),
            sequence: vec![],
            overlap_start: 0,
            overlap_end: 0,
        };
        let variant_candidate = Pileup {
            region: segment.range.clone(),
            context: SequenceContext::default(),
            pos: 1002, // Corresponds to index in the segment
            reads: bases,
            reference_base: Base::T, // Assume T is the reference base at this position
        };

        let alleles = variant_candidate.alleles();
        assert_eq!(alleles[0], Base::T, "Reference base should be first");

        assert_debug_snapshot!(alleles, @r"
        [
            T,
            A,
            C,
        ]
        ");
    }
}
