use crate::{
    call::process::PileupMappingParams, metrics::PileupMetrics, sequence::RastairReadExtras,
};
use crate::{
    call::variant_calling::ReadDirection, sequence::ChunkRegion, utils::Base, vcf::SequenceContext,
};
use color_eyre::eyre::{ContextCompat as _, Result, ensure};
use seqair::bam::record_store::RecordStore;
use seqair_types::SmolStr;

mod read;
pub use read::*;

/// A hand-built pileup column: the reads at one position, for tests to turn
/// into [`PileupMetrics`](crate::metrics::PileupMetrics) without a BAM.
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

impl Pileup {
    /// Chromosome name of the segment
    pub fn contig(&self) -> SmolStr {
        self.region.contig.clone()
    }

    /// Reference base right before the variant position
    pub fn ref_before(&self) -> Option<Base> {
        self.context.before_1
    }

    /// Reference base right after the variant position
    pub fn ref_after(&self) -> Option<Base> {
        self.context.after_1
    }
}

impl Pileup {
    /// The column as production builds it: each read becomes a record in a seqair
    /// store, and the engine's column at `pos` goes through
    /// [`ColumnDraft::accumulate`](crate::call::pileup::from_seqair::ColumnDraft::accumulate) and
    /// `finish`; a column no read reaches is
    /// [`ColumnDraft::uncovered_cpg`](crate::call::pileup::from_seqair::ColumnDraft::uncovered_cpg).
    ///
    /// Each read is a fragment of its own, at its `position` in a read of its
    /// length, so overlap dedup never applies; a missing `before_base` or
    /// `after_base` is an `N` there, which no adjacency count takes. The column
    /// takes every read the fixture lists: base quality and read-end masking are
    /// not applied. The column stands alone, so its region entropy is 0.
    pub fn metrics(&self, segment: &crate::sequence::Segment) -> Result<PileupMetrics> {
        self.metrics_with_entropy(segment, 0.0)
    }

    /// [`Self::metrics`], with the region entropy of a column among others.
    pub fn metrics_with_entropy(
        &self,
        segment: &crate::sequence::Segment,
        region_entropy: f64,
    ) -> Result<PileupMetrics> {
        use crate::call::pileup::from_seqair::{ColumnDraft, ColumnScratch};
        use seqair::bam::pileup::PileupEngine;
        use seqair_types::Pos0;

        let mut params = PileupMappingParams::default();
        params.variant_calling.quality.min_baseq = seqair_types::BaseQuality::from_byte(0);
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
            None => {
                ColumnDraft::uncovered_cpg(crate::call::test_helpers::pos0(self.pos), &segment)?
                    .wrap_err("a column without reads has metrics only in a CpG")?
            }
        };
        let metrics = draft.finish(region_entropy)?;
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
            if read.direction == ReadDirection::Reverse {
                flags |= 0x10;
            }
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
                &mut FixtureExtras { strand: read.strand },
            )?;
        }
        Ok(store)
    }
}

/// The strand of the read being pushed; nothing masked.
#[derive(Clone)]
struct FixtureExtras {
    strand: seqair_types::Strand,
}

impl seqair::bam::record_store::CustomizeRecordStore for FixtureExtras {
    type Extra = RastairReadExtras;

    fn compute(
        &mut self,
        rec: &seqair::bam::record_store::SlimRecord,
        _store: &RecordStore<RastairReadExtras>,
    ) -> RastairReadExtras {
        RastairReadExtras { strand: Some(self.strand), mask: Some(0..rec.seq_len) }
    }
}
