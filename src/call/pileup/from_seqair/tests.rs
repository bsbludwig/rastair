//! Fixtures for the column and cluster tests: reads pushed into a seqair
//! store and walked by the pileup engine, as `Readers::pileup` does.

use super::{candidates::*, cluster::*, column::*};
use crate::call::pileup::indels::IndelData;
use crate::call::pileup::indels::{IndelAllele, IndelAlleleCounts, IndelCounts};
use crate::call::variant_calling::ReadMaskParams;
use crate::sequence::RastairRecordFilter;
use crate::sequence::{ChunkRegion, Region};
use crate::utils::default;
use crate::{
    call::{
        pair_hmm::{self, Allele, Event, HaplotypeScores, HaplotypeScratch, MAX_EVENTS},
        phasing::FragmentId,
        process::PileupMappingParams,
    },
    metrics::PileupMetrics,
    sequence::{RastairReadExtras, Segment},
};
use color_eyre::eyre::Result;
use seqair::bam::DepthCap;
use seqair::bam::cigar::{CigarOp, CigarOpType};
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::{CustomizeRecordStore, FilterRawFields, SlimRecord};
use seqair::bam::{BamHeader, MateFields};
use seqair::bam::{
    pileup::{AlignmentView, PileupColumn},
    record_store::RecordStore,
};
use seqair::reader::{ResolveTid, Tid};
use seqair_types::{BamFlags, Pos0};
use seqair_types::{Base, SmallVec, Strand};
use std::num::NonZeroU32;
use std::{ops::Range, rc::Rc};

mod cluster;
mod column;

/// The only reference the test reads are on; mates are only linked when
/// both records name the same one.
fn tid0() -> Tid {
    0u32.resolve_tid(&BamHeader::synthetic(1, 1_000)).unwrap()
}

/// A mate on the test reference at `pos`; a negative `pos` means no mate position.
fn mate_at(pos: i32) -> MateFields {
    MateFields {
        tid: Some(tid0()),
        pos: u32::try_from(pos).ok().and_then(Pos0::new),
        template_len: None,
    }
}

/// The parts of `RastairRecordFilter` these tests need: the MAPQ filter,
/// strand from the flags, and the read-end mask window.
#[derive(Default, Clone)]
struct TestExtras {
    masking: ReadMaskParams,
    min_mapq: u8,
}

impl TestExtras {
    fn new(params: &PileupMappingParams) -> Self {
        Self { masking: params.read_masking.clone(), min_mapq: params.quality.min_mapq }
    }
}

impl CustomizeRecordStore for TestExtras {
    type Extra = RastairReadExtras;

    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        fields.mapq >= self.min_mapq
    }

    fn compute(
        &mut self,
        rec: &SlimRecord,
        _store: &RecordStore<RastairReadExtras>,
    ) -> RastairReadExtras {
        let strand = Strand::from(rec.flags);
        RastairReadExtras {
            strand,
            mask: self
                .masking
                .keep_range(strand, rec.flags.is_reverse(), rec.seq_len)
                .unwrap_or(0..0),
        }
    }
}

/// A read to push into a test store: everything the overlap dedup looks at.
#[derive(Clone, Debug)]
struct TestRead {
    qname: Vec<u8>,
    pos: u32,
    flags: u16,
    bases: Vec<Base>,
    quals: Vec<u8>,
    mapq: u8,
    cigar: Vec<CigarOp>,
    mate_pos: i32,
}

impl TestRead {
    /// A plain `<len>M` read whose every base is `base`.
    fn matching(qname: &[u8], pos: u32, len: usize, base: Base, flags: u16) -> Self {
        Self {
            qname: qname.to_vec(),
            pos,
            flags,
            bases: vec![base; len],
            quals: vec![40; len],
            mapq: 60,
            cigar: vec![CigarOp::new(CigarOpType::Match, len as u32)],
            mate_pos: -1,
        }
    }

    fn with_base_at(mut self, offset: usize, base: Base) -> Self {
        if let Some(slot) = self.bases.get_mut(offset) {
            *slot = base;
        }
        self
    }

    fn with_qual(mut self, qual: u8) -> Self {
        self.quals = vec![qual; self.bases.len()];
        self
    }

    fn with_mapq(mut self, mapq: u8) -> Self {
        self.mapq = mapq;
        self
    }

    fn ref_span(&self) -> u32 {
        self.cigar
            .iter()
            .filter(|op| matches!(op.op_type(), CigarOpType::Match | CigarOpType::Deletion))
            .map(|op| op.len())
            .sum()
    }
}

/// Push `reads` as one mate pair per qname, wire up the mate fields, and
/// link — what `Readers::pileup` does for a real fetch.
fn store_of(reads: &[TestRead], params: &PileupMappingParams) -> RecordStore<RastairReadExtras> {
    let mut extras = TestExtras::new(params);
    let mut store = RecordStore::<RastairReadExtras>::new();
    for read in reads {
        let mate_pos = reads
            .iter()
            .find(|other| other.qname == read.qname && other.pos != read.pos)
            .map_or(read.mate_pos, |other| other.pos as i32);
        let end = read.pos + read.ref_span().max(1) - 1;
        store
            .push_fields(
                Pos0::new(read.pos).unwrap(),
                Pos0::new(end).unwrap(),
                BamFlags::from(read.flags),
                read.mapq,
                read.bases.len() as u32,
                0,
                &read.qname,
                &read.cigar,
                &read.bases,
                &read.quals,
                &[],
                Some(tid0()),
                mate_at(mate_pos),
                &mut extras,
            )
            .unwrap();
    }
    // Linking (and sorting) is `prepare_for_pileup`'s job at the point the
    // engine is built, so reads may be listed in any order.
    store
}

/// The overlapping-pair rule by read name: group the column's surviving
/// alignments by qname and, on the second one, drop whichever
/// `resolve_pair` chose. An independent oracle — it reaches the answer by
/// matching names, which is exactly what the mate links do not do.
fn kept_by_name_collector(
    column: &PileupColumn<'_, RastairReadExtras>,
    params: &PileupMappingParams,
) -> Vec<usize> {
    let passing: Vec<usize> = column
        .alignments()
        .enumerate()
        .filter_map(|(idx, view)| {
            let baseq = view.qual()?.get()?;
            let _base = view.base()?;
            view.qpos()?;
            if baseq < params.quality.min_baseq {
                return None;
            }
            if !passes_read_masking(&view) {
                return None;
            }
            Some(idx)
        })
        .collect();

    let alignments: Vec<AlignmentView<'_, '_, RastairReadExtras>> = column.alignments().collect();
    let mut first_by_name: Vec<(&[u8], usize)> = Vec::new();
    let mut removed: Vec<usize> = Vec::new();
    for &idx in &passing {
        let view = &alignments[idx];
        let name = view.qname();
        match first_by_name.iter().find(|(seen, _)| *seen == name) {
            None => first_by_name.push((name, idx)),
            Some(&(_, other)) => {
                let this_base = view.base();
                let other_base = alignments[other].base();
                if this_base == other_base || view.flags.is_second_in_template() {
                    removed.push(idx);
                } else {
                    removed.push(other);
                }
            }
        }
    }
    passing.into_iter().filter(|idx| !removed.contains(idx)).collect()
}

/// Run both implementations over every column of `reads` and require the
/// same surviving observations: total depth, and per allele the depth and
/// the OT/OB split.
fn assert_same_as_name_collector(reads: &[TestRead], seq: &[u8], params: &PileupMappingParams) {
    let seg = segment(seq);
    let store = store_of(reads, &params);
    let last = seq.len().saturating_sub(1) as u32;
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, last));

    let mut scratch = ColumnScratch::default();
    let mut columns = 0;
    while let Some(col) = engine.pileups() {
        let pos = col.pos().as_u64() as usize;
        let expected = kept_by_name_collector(&col, params);

        let alignments: Vec<AlignmentView<'_, '_, RastairReadExtras>> = col.alignments().collect();
        let mut expected_depth = 0u32;
        let mut expected_by_base: Vec<(Base, u32, u32)> = Vec::new();
        for &idx in &expected {
            let view = &alignments[idx];
            let Some(base) = view.base() else { continue };
            expected_depth += 1;
            let strand = view.extra().strand;
            let slot = match expected_by_base.iter_mut().find(|(b, _, _)| *b == base) {
                Some(slot) => slot,
                None => {
                    expected_by_base.push((base, 0, 0));
                    expected_by_base.last_mut().unwrap()
                }
            };
            match strand {
                Strand::OT => slot.1 += 1,
                Strand::OB => slot.2 += 1,
                Strand::Unknown => {}
            }
        }

        let pm = PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap();
        assert_eq!(
            pm.pos_metrics.depth, expected_depth,
            "pos {pos}: depth differs from the name-collector rule"
        );
        for (base, ot, ob) in expected_by_base {
            let metrics = if base == pm.reference_base {
                &pm.ref_metrics
            } else {
                pm.alts
                    .iter()
                    .find(|alt| alt.base == base)
                    .map(|alt| &alt.metrics)
                    .unwrap_or_else(|| panic!("pos {pos}: no alt for {base}"))
            };
            assert_eq!(metrics.depth, ot + ob, "pos {pos}, {base}: allele depth");
            assert_eq!(metrics.strand_count.ot, ot, "pos {pos}, {base}: OT count");
            assert_eq!(metrics.strand_count.ob, ob, "pos {pos}, {base}: OB count");
        }
        columns += 1;
    }
    assert!(columns > 0, "no columns produced");
}

/// The closed position span seqair's pileup engine and readers take.
fn span(first: u32, last: u32) -> core::range::RangeInclusive<Pos0> {
    (Pos0::new(first).unwrap()..=Pos0::new(last).unwrap()).into()
}

fn segment(seq: &[u8]) -> Rc<Segment> {
    segment_at(0, seq)
}

/// A segment placed at `start` on the contig. Everything a segment resolves
/// against the reference goes through `pos - start`, so a fixture at 0
/// cannot tell a genomic position from a read offset — see
/// [`deletion_allele_reads_the_reference_at_the_deletion_site`].
fn segment_at(start: u64, seq: &[u8]) -> Rc<Segment> {
    let end = start + seq.len().saturating_sub(1) as u64;
    Rc::new(Segment {
        range: std::sync::Arc::new(ChunkRegion {
            region: Region { contig: "chr1".into(), start, end },
            last_position: start + seq.len() as u64,
            overlap_start: 0,
            overlap_end: 0,
        }),
        sequence: seq.to_vec(),
        overlap_start: 0,
        overlap_end: 0,
    })
}

const OT_FIRST: u16 = 99; // paired, proper, mate reverse, first in template
const OT_SECOND: u16 = 147; // paired, proper, reverse, second in template
const OB_FIRST: u16 = 83; // paired, proper, reverse, first in template
const OB_SECOND: u16 = 163; // paired, proper, mate reverse, second in template

/// `ACGTACGTACGTACGTACGT`, long enough for two 8bp reads to overlap.
const REF: &[u8] = b"ACGTACGTACGTACGTACGT";

fn dedup_params() -> PileupMappingParams {
    PileupMappingParams::default()
}

// ── column reads ────────────────────────────────────────────────────────

/// A read whose bases are exactly the reference over `[pos, pos + len)`, so
/// that `with_base_at` introduces the only mismatch in the column.
fn ref_read(qname: &[u8], pos: u32, len: usize, flags: u16) -> TestRead {
    let start = pos as usize;
    let mut read = TestRead::matching(qname, pos, len, Base::A, flags);
    read.bases = REF
        .get(start..start + len)
        .expect("read lies inside the reference")
        .iter()
        .map(|&b| Base::from(b))
        .collect();
    read
}

/// Build the column at `pos`, handing back the scratch too so a test can see
/// the counters that live for a whole region.
fn phase_column(
    reads: &[TestRead],
    params: &PileupMappingParams,
    pos: u32,
) -> (PileupMetrics, ColumnScratch) {
    let seg = segment(REF);
    let store = store_of(reads, params);
    let last = u32::try_from(REF.len() - 1).expect("a short reference");
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, last));
    let mut scratch = ColumnScratch::default();
    let mut metrics = None;
    while let Some(col) = engine.pileups() {
        if col.pos() == Pos0::new(pos).unwrap() {
            metrics =
                Some(PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap());
        }
    }
    (metrics.expect("column at the requested position"), scratch)
}

/// Non-repetitive reference for the pair-HMM tests, placed off zero.
const HAPLOTYPE_START: u32 = 1000;
const HAPLOTYPE_SEQ: &[u8] = b"ACGTTAGCATCGGATCCATGACTGATCGGCTAGCTTACGGATCCTAGGCTAAGCTTGACT";
/// The column the pair-HMM tests' insertions sit at.
const HAPLOTYPE_ANCHOR: u32 = HAPLOTYPE_START + 29;

/// `<len>M` at `pos` copying [`HAPLOTYPE_SEQ`], its own fragment.
fn haplotype_read(qname: &[u8], pos: u32, len: usize, flags: u16) -> TestRead {
    let from = (pos - HAPLOTYPE_START) as usize;
    let mut read = TestRead::matching(qname, pos, len, Base::A, flags);
    read.bases = HAPLOTYPE_SEQ[from..from + len].iter().map(|&b| Base::from(b)).collect();
    read
}

/// A read over offsets 5..55 of [`HAPLOTYPE_SEQ`] with `inserted` after the
/// anchor, 25 aligned bases either side.
fn inserting_read(qname: &[u8], flags: u16, inserted: &[u8]) -> TestRead {
    let bases: Vec<Base> = HAPLOTYPE_SEQ[5..30]
        .iter()
        .chain(inserted)
        .chain(&HAPLOTYPE_SEQ[30..55])
        .map(|&b| Base::from(b))
        .collect();
    TestRead {
        qname: qname.to_vec(),
        pos: HAPLOTYPE_START + 5,
        flags,
        quals: vec![40; bases.len()],
        bases,
        mapq: 60,
        cigar: vec![
            CigarOp::new(CigarOpType::Match, 25),
            CigarOp::new(CigarOpType::Insertion, u32::try_from(inserted.len()).unwrap()),
            CigarOp::new(CigarOpType::Match, 25),
        ],
        mate_pos: -1,
    }
}

/// The column at `pos` over [`HAPLOTYPE_SEQ`].
fn haplotype_column(reads: &[TestRead], params: &PileupMappingParams, pos: u32) -> PileupMetrics {
    let seg = segment_at(u64::from(HAPLOTYPE_START), HAPLOTYPE_SEQ);
    let store = store_of(reads, params);
    let last = HAPLOTYPE_START + u32::try_from(HAPLOTYPE_SEQ.len()).unwrap() - 1;
    let mut engine =
        PileupEngine::new(store.prepare_for_pileup().input, span(HAPLOTYPE_START, last));
    let mut scratch = ColumnScratch::default();
    let mut metrics = None;
    while let Some(col) = engine.pileups() {
        if col.pos().as_u64() == u64::from(pos) {
            metrics =
                Some(PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap());
        }
    }
    metrics.expect("the column")
}

fn fragment_of(qname: &[u8]) -> Option<FragmentId> {
    FragmentId::new(seqair::bam::record_store::qname_hash(qname))
}

impl PileupMetrics {
    /// Accumulate and finish one column in one go.
    ///
    /// Production splits these two halves so a column can be rejected between
    /// them; the tests want every column, and this is the form they ask for.
    pub(crate) fn from_seqair(
        column: &PileupColumn<'_, RastairReadExtras>,
        segment: Rc<Segment>,
        params: &PileupMappingParams,
        scratch: &mut ColumnScratch,
    ) -> Result<PileupMetrics> {
        let draft = ColumnDraft::accumulate(column, Rc::clone(&segment), params, scratch)?;
        let mut metrics = draft.finish()?;
        let pos = column.pos();
        let spanning = column.records_overlapping(pos, pos);
        let context = ClusterContext::new(column.store(), &segment, params, &[]);
        let column = std::slice::from_mut(&mut metrics);
        let events = cluster_candidates(column.iter(), params.indel_floor);
        score_cluster(column, &[0], events, ClusterPass::Indels, spanning, &context, scratch);
        flush_clusters(column, scratch);
        Ok(metrics)
    }
}
