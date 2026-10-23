use super::*;

/// Three reads insert `GA`, three do not; all six are scored, each
/// under its own fragment.
#[test]
fn an_indel_column_is_scored_with_fragment_ids() {
    let mut reads: Vec<TestRead> =
        (0..3).map(|i| inserting_read(format!("ins{i}").as_bytes(), 0, b"GA")).collect();
    reads.extend(
        (0..3).map(|i| haplotype_read(format!("ref{i}").as_bytes(), HAPLOTYPE_START + 5, 50, 16)),
    );
    let scoring = PileupMappingParams { call_indels: true, ..default() };

    let pm = haplotype_column(&reads, &scoring, HAPLOTYPE_ANCHOR);
    let data = pm.indel_data.as_deref().expect("an insertion column");
    let scores = &data.scored.as_ref().expect("scores").haplotypes;
    let inserted = IndelAllele::Insertion([Base::G, Base::A].into_iter().collect());
    assert_eq!(scores.cluster().events().len(), 1);
    assert_eq!(scores.cluster().events()[0].allele, inserted.into());
    let mut ids: Vec<_> = scores.reads.iter().map(|r| r.fragment).collect();
    ids.sort_unstable();
    let mut expected: Vec<_> = reads.iter().map(|r| fragment_of(&r.qname)).collect();
    expected.sort_unstable();
    assert_eq!(ids, expected);
    assert!(ids.iter().all(Option::is_some));
    let (low, high) = scores.verdict.best.haplotypes();
    assert_eq!((low.index(), high.index()), (0, 1));
}

/// An allele below `--min-indel-ao` is not a candidate, so a column with
/// no other has nothing to score.
#[test]
fn an_allele_below_the_minimum_support_is_not_scored() {
    let mut reads = vec![inserting_read(b"ins", 0, b"GA")];
    reads.extend(
        (0..3).map(|i| haplotype_read(format!("ref{i}").as_bytes(), HAPLOTYPE_START + 5, 50, 16)),
    );
    let params = PileupMappingParams { call_indels: true, ..default() };
    let pm = haplotype_column(&reads, &params, HAPLOTYPE_ANCHOR);
    assert!(pm.indel_data.expect("an insertion column").scored.is_none());
}

fn deleted(bases: &[u8]) -> IndelAllele {
    IndelAllele::Deletion(bases.iter().map(|&b| Base::from(b)).collect())
}

/// A column at `pos` whose reads carry `alleles`, each on `total` reads.
fn indel_column(pos: u32, alleles: &[(IndelAllele, u32)]) -> PileupMetrics {
    let mut column = crate::call::phasing::fixtures::column(pos, Base::A, &[], None, &[]);
    let alleles = alleles
        .iter()
        .map(|(allele, total)| IndelAlleleCounts { allele: allele.clone(), total: *total })
        .collect();
    let counts = IndelCounts { alleles };
    column.indel_data = Some(Box::new(IndelData { counts, scored: None }));
    column
}

fn events_of(columns: &[PileupMetrics]) -> Vec<(u32, IndelAllele)> {
    cluster_candidates(columns.iter(), CandidateFloor { min_reads: 2 })
        .into_iter()
        .filter_map(|e| Some((e.anchor, e.allele.as_indel()?.clone())))
        .collect()
}

/// Most reads first, then the earlier anchor, then insertions before
/// deletions; an allele below `--min-indel-ao` is no event, and the
/// fourth is left out.
#[test]
fn a_cluster_takes_its_best_supported_events_in_order() {
    let insertion = IndelAllele::Insertion([Base::A].into_iter().collect());
    let columns = [
        indel_column(100, &[(deleted(b"T"), 5), (insertion.clone(), 5), (deleted(b"TT"), 1)]),
        indel_column(105, &[(deleted(b"G"), 5)]),
        indel_column(110, &[(deleted(b"TT"), 8)]),
    ];
    assert_eq!(
        events_of(&columns),
        [(110, deleted(b"TT")), (100, insertion), (100, deleted(b"T"))]
    );
}

/// A third allele at one anchor could not share a pair with the other
/// two, so the next anchor's event takes its place.
#[test]
fn a_cluster_takes_at_most_two_events_per_anchor() {
    let columns = [
        indel_column(100, &[(deleted(b"T"), 9), (deleted(b"TT"), 8), (deleted(b"TTT"), 7)]),
        indel_column(120, &[(deleted(b"G"), 3)]),
    ];
    assert_eq!(
        events_of(&columns),
        [(100, deleted(b"T")), (100, deleted(b"TT")), (120, deleted(b"G"))]
    );
}

/// A column with `alt` on `ot`/`ob` reads and the reference on
/// `ref_ot`/`ref_ob`, depth their sum.
fn snv_column(
    pos: u32,
    reference: Base,
    alt: Base,
    [ot, ob]: [u32; 2],
    [ref_ot, ref_ob]: [u32; 2],
) -> PileupMetrics {
    let mut column = mismatch_column(pos, reference, alt, ot, ob);
    column.ref_metrics.strand_count.ot = ref_ot;
    column.ref_metrics.strand_count.ob = ref_ob;
    column.ref_metrics.depth = ref_ot + ref_ob;
    column.pos_metrics.depth = ot + ob + ref_ot + ref_ob;
    column
}

fn admitted(column: &PileupMetrics) -> Vec<Base> {
    snv_events(column)
        .iter()
        .filter_map(|e| match e.allele {
            Allele::Base(b) => Some(b),
            Allele::Indel(_) => None,
        })
        .collect()
}

/// Two reads and a het the informative strand allows: a `T` over a `C`
/// counts on OB only, where TAPS converts nothing, so six OT `T`s are a
/// methylated `CpG` and no event, while three OB `T`s among twelve are.
#[test]
fn an_snv_event_needs_two_reads_and_a_possible_het() {
    let (c, t, a) = (Base::C, Base::T, Base::A);
    assert_eq!(admitted(&snv_column(10, c, t, [6, 0], [0, 12])), [] as [Base; 0]);
    assert_eq!(admitted(&snv_column(10, c, t, [0, 3], [6, 9])), [t]);
    assert_eq!(admitted(&snv_column(10, c, a, [1, 0], [6, 6])), [] as [Base; 0], "one read");
    assert_eq!(admitted(&snv_column(10, c, a, [1, 1], [6, 6])), [a]);
    // 3 of 70 is below `CandidateFloor::SNV`.
    assert_eq!(admitted(&snv_column(10, c, a, [2, 1], [33, 34])), [] as [Base; 0]);
    let mut unknown = snv_column(10, c, a, [3, 3], [3, 3]);
    unknown.reference_base = Base::Unknown;
    assert!(snv_events(&unknown).is_empty());
}

fn snv_at(pos: u32) -> SmallVec<Event, 2> {
    [Event { anchor: pos, allele: Allele::Base(Base::G) }].into_iter().collect()
}

/// Indel columns chain while each is within `CLUSTER_GAP` of the last.
#[test]
fn indel_columns_chain_within_the_cluster_gap() {
    let columns = [
        indel_column(100, &[(deleted(b"T"), 5)]),
        indel_column(115, &[(deleted(b"G"), 4)]),
        indel_column(135, &[(deleted(b"G"), 4)]),
        indel_column(200, &[(deleted(b"TT"), 3)]),
    ];
    let chains: Vec<&[usize]> = vec![&[0, 1, 2], &[3]];
    assert_eq!(cluster_chains(&columns, &[0, 1, 2, 3]), chains);
}

/// An SNV near an indel chain is scored with the chain's indel events as
/// context; the rest chain among themselves, cut where a fourth event
/// would not fit, and one that would not fit beside the indels goes with
/// the rest.
#[test]
fn snvs_join_a_nearby_indel_chain_or_chain_on_their_own() {
    let columns: Vec<PileupMetrics> = [90, 100, 110, 125, 135, 140, 145, 150, 400]
        .into_iter()
        .map(|pos| {
            if pos == 100 {
                indel_column(100, &[(deleted(b"T"), 5), (deleted(b"TT"), 4)])
            } else {
                mismatch_column(pos, Base::A, Base::G, 2, 2)
            }
        })
        .collect();
    let snvs: Vec<_> =
        [0usize, 2, 3, 4, 5, 6, 7, 8].into_iter().map(|i| (i, snv_at(columns[i].pos))).collect();
    let clusters = snv_clusters(&columns, &[1], &snvs, CandidateFloor { min_reads: 2 });
    let summary: Vec<(Vec<u32>, usize)> = clusters
        .iter()
        .map(|c| (c.columns.iter().map(|&i| columns[i].pos).collect(), c.events.len()))
        .collect();
    assert_eq!(
        summary,
        [(vec![90, 100], 3), (vec![110, 125, 135], 3), (vec![140, 145, 150], 3), (vec![400], 1),]
    );
    assert!(clusters[0].events[0].allele.as_indel().is_some(), "indel events rank first");
    assert!(clusters[0].events[1].allele.as_indel().is_some());
}

/// A column with `ot` and `ob` reads of `alt` over `reference`.
fn mismatch_column(pos: u32, reference: Base, alt: Base, ot: u32, ob: u32) -> PileupMetrics {
    let mut column = crate::call::phasing::fixtures::column(
        pos,
        reference,
        &[(alt, crate::metrics::AltCall::Uncalled)],
        None,
        &[],
    );
    for alt in &mut column.alts {
        alt.metrics.strand_count.ot = ot;
        alt.metrics.strand_count.ob = ob;
        alt.metrics.depth = ot + ob;
    }
    column
}

/// Five reads of a base TAPS cannot explain mask a column: a `T` over a
/// `C` only on the bottom strand, an `A` over a `G` only on the top, any
/// other base on either; the conversions TAPS makes do not.
#[test]
fn the_window_mask_takes_mismatches_taps_does_not_explain() {
    use Base::{A, C, G, T};
    let columns = [
        mismatch_column(10, C, T, 0, 5),
        mismatch_column(11, C, T, 5, 0),
        mismatch_column(12, G, A, 5, 0),
        mismatch_column(13, G, A, 0, 5),
        mismatch_column(14, A, C, 3, 2),
        mismatch_column(15, A, C, 2, 2),
    ];
    assert_eq!(WindowMask::new(&columns).0, [10, 12, 14]);
}

/// The scores of a 1 bp deletion at [`HAPLOTYPE_ANCHOR`] under `mask`.
fn deletion_scores(reads: &[TestRead], mask: WindowMask) -> Option<HaplotypeScores> {
    let seg = segment_at(u64::from(HAPLOTYPE_START), HAPLOTYPE_SEQ);
    let params = PileupMappingParams { call_indels: true, ..default() };
    let store = store_of(reads, &params);
    let last = HAPLOTYPE_START + u32::try_from(HAPLOTYPE_SEQ.len()).unwrap() - 1;
    let mut engine =
        PileupEngine::new(store.prepare_for_pileup().input, span(HAPLOTYPE_START, last));
    let events: SmallVec<Event, MAX_EVENTS> =
        [Event { anchor: HAPLOTYPE_ANCHOR, allele: deleted(b"T").into() }].into_iter().collect();
    while let Some(col) = engine.pileups() {
        let here = col.pos();
        if here.as_u64() == u64::from(HAPLOTYPE_ANCHOR) {
            let spanning = col.records_overlapping(here, here);
            let scratch = &mut HaplotypeScratch::default();
            let context =
                ClusterContext { store: col.store(), segment: &seg, params: &params, mask };
            return haplotype_scores(&context, spanning, events, scratch);
        }
    }
    None
}

/// Reads ending on an SNV after the anchor whose base is the one after it:
/// the reference charges them a mismatch, the deletion of the SNV's base
/// explains them. Masked, the SNV's position matches any base, and the
/// deletion has nothing left to explain.
#[test]
fn a_masked_position_matches_any_base() {
    let snv = HAPLOTYPE_ANCHOR + 1;
    let offset = (snv - HAPLOTYPE_START) as usize;
    assert_eq!(HAPLOTYPE_SEQ[offset + 1], b'A', "the base after the SNV");
    let mut reads: Vec<TestRead> = (0..6)
        .map(|i| {
            let start = HAPLOTYPE_START + 5 + i;
            let len = (snv - start + 1) as usize;
            haplotype_read(format!("snv{i}").as_bytes(), start, len, 0)
                .with_base_at(len - 1, Base::A)
        })
        .collect();
    reads.extend(
        (0..4).map(|i| haplotype_read(format!("ref{i}").as_bytes(), HAPLOTYPE_START + 5, 50, 16)),
    );
    let event = pair_hmm::EventIdx::new(0).unwrap();

    let unmasked = deletion_scores(&reads, WindowMask::default()).expect("scores");
    assert_eq!(unmasked.verdict.copies_of(event), 1, "{:?}", unmasked.verdict);
    let masked = deletion_scores(&reads, WindowMask(vec![snv])).expect("scores");
    assert_eq!(masked.verdict.copies_of(event), 0, "{:?}", masked.verdict);
}

/// A read at `pos` of [`HAPLOTYPE_SEQ`]'s bases, clipped `lead` bases
/// before and `trail` after `aligned` aligned ones.
fn clipped_read(qname: &[u8], pos: u32, lead: u32, aligned: u32, trail: u32) -> TestRead {
    // The clipped bases are the reference's where it has them; which bases
    // a clip off the segment carries does not matter to its placement.
    let from = usize::try_from(pos.saturating_sub(HAPLOTYPE_START).saturating_sub(lead)).unwrap();
    let len = usize::try_from(lead + aligned + trail).unwrap();
    let mut cigar = vec![CigarOp::new(CigarOpType::Match, aligned)];
    if lead > 0 {
        cigar.insert(0, CigarOp::new(CigarOpType::SoftClip, lead));
    }
    if trail > 0 {
        cigar.push(CigarOp::new(CigarOpType::SoftClip, trail));
    }
    let bases: Vec<Base> =
        HAPLOTYPE_SEQ.iter().cycle().skip(from).take(len).map(|&b| Base::from(b)).collect();
    TestRead {
        cigar,
        quals: vec![40; bases.len()],
        bases,
        ..TestRead::matching(qname, pos, len, Base::A, 0)
    }
}

/// Each read the pair-HMM would score at `pos` of [`HAPLOTYPE_SEQ`] placed
/// at `start`: its placed span, the stored bases it scores, and its name.
fn placed_reads_at(
    reads: &[TestRead],
    start: u32,
    pos: u32,
) -> Vec<(u64, u64, Range<usize>, Vec<u8>)> {
    let seg = segment_at(u64::from(start), HAPLOTYPE_SEQ);
    let params = PileupMappingParams::default();
    let store = store_of(reads, &params);
    let last = start + u32::try_from(HAPLOTYPE_SEQ.len()).unwrap() - 1;
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(start, last));
    let mut placed = Vec::new();
    while let Some(col) = engine.pileups() {
        if col.pos().as_u64() == u64::from(pos) {
            let here = col.pos();
            placed =
                spanning_reads(col.store(), col.records_overlapping(here, here), &seg, &params)
                    .iter()
                    .map(|r| (r.start, r.end, r.bases.clone(), r.rec.qname().to_vec()))
                    .collect();
        }
    }
    placed
}

/// A soft clip the segment's reference reaches is scored with the read,
/// on either end; the last reference base is inside, `end` being inclusive.
#[test]
fn a_clip_that_fits_the_segment_is_placed_with_the_read() {
    let lead = clipped_read(b"lead", HAPLOTYPE_START + 10, 3, 20, 0);
    assert_eq!(
        placed_reads_at(&[lead], HAPLOTYPE_START, HAPLOTYPE_START + 15),
        [(1007, 1029, 0..23, b"lead".to_vec())]
    );
    // Aligned to 1037..=1056, the trailing clip ends on the last base, 1059.
    let trail = clipped_read(b"trail", HAPLOTYPE_START + 37, 0, 20, 3);
    assert_eq!(
        placed_reads_at(&[trail], HAPLOTYPE_START, HAPLOTYPE_START + 45),
        [(1037, 1059, 0..23, b"trail".to_vec())]
    );
}

/// A clip that would run off the segment is dropped, and the read is
/// scored on its aligned bases alone.
#[test]
fn a_clip_outside_the_segment_leaves_the_aligned_bases() {
    let lead = clipped_read(b"lead", HAPLOTYPE_START + 1, 3, 20, 0);
    assert_eq!(
        placed_reads_at(&[lead], HAPLOTYPE_START, HAPLOTYPE_START + 5),
        [(1001, 1020, 3..23, b"lead".to_vec())]
    );
    let trail = clipped_read(b"trail", HAPLOTYPE_START + 38, 0, 20, 3);
    assert_eq!(
        placed_reads_at(&[trail], HAPLOTYPE_START, HAPLOTYPE_START + 45),
        [(1038, 1057, 0..20, b"trail".to_vec())]
    );
}

/// At the contig start a clip can reach before position 0; the read is
/// then placed on its aligned bases, not shifted onto 0.
#[test]
fn a_clip_before_the_contig_start_leaves_the_aligned_bases() {
    let lead = clipped_read(b"lead", 1, 3, 20, 0);
    assert_eq!(placed_reads_at(&[lead], 0, 5), [(1, 20, 3..23, b"lead".to_vec())]);
}

/// One read per fragment, whichever mate reaches further, and a read
/// without a mate is its own fragment.
#[test]
fn a_fragment_is_scored_once_from_the_mate_reaching_further() {
    let mate = |qname: &[u8], pos: u32, len: u32, flags: u16| TestRead {
        flags,
        ..clipped_read(qname, pos, 0, len, 0)
    };
    let reads = [
        mate(b"right", HAPLOTYPE_START, 20, 99),
        mate(b"right", HAPLOTYPE_START + 10, 30, 147),
        mate(b"left", HAPLOTYPE_START + 2, 40, 99),
        mate(b"left", HAPLOTYPE_START + 12, 20, 147),
        mate(b"solo", HAPLOTYPE_START + 5, 20, 0),
    ];
    let mut placed = placed_reads_at(&reads, HAPLOTYPE_START, HAPLOTYPE_START + 15);
    placed.sort_by(|a, b| a.3.cmp(&b.3));
    assert_eq!(
        placed,
        [
            (1002, 1041, 0..40, b"left".to_vec()),
            (1010, 1039, 0..30, b"right".to_vec()),
            (1005, 1024, 0..20, b"solo".to_vec()),
        ]
    );
}
