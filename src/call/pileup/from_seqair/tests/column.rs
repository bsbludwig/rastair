use super::*;

/// Overlap dedup keeps one read per fragment, and the rule can keep the
/// mate that does *not* span the indel. The fragment's indel evidence must
/// survive that, through the mate.
///
/// Fixture: mates agree on the anchor base, so the rule drops the right
/// mate — which is the one carrying the deletion.
#[test]
fn a_dropped_mate_does_not_take_the_fragments_indel_with_it() {
    // Segment offsets:  0123456789
    let seq = b"AAAACGTTGG";
    let seg = segment(seq);
    let params = PileupMappingParams { call_indels: true, ..default() };

    // Left mate, 6M at 0, no indel. Right mate, 2M2D4M at 2: it deletes the
    // segment's C,G at offsets 4-5, anchored at position 3.
    let reads = [
        TestRead {
            qname: b"pair".to_vec(),
            pos: 0,
            flags: 99,
            bases: vec![Base::A, Base::A, Base::A, Base::A, Base::C, Base::G],
            quals: vec![40; 6],
            mapq: 60,
            cigar: vec![CigarOp::new(CigarOpType::Match, 6)],
            mate_pos: 2,
        },
        TestRead {
            qname: b"pair".to_vec(),
            pos: 2,
            flags: 147,
            bases: vec![Base::A, Base::A, Base::T, Base::T, Base::G, Base::G],
            quals: vec![40; 6],
            mapq: 60,
            cigar: vec![
                CigarOp::new(CigarOpType::Match, 2),
                CigarOp::new(CigarOpType::Deletion, 2),
                CigarOp::new(CigarOpType::Match, 4),
            ],
            mate_pos: 0,
        },
    ];

    let alleles_at_anchor = |params: &PileupMappingParams| {
        let store = store_of(&reads, &params);
        let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 9));
        let mut scratch = ColumnScratch::default();
        let mut out = Vec::new();
        while let Some(col) = engine.pileups() {
            if col.pos().as_u64() != 3 {
                continue;
            }
            let pm = PileupMetrics::from_seqair(&col, seg.clone(), params, &mut scratch).unwrap();
            out = pm
                .indel_data
                .map(|d| d.counts.alleles.iter().map(|a| (a.allele.clone(), a.total)).collect())
                .unwrap_or_default();
        }
        out
    };

    let deletion = IndelAllele::Deletion([Base::C, Base::G].into_iter().collect());
    assert_eq!(
        alleles_at_anchor(&params),
        vec![(deletion, 1)],
        "the deduplicated column must still report the dropped mate's deletion"
    );
}

/// A deletion's REF bases must come from the reference *at the deletion*,
/// which is `genomic_pos + 1 - segment.start` into the segment's sequence.
///
/// Reading it at `view.qpos()`, a read-local offset, instead of the genomic
/// `pos` saturates the subtraction to 0 and returns the segment's opening
/// bases — real-looking sequence from the wrong locus. A segment starting
/// at 0 cannot see this, because there the clamp and the correct answer
/// coincide for a read at offset 0.
#[test]
fn deletion_allele_reads_the_reference_at_the_deletion_site() {
    const START: u64 = 1000;
    // Segment sequence, offset:  0123456789
    //                            AAAACGTTGG
    // A read matches 4 bases from 1000, deletes CG at 1004-1005, matches on.
    // The anchor column is 1003 (the base before the deletion), so the
    // deleted bases are the segment's offsets 4..6 = C,G — and *not* its
    // opening A,A, which is what the read-offset bug produced.
    let seq = b"AAAACGTTGG";
    let seg = segment_at(START, seq);
    let params = PileupMappingParams { call_indels: true, ..default() };

    let mut extras = TestExtras::new(&params);
    let mut store = RecordStore::<RastairReadExtras>::new();
    store
        .push_fields(
            Pos0::new(START as u32).unwrap(),
            Pos0::new(START as u32 + 9).unwrap(),
            BamFlags::from(99u16),
            60,
            8,
            0,
            b"deleter",
            &[
                CigarOp::new(CigarOpType::Match, 4),
                CigarOp::new(CigarOpType::Deletion, 2),
                CigarOp::new(CigarOpType::Match, 4),
            ],
            &[Base::A, Base::A, Base::A, Base::A, Base::T, Base::T, Base::G, Base::G],
            &[40u8; 8],
            &[],
            Some(tid0()),
            MateFields::default(),
            &mut extras,
        )
        .unwrap();

    let mut engine =
        PileupEngine::new(store.prepare_for_pileup().input, span(START as u32, START as u32 + 9));
    let mut allele = None;
    while let Some(col) = engine.pileups() {
        if col.pos().as_u64() != START + 3 {
            continue;
        }
        let metrics =
            PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut ColumnScratch::default())
                .unwrap();
        let data = metrics.indel_data.expect("the anchor column carries the deletion");
        let [counts] = data.counts.alleles.as_ref() else {
            panic!("one deletion allele: {:?}", data.counts.alleles)
        };
        allele = Some(counts.allele.clone());
    }

    assert_eq!(
        allele,
        Some(IndelAllele::Deletion([Base::C, Base::G].into_iter().collect())),
        "deleted bases must be read at the deletion site, not at the segment start"
    );
}

// ── differential: mate links vs. pairing by read name ───────────────────

/// Mates that agree on every base: the later one goes, whichever way round
/// the pair is.
#[test]
fn dedup_matches_the_old_rule_when_mates_agree() {
    for (a_flags, b_flags) in [(OT_FIRST, OT_SECOND), (OB_SECOND, OB_FIRST)] {
        let reads = vec![
            TestRead::matching(b"pair", 2, 8, Base::A, a_flags),
            TestRead::matching(b"pair", 6, 8, Base::A, b_flags),
        ];
        assert_same_as_name_collector(&reads, REF, &dedup_params());
    }
}

/// Mates that disagree at the overlapped base. Which one survives depends
/// on which is second in the template — the case the single-pass rewrite
/// had to get right without seeing the pair as a unit.
#[test]
fn dedup_matches_the_old_rule_when_mates_disagree() {
    for (a_flags, b_flags) in [(OT_FIRST, OT_SECOND), (OT_SECOND, OT_FIRST)] {
        let reads = vec![
            TestRead::matching(b"pair", 2, 8, Base::A, a_flags),
            TestRead::matching(b"pair", 6, 8, Base::A, b_flags).with_base_at(0, Base::G),
        ];
        assert_same_as_name_collector(&reads, REF, &dedup_params());
    }
}

/// A mate that fails base quality, mapping quality, or read-end masking at
/// the shared column never forms a pair there, so the other one survives
/// even though both cover the position.
#[test]
fn dedup_matches_the_old_rule_when_one_mate_is_filtered() {
    let params = dedup_params();
    let low_baseq = vec![
        TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
        TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND).with_qual(2),
    ];
    assert_same_as_name_collector(&low_baseq, REF, &params);

    let low_mapq = vec![
        TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST).with_mapq(0),
        TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND),
    ];
    assert_same_as_name_collector(&low_mapq, REF, &params);

    // Mask the first two bases of every read: at the overlap's left edge
    // one mate is masked out while the other is not.
    let mut masked = dedup_params();
    masked.variant_calling.read_masking =
        ReadMaskParams::new("2,0,2,0".parse().unwrap(), "2,0,2,0".parse().unwrap());
    let reads = vec![
        TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
        TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND),
    ];
    assert_same_as_name_collector(&reads, REF, &masked);
}

/// A deletion in one mate: at the deleted positions it has no base, so no
/// pair forms and the other mate stands alone.
#[test]
fn dedup_matches_the_old_rule_across_a_deletion() {
    let mut with_del = TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND);
    with_del.cigar = vec![
        CigarOp::new(CigarOpType::Match, 2),
        CigarOp::new(CigarOpType::Deletion, 2),
        CigarOp::new(CigarOpType::Match, 6),
    ];
    with_del.bases = vec![Base::A; 8];
    with_del.quals = vec![40; 8];
    let reads = vec![TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST), with_del];
    assert_same_as_name_collector(&reads, REF, &dedup_params());
}

/// An insertion in one mate — the anchor base is still a normal
/// observation, and the pair must resolve there like any other column.
#[test]
fn dedup_matches_the_old_rule_across_an_insertion() {
    let mut with_ins = TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND);
    with_ins.cigar = vec![
        CigarOp::new(CigarOpType::Match, 2),
        CigarOp::new(CigarOpType::Insertion, 2),
        CigarOp::new(CigarOpType::Match, 6),
    ];
    let reads = vec![TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST), with_ins];
    assert_same_as_name_collector(&reads, REF, &dedup_params());
}

/// The exact edges of the overlap: the pair meets on one base only, at the
/// last position of the left mate. This is where the half-open/inclusive
/// mix-up in the overlap interval showed up on real data.
#[test]
fn dedup_matches_the_old_rule_at_the_overlap_boundaries() {
    // Left mate covers 2..=9, right mate 9..=16: they share exactly base 9.
    let reads = vec![
        TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
        TestRead::matching(b"pair", 9, 8, Base::A, OT_SECOND),
    ];
    assert_same_as_name_collector(&reads, REF, &dedup_params());

    // And one base further apart: no shared position at all.
    let disjoint = vec![
        TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
        TestRead::matching(b"pair", 10, 8, Base::A, OT_SECOND),
    ];
    assert_same_as_name_collector(&disjoint, REF, &dedup_params());
}

/// Several pairs at one column, interleaved, some agreeing and some not —
/// the pending-drop bookkeeping has to keep them apart.
#[test]
fn dedup_matches_the_old_rule_for_interleaved_pairs() {
    let reads = vec![
        TestRead::matching(b"p1", 2, 8, Base::A, OT_FIRST),
        TestRead::matching(b"p2", 3, 8, Base::A, OB_FIRST),
        TestRead::matching(b"p3", 4, 8, Base::A, OT_SECOND),
        TestRead::matching(b"p1", 6, 8, Base::A, OT_SECOND).with_base_at(0, Base::G),
        TestRead::matching(b"p2", 7, 8, Base::A, OB_SECOND),
        TestRead::matching(b"p3", 8, 8, Base::A, OT_FIRST).with_base_at(0, Base::T),
    ];
    assert_same_as_name_collector(&reads, REF, &dedup_params());
}

/// `--max-coverage` counts reads that pass the record filter. Two MAPQ 0
/// reads arrive first at the column; with a cap of 2 they must not take
/// the two slots, or the column comes out empty.
#[test]
fn coverage_cap_counts_filtered_reads_only() {
    let reads = [(b"junk1", 0u8), (b"junk2", 0), (b"good1", 60), (b"good2", 60), (b"good3", 60)];
    let params = dedup_params();
    let filter = RastairRecordFilter {
        read_flags: params.read_flags.clone(),
        min_mapq: params.quality.min_mapq,
        ..Default::default()
    };
    let mut capped = DepthCap::with_cap(filter, NonZeroU32::new(2));
    let mut store = RecordStore::<RastairReadExtras>::new();
    for (name, mapq) in reads {
        store
            .push_fields(
                Pos0::new(2).unwrap(),
                Pos0::new(9).unwrap(),
                BamFlags::from(OT_FIRST),
                mapq,
                8,
                0,
                name,
                &[CigarOp::new(CigarOpType::Match, 8)],
                &[Base::A; 8],
                &[40u8; 8],
                &[],
                Some(tid0()),
                mate_at(500),
                &mut capped,
            )
            .unwrap();
    }
    let seg = segment(REF);
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 9));
    let mut scratch = ColumnScratch::default();
    let mut depth = None;
    while let Some(col) = engine.pileups() {
        if col.pos().as_u64() != 5 {
            continue;
        }
        let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut scratch).unwrap();
        depth = Some(pm.pos_metrics.depth);
    }
    assert_eq!(depth, Some(2), "the cap applies to reads that pass the MAPQ filter");
}

/// Two mates reading different bases at one column keep one read: the
/// fragment is one observation, whichever mate is wrong.
#[test]
fn a_fragment_whose_mates_disagree_keeps_one_read() {
    let disagree = [
        ref_read(b"pair", 2, 8, OT_FIRST),
        ref_read(b"pair", 6, 8, OT_SECOND).with_base_at(1, Base::G),
    ];
    let (pm, _) = phase_column(&disagree, &dedup_params(), 7);
    assert_eq!(pm.pos_metrics.depth, 1);
}
