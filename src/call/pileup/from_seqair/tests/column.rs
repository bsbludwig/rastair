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
    let params = PileupMappingParams { engine: Engine::PairHmm, ..default() };

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
        vec![(deletion.clone(), 1)],
        "the deduplicated column must still report the dropped mate's deletion"
    );

    // With dedup off both mates are their own fragment, and only one of
    // them spans the deletion — so the count is unchanged, not doubled.
    let mut kept = PileupMappingParams { engine: Engine::PairHmm, ..default() };
    kept.variant_calling.keep_overlapping_reads = true;
    assert_eq!(
        alleles_at_anchor(&kept),
        vec![(deletion, 1)],
        "keeping both mates must not turn one fragment's indel into two"
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
    let params = PileupMappingParams { engine: Engine::PairHmm, ..default() };

    let mut extras = TestExtras::new(&params);
    let mut store = RecordStore::<RastairReadExtras>::new();
    store
        .push_fields(
            Pos0::new(START as u32).unwrap(),
            Pos0::new(START as u32 + 9).unwrap(),
            BamFlags::from(99u16),
            MappingQuality::from_byte(60),
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

/// A read aligned to the G of a CpG with its leading base (a methylated T
/// over the C) soft-clipped is rescued into the pileup at the C: with the
/// engine overhang on, the T appears as an OT alt at the CpG-C; with it off,
/// the C position has no column at all.
#[test]
fn rescues_soft_clipped_cpg_partner() {
    // Reference: T T C G T T — CpG is C@2 / G@3.
    let seg = segment(b"TTCGTT");
    let params = PileupMappingParams::default();
    let mut extras = TestExtras::new(&params);

    let mut build_store = || {
        let mut store = RecordStore::<RastairReadExtras>::new();
        // 1S 3M at pos 3: clip base T over ref C@2, aligned G,T,T over 3,4,5.
        // flags 99 = paired/proper/mate-reverse/first → OT.
        store
            .push_fields(
                Pos0::new(3).unwrap(),
                Pos0::new(5).unwrap(),
                BamFlags::from(99u16),
                MappingQuality::from_byte(60),
                3,
                0,
                b"clipped",
                &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
                &[Base::T, Base::G, Base::T, Base::T],
                &[40u8; 4],
                &[],
                Some(tid0()),
                MateFields::default(),
                &mut extras,
            )
            .unwrap();
        store
    };

    let mut metrics_at = |overhang: u32| -> Option<PileupMetrics> {
        let mut engine = PileupEngine::new(build_store().prepare_for_pileup().input, span(0, 5));
        engine.set_soft_clip_overhang(overhang);
        let mut out = None;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(2).unwrap() {
                out = Some(
                    PileupMetrics::from_seqair(
                        &col,
                        seg.clone(),
                        &params,
                        &mut ColumnScratch::default(),
                    )
                    .unwrap(),
                );
            }
        }
        out
    };

    // Overhang off: nothing covers position 2, so no column is produced.
    assert!(metrics_at(0).is_none(), "no rescue without overhang");

    // Overhang on: the clipped T is rescued as an OT alt at the CpG-C.
    let pm = metrics_at(1).expect("CpG-C column emitted via soft-clip rescue");
    assert_eq!(pm.reference_base, Base::C);
    let t = pm.alt(Base::T).expect("rescued T alt present at CpG-C");
    assert_eq!(t.strand_count.ot, 1, "rescued methylation read counted on OT");
    assert_eq!(t.strand_count.ob, 0);
}

#[test]
fn rescues_ob_strand_cpg_partner() {
    // Reference: T T C G T T — CpG is C@2 / G@3.
    let seg = segment(b"TTCGTT");
    let params = PileupMappingParams::default();
    let mut extras = TestExtras::new(&params);

    let mut build_store = || {
        let mut store = RecordStore::<RastairReadExtras>::new();
        // 3M 1S at pos 0: aligned T,T,C over ref 0,1,2; clip base A over ref
        // G@3. flag 83 = paired/proper/reverse/first → OB.
        store
            .push_fields(
                Pos0::new(0).unwrap(),
                Pos0::new(2).unwrap(),
                BamFlags::from(83u16),
                MappingQuality::from_byte(60),
                3,
                0,
                b"clipped",
                &[CigarOp::new(CigarOpType::Match, 3), CigarOp::new(CigarOpType::SoftClip, 1)],
                &[Base::T, Base::T, Base::C, Base::A],
                &[40u8; 4],
                &[],
                Some(tid0()),
                MateFields::default(),
                &mut extras,
            )
            .unwrap();
        store
    };

    let mut metrics_at = |overhang: u32| -> Option<PileupMetrics> {
        let mut engine = PileupEngine::new(build_store().prepare_for_pileup().input, span(0, 5));
        engine.set_soft_clip_overhang(overhang);
        let mut out = None;
        while let Some(col) = engine.pileups() {
            if col.pos() == Pos0::new(3).unwrap() {
                out = Some(
                    PileupMetrics::from_seqair(
                        &col,
                        seg.clone(),
                        &params,
                        &mut ColumnScratch::default(),
                    )
                    .unwrap(),
                );
            }
        }
        out
    };

    assert!(metrics_at(0).is_none(), "no rescue without overhang");

    let pm = metrics_at(1).expect("CpG-G column emitted via soft-clip rescue");
    assert_eq!(pm.reference_base, Base::G);
    let a = pm.alt(Base::A).expect("rescued A alt present at CpG-G");
    assert_eq!(a.strand_count.ob, 1, "rescued methylation read counted on OB");
    assert_eq!(a.strand_count.ot, 0);
}

/// Strand half of the gate: an OT read whose trailing clip lands on the G of
/// a CpG must *not* be rescued — OB methylation evidence (G→A) cannot come
/// from an OT read, so this is a plain end-of-read mismatch.
#[test]
fn does_not_rescue_ot_clip_over_ref_g() {
    // Reference: T T C G T T — G@3 is the CpG-G, before_1 = C@2.
    let seg = segment(b"TTCGTT");
    let params = PileupMappingParams::default();
    let mut extras = TestExtras::new(&params);

    let mut store = RecordStore::<RastairReadExtras>::new();
    // 3M 1S at pos 0, flag 99 → OT. Trailing clip A projects onto ref G@3.
    store
        .push_fields(
            Pos0::new(0).unwrap(),
            Pos0::new(2).unwrap(),
            BamFlags::from(99u16),
            MappingQuality::from_byte(60),
            3,
            0,
            b"clipped",
            &[CigarOp::new(CigarOpType::Match, 3), CigarOp::new(CigarOpType::SoftClip, 1)],
            &[Base::T, Base::T, Base::C, Base::A],
            &[40u8; 4],
            &[],
            Some(tid0()),
            MateFields::default(),
            &mut extras,
        )
        .unwrap();

    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
    engine.set_soft_clip_overhang(1);
    while let Some(col) = engine.pileups() {
        if col.pos() == Pos0::new(3).unwrap() {
            let pm = PileupMetrics::from_seqair(
                &col,
                seg.clone(),
                &params,
                &mut ColumnScratch::default(),
            )
            .unwrap();
            assert!(pm.alt(Base::A).is_none(), "OT clip over ref G must not be rescued");
            assert_eq!(pm.pos_metrics.depth, 0);
        }
    }
}

/// The same clipped base over a non-CpG C (no following G) is *not* rescued:
/// it is a plain end-of-read mismatch, not a methylation partner.
#[test]
fn does_not_rescue_outside_cpg_context() {
    // Reference: T T C A T T — C@2 is followed by A, so not a CpG.
    let seg = segment(b"TTCATT");
    let params = PileupMappingParams::default();
    let mut extras = TestExtras::new(&params);

    let mut store = RecordStore::<RastairReadExtras>::new();
    store
        .push_fields(
            Pos0::new(3).unwrap(),
            Pos0::new(5).unwrap(),
            BamFlags::from(99u16),
            MappingQuality::from_byte(60),
            3,
            0,
            b"clipped",
            &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
            &[Base::T, Base::A, Base::T, Base::T],
            &[40u8; 4],
            &[],
            Some(tid0()),
            MateFields::default(),
            &mut extras,
        )
        .unwrap();

    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
    engine.set_soft_clip_overhang(1);
    while let Some(col) = engine.pileups() {
        if col.pos() == Pos0::new(2).unwrap() {
            let pm = PileupMetrics::from_seqair(
                &col,
                seg.clone(),
                &params,
                &mut ColumnScratch::default(),
            )
            .unwrap();
            // The soft-clip view exists but is gated out: no T alt, no depth.
            assert!(pm.alt(Base::T).is_none(), "non-CpG clip must not be rescued");
            assert_eq!(pm.pos_metrics.depth, 0);
        }
    }
}

/// A clipped base sitting on the right strand at a CpG partner but whose
/// *observed* base is not bisulfite-relevant (here a C→G mismatch over the
/// CpG-C on OT) must not be rescued: it is a fringe SNP/error, not
/// methylation evidence, and rescuing it would feed noise into variant
/// calling.
#[test]
fn does_not_rescue_non_bisulfite_clip_base() {
    // Reference: T T C G T T — CpG is C@2 / G@3.
    let seg = segment(b"TTCGTT");
    let params = PileupMappingParams::default();
    let mut extras = TestExtras::new(&params);

    let mut store = RecordStore::<RastairReadExtras>::new();
    // 1S 3M at pos 3, flag 99 → OT. Clip base G (not T/C) projects onto C@2.
    store
        .push_fields(
            Pos0::new(3).unwrap(),
            Pos0::new(5).unwrap(),
            BamFlags::from(99u16),
            MappingQuality::from_byte(60),
            3,
            0,
            b"clipped",
            &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
            &[Base::G, Base::G, Base::T, Base::T],
            &[40u8; 4],
            &[],
            Some(tid0()),
            MateFields::default(),
            &mut extras,
        )
        .unwrap();

    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
    engine.set_soft_clip_overhang(1);
    while let Some(col) = engine.pileups() {
        if col.pos() == Pos0::new(2).unwrap() {
            let pm = PileupMetrics::from_seqair(
                &col,
                seg.clone(),
                &params,
                &mut ColumnScratch::default(),
            )
            .unwrap();
            assert!(pm.alt(Base::G).is_none(), "non-bisulfite fringe clip must not be rescued");
            assert_eq!(pm.pos_metrics.depth, 0);
        }
    }
}

/// A paired read whose two mates both land on the same CpG-C — one via an
/// aligned base, the other via a rescued soft-clip partner — must be counted
/// once. The engine presents both (depth 2); `from_seqair`'s overlapping-pair
/// dedup collapses them to a single observation, so the rescued fringe base
/// does not double-count the molecule.
#[test]
fn rescued_partner_is_deduped_against_mate() {
    // Reference: T T C G T T — CpG is C@2 / G@3.
    let seg = segment(b"TTCGTT");
    let params = PileupMappingParams::default();
    let mut extras = TestExtras::new(&params);

    let mut store = RecordStore::<RastairReadExtras>::new();
    // Mate A: first in template, forward (flag 99 → OT). 3M at pos 2 covers
    // the CpG-C with an aligned C.
    store
        .push_fields(
            Pos0::new(2).unwrap(),
            Pos0::new(4).unwrap(),
            BamFlags::from(99u16),
            MappingQuality::from_byte(60),
            3,
            0,
            b"pair",
            &[CigarOp::new(CigarOpType::Match, 3)],
            &[Base::C, Base::G, Base::T],
            &[40u8; 3],
            &[],
            Some(tid0()),
            mate_at(3),
            &mut extras,
        )
        .unwrap();
    // Mate B: second in template, reverse (flag 147 → OT). 1S 3M at pos 3,
    // its clipped T projecting back onto the same CpG-C.
    store
        .push_fields(
            Pos0::new(3).unwrap(),
            Pos0::new(5).unwrap(),
            BamFlags::from(147u16),
            MappingQuality::from_byte(60),
            3,
            0,
            b"pair",
            &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
            &[Base::T, Base::G, Base::T, Base::T],
            &[40u8; 4],
            &[],
            Some(tid0()),
            mate_at(2),
            &mut extras,
        )
        .unwrap();

    // `Readers::pileup` links mates after fetching; this test drives the
    // engine directly, so it links by hand.
    let stats = store.link_mates();
    assert_eq!(stats.pairs, 1, "the fixture's mates must link");
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
    engine.set_soft_clip_overhang(1);

    let mut checked = false;
    while let Some(col) = engine.pileups() {
        if col.pos() == Pos0::new(2).unwrap() {
            // The engine presents both the aligned mate and the rescued clip.
            assert_eq!(col.depth(), 2, "both mates present at the CpG-C before dedup");
            let pm = PileupMetrics::from_seqair(
                &col,
                seg.clone(),
                &params,
                &mut ColumnScratch::default(),
            )
            .unwrap();
            assert_eq!(pm.pos_metrics.depth, 1, "rescued partner deduped against its mate");
            checked = true;
        }
    }
    assert!(checked, "CpG-C column must be produced");
}

/// Read-end masking must be applied consistently across both `from_seqair`
/// passes. A rescued soft-clip CpG partner bypasses masking during counting
/// (it is a fringe base by construction); if the *dedup* pass still applies
/// masking to it, the rescued view is dropped from the overlapping-pair
/// collector while its aligned mate is not — so the pair is never detected
/// and the molecule is counted twice. With masking that targets only the
/// clipped mate, the deduped depth must still be 1, not 2.
#[test]
fn rescued_partner_masking_is_consistent_across_passes() {
    // Reference: T T C G T T — CpG is C@2 / G@3.
    let seg = segment(b"TTCGTT");
    // Mask one base from the 3' end of OT reverse reads (the clipped mate B);
    // OT forward (mate A) is left untouched.
    let mut params = PileupMappingParams::default();
    params.variant_calling.read_masking =
        ReadMaskParams::new("0,0,0,1".parse().unwrap(), ReadMaskSetting::default());
    let mut extras = TestExtras::new(&params);

    let mut store = RecordStore::<RastairReadExtras>::new();
    // Mate A: first in template, forward (flag 99 → OT). 3M at pos 2 covers
    // the CpG-C with an aligned C; OT-forward masking is zero so it survives.
    store
        .push_fields(
            Pos0::new(2).unwrap(),
            Pos0::new(4).unwrap(),
            BamFlags::from(99u16),
            MappingQuality::from_byte(60),
            3,
            0,
            b"pair",
            &[CigarOp::new(CigarOpType::Match, 3)],
            &[Base::C, Base::G, Base::T],
            &[40u8; 3],
            &[],
            Some(tid0()),
            mate_at(3),
            &mut extras,
        )
        .unwrap();
    // Mate B: second in template, reverse (flag 147 → OT). 1S 3M at pos 3,
    // its clipped T (read pos 0) projecting onto the same CpG-C. The OT
    // reverse mask rejects read position 0.
    store
        .push_fields(
            Pos0::new(3).unwrap(),
            Pos0::new(5).unwrap(),
            BamFlags::from(147u16),
            MappingQuality::from_byte(60),
            3,
            0,
            b"pair",
            &[CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)],
            &[Base::T, Base::G, Base::T, Base::T],
            &[40u8; 4],
            &[],
            Some(tid0()),
            mate_at(2),
            &mut extras,
        )
        .unwrap();

    // `Readers::pileup` links mates after fetching; this test drives the
    // engine directly, so it links by hand.
    let stats = store.link_mates();
    assert_eq!(stats.pairs, 1, "the fixture's mates must link");
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 5));
    engine.set_soft_clip_overhang(1);

    let mut checked = false;
    while let Some(col) = engine.pileups() {
        if col.pos() == Pos0::new(2).unwrap() {
            let pm = PileupMetrics::from_seqair(
                &col,
                seg.clone(),
                &params,
                &mut ColumnScratch::default(),
            )
            .unwrap();
            assert_eq!(
                pm.pos_metrics.depth, 1,
                "rescued partner deduped against its mate even with read-end masking active"
            );
            checked = true;
        }
    }
    assert!(checked, "CpG-C column must be produced");
}

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

/// With `--rescue-soft-clip-cpg`, a clipped base is projected onto a column
/// outside its own alignment. It still belongs to the same molecule as its
/// mate's aligned base there, so it must dedup against it — which is why
/// the mate-overlap interval is widened by the overhang.
#[test]
fn dedup_matches_the_old_rule_for_a_rescued_soft_clip() {
    // REF has a CpG at 5 (C) / 6 (G). The left mate covers it directly; the
    // right mate's clipped T projects back onto the C.
    let mut clipped = TestRead::matching(b"pair", 6, 3, Base::G, OT_SECOND);
    clipped.cigar =
        vec![CigarOp::new(CigarOpType::SoftClip, 1), CigarOp::new(CigarOpType::Match, 3)];
    clipped.bases = vec![Base::T, Base::G, Base::T, Base::A];
    clipped.quals = vec![40; 4];

    let reads = vec![TestRead::matching(b"pair", 5, 3, Base::C, OT_FIRST), clipped];
    let params = PileupMappingParams { rescue_soft_clip_cpg: true, ..dedup_params() };
    assert_same_as_name_collector(&reads, REF, &params);

    // And the column is genuinely one where dedup has to fire: without the
    // widened interval both mates would be counted at the CpG-C.
    let seg = segment(REF);
    let store = store_of(&reads, &params);
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 19));
    engine.set_soft_clip_overhang(1);
    let mut scratch = ColumnScratch::default();
    let mut depth_at_cpg = None;
    while let Some(col) = engine.pileups() {
        if col.pos() == Pos0::new(5).unwrap() {
            assert_eq!(col.depth(), 2, "engine presents both the aligned base and the clip");
            let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut scratch).unwrap();
            depth_at_cpg = Some(pm.pos_metrics.depth);
        }
    }
    assert_eq!(depth_at_cpg, Some(1), "the rescued clip deduped against its mate");
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
                MappingQuality::from_byte(mapq),
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

/// With `--keep-overlapping-reads` no dedup happens at all, so every
/// filtered observation survives — including both halves of a pair.
#[test]
fn keeping_overlapping_reads_keeps_both_mates() {
    let reads = vec![
        TestRead::matching(b"pair", 2, 8, Base::A, OT_FIRST),
        TestRead::matching(b"pair", 6, 8, Base::A, OT_SECOND),
    ];
    let mut params = dedup_params();
    params.variant_calling.keep_overlapping_reads = true;
    let seg = segment(REF);
    let store = store_of(&reads, &params);
    let mut engine = PileupEngine::new(store.prepare_for_pileup().input, span(0, 19));
    let mut scratch = ColumnScratch::default();
    let mut overlap_depth = None;
    while let Some(col) = engine.pileups() {
        if col.pos() == Pos0::new(7).unwrap() {
            let pm = PileupMetrics::from_seqair(&col, seg.clone(), &params, &mut scratch).unwrap();
            overlap_depth = Some(pm.pos_metrics.depth);
        }
    }
    assert_eq!(overlap_depth, Some(2), "both mates must survive inside the overlap");
}

/// Two mates reading different bases at one column keep one read: the
/// fragment is one observation, whichever mate is wrong.
#[test]
fn a_fragment_whose_mates_disagree_keeps_one_read() {
    let disagree = [
        ref_read(b"pair", 2, 8, OT_FIRST),
        ref_read(b"pair", 6, 8, OT_SECOND).with_base_at(1, Base::G),
    ];
    let (pm, _) = column_with_scratch(&disagree, &dedup_params(), 7);
    assert_eq!(pm.pos_metrics.depth, 1);
}
