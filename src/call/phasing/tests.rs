use super::{
    fixtures::{backing, obs},
    *,
};
use crate::call::phasing::{
    DEFAULT_MIN_SWITCH_QUALITY,
    fixtures::{call, het},
};
use seqair_types::Base::*;

#[test]
fn a_zero_hash_is_no_fragment() {
    assert!(FragmentId::new(0).is_none());
    assert!(FragmentId::new(7).is_some());
}

#[test]
fn all_reference_column_is_not_a_candidate() {
    let column = [obs(1, A, Strand::OT), obs(2, A, Strand::OB), obs(3, A, Strand::OT)];
    assert!(candidate_observations(&column, A).is_none());
}

/// The same column with T on the *unconverted* strand is SNP evidence.
#[test]
fn t_on_the_bottom_strand_at_a_reference_c_is_a_candidate() {
    let column = [
        obs(1, T, Strand::OT),
        obs(2, T, Strand::OB),
        obs(3, T, Strand::OB),
        obs(4, C, Strand::OT),
    ];
    let kept = candidate_observations(&column, C).expect("candidate");
    assert_eq!(kept.len(), column.len(), "the ambiguous observations are kept too");
}

/// Converted reads are kept: outside a CpG they are a second allele, and only
/// the phasing gate knows the context.
#[test]
fn converted_reads_are_a_candidate() {
    let column = [obs(1, T, Strand::OT), obs(2, T, Strand::OT), obs(3, C, Strand::OB)];
    assert!(candidate_observations(&column, C).is_some());

    let bottom = [obs(1, A, Strand::OB), obs(2, A, Strand::OB), obs(3, G, Strand::OT)];
    assert!(candidate_observations(&bottom, G).is_some());
}

#[test]
fn one_converted_read_is_not_a_candidate() {
    let column = [obs(1, T, Strand::OT), obs(2, C, Strand::OT), obs(3, C, Strand::OB)];
    assert!(candidate_observations(&column, C).is_none());
}

#[test]
fn evidence_must_land_on_one_alternative_base() {
    let spread = [obs(1, C, Strand::OT), obs(2, G, Strand::OT), obs(3, A, Strand::OB)];
    assert!(
        candidate_observations(&spread, A).is_none(),
        "one mismatch on each of two bases is noise, not an allele"
    );

    let agreeing = [obs(1, C, Strand::OT), obs(2, C, Strand::OT), obs(3, A, Strand::OB)];
    assert!(candidate_observations(&agreeing, A).is_some());
}

#[test]
fn a_single_mismatch_is_not_enough() {
    let column = [obs(1, G, Strand::OT), obs(2, A, Strand::OB)];
    assert!(candidate_observations(&column, A).is_none());
}

#[test]
fn unknown_bases_are_not_evidence() {
    let column = [obs(1, Unknown, Strand::OT), obs(2, Unknown, Strand::OB)];
    assert!(candidate_observations(&column, A).is_none());
}

#[test]
fn an_unknown_reference_base_is_never_a_candidate() {
    let column = [obs(1, A, Strand::OT), obs(2, A, Strand::OB), obs(3, C, Strand::OT)];
    assert!(candidate_observations(&column, Unknown).is_none());
}

#[test]
fn a_candidate_always_carries_observations() {
    assert!(candidate_observations(&[], C).is_none());
}

#[test]
fn unknown_strand_reads_count_as_evidence() {
    let column = [obs(1, T, Strand::Unknown), obs(2, T, Strand::Unknown)];
    assert!(candidate_observations(&column, C).is_some());
}

// ── phase_segment: anchoring on what the segment emits ─────────────────────

/// Three `A/G` hets at 1000, 1040 and 1080, linked by two templates. Which
/// base each template shows at each site is given per site, in that order.
fn three_hets(template_1: [Base; 3], template_2: [Base; 3]) -> Vec<PileupMetrics> {
    [1000, 1040, 1080]
        .into_iter()
        .zip(template_1)
        .zip(template_2)
        .map(|((pos, first), second)| {
            het(pos, A, G, &[obs(1, first, Strand::OT), obs(2, second, Strand::OT)])
        })
        .collect()
}

fn calls(records: &[PileupMetrics]) -> Vec<Option<PhaseCall>> {
    records.iter().map(|record| call_of(record, SiteKind::Snv)).collect()
}

#[test]
fn a_block_covering_the_whole_segment_is_anchored_on_its_first_site() {
    let mut records = three_hets([A, A, G], [G, G, A]);

    phase_segment(&mut records, |_| true, None, DEFAULT_MIN_SWITCH_QUALITY);

    assert_eq!(
        calls(&records),
        [Some(call(1001, 0, 31)), Some(call(1001, 0, 31)), Some(call(1001, 1, 31))]
    );
}

/// The solver sees every site, the overlap's included, but `PS` has to name a
/// record this segment writes, so the block is anchored on its first emitted
/// site.
#[test]
fn a_block_is_anchored_on_its_first_emitted_site() {
    let mut records = three_hets([A, A, A], [G, G, G]);

    phase_segment(&mut records, |pos| pos >= 1040, None, DEFAULT_MIN_SWITCH_QUALITY);

    assert_eq!(calls(&records), [None, Some(call(1041, 0, 31)), Some(call(1041, 0, 31))]);
}

/// Anchoring on a site where haplotype 1 happens to carry the alt would open
/// the block on `1|0`; the whole block is flipped instead, or the flip would
/// be a switch error.
#[test]
fn a_block_opens_on_the_lower_allele_of_its_anchor() {
    let mut records = three_hets([A, G, A], [G, A, G]);

    phase_segment(&mut records, |pos| pos >= 1040, None, DEFAULT_MIN_SWITCH_QUALITY);

    assert_eq!(calls(&records), [None, Some(call(1041, 0, 31)), Some(call(1041, 1, 31))]);
}

/// Nothing the lone survivor links to is emitted here; the neighbouring
/// segment phases the block.
#[test]
fn a_block_with_one_emitted_site_loses_its_phase() {
    let mut records = three_hets([A, A, A], [G, G, G]);

    phase_segment(&mut records, |pos| pos >= 1080, None, DEFAULT_MIN_SWITCH_QUALITY);

    assert_eq!(calls(&records), [None, None, None]);
}

#[test]
fn a_site_no_fragment_links_stays_unphased() {
    let mut records = three_hets([A, A, A], [G, G, G]);
    if let Some(lonely) = records.get_mut(2) {
        let mut observations = vec![obs(7, A, Strand::OT), obs(8, G, Strand::OT)];
        observations.extend(backing(1080, [A, G]));
        lonely.phase = Some(Box::new(ColumnPhase::Observed(observations.into())));
    }

    phase_segment(&mut records, |_| true, None, DEFAULT_MIN_SWITCH_QUALITY);

    assert_eq!(calls(&records), [Some(call(1001, 0, 31)), Some(call(1001, 0, 31)), None]);
}

/// Nothing reads the observations after the segment is phased.
#[test]
fn phasing_drops_the_observations() {
    let mut records = three_hets([A, A, A], [G, G, G]);

    phase_segment(&mut records, |pos| pos >= 1080, None, DEFAULT_MIN_SWITCH_QUALITY);

    assert!(records.iter().all(|r| r.phase.is_none()));
}

/// An indel column carries no base observations, so phasing has to give it a
/// place for its call.
#[test]
fn a_phased_indel_column_gets_a_phase_of_its_own() {
    use super::fixtures::{column, confirming, het, insertion, scored, with_indels};
    use crate::call::variant_calling::GenotypeTag;
    const NO: f32 = f32::NEG_INFINITY;

    let mut columns = [
        het(1000, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]),
        with_indels(
            column(1040, A, &[], Some(GenotypeTag::HomRef), &[]),
            &[insertion(&[T, T])],
            &[
                [scored(1, [-4.0, -1.0, NO]), scored(2, [-1.0, -4.0, NO])].as_slice(),
                &confirming([0, 1]),
            ]
            .concat(),
        ),
    ];

    phase_segment(&mut columns, |_| true, None, DEFAULT_MIN_SWITCH_QUALITY);

    let [snv, indel] = &columns;
    assert_eq!(call_of(snv, SiteKind::Snv).map(|c| c.first), Some(GtAllele::Ref));
    assert!(call_of(snv, SiteKind::Indel).is_none());
    assert!(call_of(indel, SiteKind::Snv).is_none());
    assert_eq!(
        call_of(indel, SiteKind::Indel).map(|c| (c.set, c.first)),
        Some((call(1001, 0, 31).set, fixtures::allele(1))),
        "fragment 1 carries the SNV's reference and the insertion"
    );
}

/// An indel record can open a block: it anchors `PS`, and the block flips to
/// open on its `0|1`.
#[test]
fn an_indel_site_first_in_a_block_anchors_and_flips_it() {
    use super::fixtures::{column, confirming, het, insertion, scored, with_indels};
    use crate::call::variant_calling::GenotypeTag;
    const NO: f32 = f32::NEG_INFINITY;

    // Fragment 1 carries the insertion and fragment 2 the reference, so
    // fragment 1's `A` and fragment 2's `G` put the SNV's alleles trans to it.
    let mut columns = [
        with_indels(
            column(1000, A, &[], Some(GenotypeTag::HomRef), &[]),
            &[insertion(&[T, T])],
            &[
                [scored(1, [-4.0, -1.0, NO]), scored(2, [-1.0, -4.0, NO])].as_slice(),
                &confirming([0, 1]),
            ]
            .concat(),
        ),
        het(1040, A, G, &[obs(1, A, Strand::OT), obs(2, G, Strand::OT)]),
    ];

    phase_segment(&mut columns, |_| true, None, DEFAULT_MIN_SWITCH_QUALITY);

    let [indel, snv] = &columns;
    let pair = |call: Option<PhaseCall>| call.map(|c| (c.set, c.first));
    assert_eq!(pair(call_of(indel, SiteKind::Indel)), Some((call(1001, 0, 31).set, GtAllele::Ref)));
    assert_eq!(
        pair(call_of(snv, SiteKind::Snv)),
        Some((call(1001, 0, 31).set, fixtures::allele(1)))
    );
}

/// The SNV record is written before the indel record of its column, so at a
/// shared position it is the one that opens on `0|1`.
#[test]
fn at_a_shared_position_the_snv_record_opens_the_block() {
    use super::fixtures::{confirming, insertion, scored, with_indels};
    const NO: f32 = f32::NEG_INFINITY;

    // Fragment 1 shows the SNV's alt and the indel's reference: trans.
    let both = with_indels(
        het(1000, A, G, &[obs(1, G, Strand::OT), obs(2, A, Strand::OT)]),
        &[insertion(&[T, T])],
        &[
            [scored(1, [-1.0, -4.0, NO]), scored(2, [-4.0, -1.0, NO])].as_slice(),
            &confirming([0, 1]),
        ]
        .concat(),
    );
    let mut columns = [both];

    phase_segment(&mut columns, |_| true, None, DEFAULT_MIN_SWITCH_QUALITY);

    let [both] = &columns;
    let first = |kind| call_of(both, kind).map(|c| c.first);
    assert_eq!(first(SiteKind::Snv), Some(GtAllele::Ref));
    assert_eq!(first(SiteKind::Indel), Some(fixtures::allele(1)));
}
