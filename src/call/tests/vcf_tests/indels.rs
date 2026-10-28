//! Indel records: `QUAL` from the allele's margin, `GT` from the called pair,
//! `GQ` and `PL` from the record's genotype likelihoods.

use crate::call::{
    pair_hmm::fixtures::*,
    phasing::{
        PhaseCall, RecordPhases,
        fixtures::{allele, call},
    },
    pileup::indels::IndelAllele,
    tests::{
        fixtures::{column, with_indels},
        utils::*,
    },
};
use seqair_types::Base::{self, *};

fn deletion(bases: &[Base]) -> IndelAllele {
    IndelAllele::Deletion(bases.iter().copied().collect())
}

/// A column at reference `A`, covered by `reads`, whose indel calls the
/// pair-HMM makes from them scored against `alleles`, as the pipeline would.
fn called(alleles: &[IndelAllele], reads: &[([f32; 3], usize)]) -> crate::metrics::PileupMetrics {
    let mut site = column(1000, A, &[], None, &[]);
    site.pos_metrics.depth = reads.iter().map(|&(_, n)| u32::try_from(n).unwrap_or(u32::MAX)).sum();
    with_indels(site, alleles, &rows(reads))
}

/// Each data line as `(REF, ALT, QUAL, FILTER, GT, GQ, PL)`.
fn records(pileup: crate::metrics::PileupMetrics, filters: RecordFilters) -> Vec<[String; 7]> {
    let format = |record: &VcfRecord, key| record.format.get(key).cloned().unwrap_or_default();
    metrics_to_vcf(&[pileup], filters)
        .expect("a VCF")
        .into_iter()
        .map(|record| {
            [
                record.r#ref.clone(),
                record.alt.join(","),
                record.qual.clone(),
                record.filters.join(";"),
                format(&record, "GT"),
                format(&record, "GQ"),
                format(&record, "PL"),
            ]
        })
        .collect()
}

#[test]
fn a_pair_hmm_het_writes_its_margin_as_qual_and_its_likelihoods_as_gq_and_pl() {
    let pileup = called(&[deletion(&[T])], &[(REF_READ, 2), (ALT1_READ, 2)]);
    let data = pileup.indel_data.as_deref().expect("indel data");
    let scored = data.scored.as_ref().expect("scores");
    let calls: Vec<_> = scored.records.scored_calls().collect();
    let [(call, score)] = calls[..] else { panic!("one call: {:?}", scored.records) };
    let qual = score.quality.as_int();
    let gq = score.genotype_quality(call.genotype).as_int();
    let pl = score.likelihoods.phred();
    assert!(qual >= 10 && gq > 0, "{qual} {gq}");
    // The written genotype is the likeliest, so its `PL` is 0 and the
    // runner-up's is at least its `GQ`.
    let [hom_ref, het, hom_alt] = pl[..] else { panic!("three genotypes: {pl:?}") };
    assert_eq!(het, 0, "{pl:?}");
    assert!(hom_ref.min(hom_alt) >= gq, "{pl:?} {gq}");

    let pl = format!("{hom_ref},{het},{hom_alt}");
    assert_eq!(
        records(pileup, RecordFilters::variants()),
        [["AT", "A", &qual.to_string(), "PASS", "0/1", &gq.to_string(), &pl].map(str::to_owned)]
    );
}

/// Two called alleles are one multi-allelic record, `1/2`.
#[test]
fn a_pair_hmm_compound_het_is_one_record() {
    let pileup = called(&[deletion(&[T]), deletion(&[T, T])], &[(ALT1_READ, 3), (ALT2_READ, 3)]);
    let [[reference, alts, _, filter, gt, _, pl]] = &records(pileup, RecordFilters::variants())[..]
    else {
        panic!("one record")
    };
    assert_eq!((reference.as_str(), filter.as_str(), gt.as_str()), ("ATT", "PASS", "1/2"));
    // Six genotypes in `GL` order, `1/2` (the fifth) the likeliest.
    let pl: Vec<i32> = pl.split(',').map(|p| p.parse().expect("an integer")).collect();
    assert_eq!(pl.len(), 6, "{pl:?}");
    assert_eq!(pl.get(4), Some(&0), "{pl:?}");
    let mut alts: Vec<_> = alts.split(',').collect();
    alts.sort_unstable();
    assert_eq!(alts, ["A", "AT"]);
}

/// An allele the pair-HMM does not call is written only under `--all`, as
/// `0/0` with its low `QUAL` and the `low_indel_qual` FILTER.
#[test]
fn a_rejected_pair_hmm_allele_is_filtered_and_only_written_under_all() {
    let rejected = || called(&[deletion(&[T])], &[(REF_READ, 8), (ALT1_READ, 1)]);
    assert!(records(rejected(), RecordFilters::variants()).is_empty());

    let [[reference, alt, qual, filter, gt, _, _]] = &records(rejected(), RecordFilters::all())[..]
    else {
        panic!("one record")
    };
    assert_eq!(
        [reference, alt, filter, gt].map(String::as_str),
        ["AT", "A", "low_indel_qual", "0/0"]
    );
    assert!(qual.parse::<f32>().is_ok_and(|qual| qual < 10.0), "{qual}");
}

/// `GT` and `PS` of each data line, `PS` with `PQ` after a colon.
fn phases(pileup: crate::metrics::PileupMetrics) -> Vec<(String, Option<String>)> {
    metrics_to_vcf(&[pileup], RecordFilters::all())
        .expect("a VCF")
        .into_iter()
        .map(|record| {
            let ps = record
                .format
                .get("PS")
                .map(|ps| format!("{ps}:{}", record.format.get("PQ").map_or("", String::as_str)));
            (record.format.get("GT").cloned().unwrap_or_default(), ps)
        })
        .collect()
}

fn with_indel_phase(
    mut pileup: crate::metrics::PileupMetrics,
    first: u8,
) -> crate::metrics::PileupMetrics {
    let phase = PhaseCall { first: allele(first), ..call(901, 0, 31) };
    pileup.phase = Some(Box::new(RecordPhases { snv: None, indel: Some(phase) }));
    pileup
}

/// Phasing made a site of the column's one passing het indel record, so the
/// phase goes on that record, and not on a filtered one beside it.
#[test]
fn a_phased_pair_hmm_het_carries_its_phase_and_a_rejected_neighbour_does_not() {
    let pileup = called(&[deletion(&[T]), deletion(&[T, T])], &[(REF_READ, 3), (ALT1_READ, 3)]);
    assert_eq!(
        phases(with_indel_phase(pileup, 1)),
        [("1|0".to_owned(), Some("901:31".to_owned())), ("0/0".to_owned(), None)]
    );
}

#[test]
fn a_phased_pair_hmm_compound_het_is_one_phased_record() {
    let pileup = called(&[deletion(&[T]), deletion(&[T, T])], &[(ALT1_READ, 3), (ALT2_READ, 3)]);
    assert_eq!(
        phases(with_indel_phase(pileup, 2)),
        [("2|1".to_owned(), Some("901:31".to_owned()))]
    );
}
