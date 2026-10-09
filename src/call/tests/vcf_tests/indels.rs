//! Indel records: `QUAL` from the allele's margin,
//! `GT` and `GC` from the pair-HMM's best pair.

use crate::call::{
    pair_hmm::fixtures::*,
    phasing::fixtures::{column, with_indels},
    pileup::indels::IndelAllele,
    tests::utils::*,
};
use seqair_types::Base::{self, *};

fn deletion(bases: &[Base]) -> IndelAllele {
    IndelAllele::Deletion(bases.iter().copied().collect())
}

/// A column at reference `A` whose indel calls the pair-HMM makes from `reads`
/// scored against `alleles`, as the pipeline would.
fn called(alleles: &[IndelAllele], reads: &[([f32; 3], usize)]) -> crate::metrics::PileupMetrics {
    with_indels(column(1000, A, &[], None, &[]), alleles, &rows(reads))
}

/// Each data line as `(REF, ALT, QUAL, FILTER, GT, GC)`.
fn records(pileup: crate::metrics::PileupMetrics, filters: RecordFilters) -> Vec<[String; 6]> {
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
                format(&record, "GC"),
            ]
        })
        .collect()
}

#[test]
fn a_pair_hmm_het_writes_its_margin_as_qual_and_its_pair_quality_as_gc() {
    let pileup = called(&[deletion(&[T])], &[(REF_READ, 2), (ALT1_READ, 2)]);
    let data = pileup.indel_data.as_deref().expect("indel data");
    let scored = data.scored.as_ref().expect("scores");
    let verdict = &scored.haplotypes.verdict;
    let [call] = scored.records.calls() else { panic!("one call: {:?}", scored.records) };
    let (qual, gc) = (call.quality.as_int(), verdict.quality.as_int());
    assert!(qual >= 10 && gc > 0, "{qual} {gc}");

    assert_eq!(
        records(pileup, RecordFilters::variants()),
        [["AT", "A", &qual.to_string(), "PASS", "0/1", &gc.to_string()].map(str::to_owned)]
    );
}

/// Two called alleles are one multi-allelic record, `1/2`.
#[test]
fn a_pair_hmm_compound_het_is_one_record() {
    let pileup = called(&[deletion(&[T]), deletion(&[T, T])], &[(ALT1_READ, 3), (ALT2_READ, 3)]);
    let [[reference, alts, _, filter, gt, _]] = &records(pileup, RecordFilters::variants())[..]
    else {
        panic!("one record")
    };
    assert_eq!((reference.as_str(), filter.as_str(), gt.as_str()), ("ATT", "PASS", "1/2"));
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

    let [[reference, alt, qual, filter, gt, _]] = &records(rejected(), RecordFilters::all())[..]
    else {
        panic!("one record")
    };
    assert_eq!(
        [reference, alt, filter, gt].map(String::as_str),
        ["AT", "A", "low_indel_qual", "0/0"]
    );
    assert!(qual.parse::<f32>().is_ok_and(|qual| qual < 10.0), "{qual}");
}
