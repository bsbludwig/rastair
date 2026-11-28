use super::*;
use Base::{A, C, G, T};
use bcf::Read as _;
use hegel::{TestCase, generators as gs};
use rust_htslib::bam::record::{Cigar, CigarString};

fn bases() -> impl gs::PrintableGenerator<Base> {
    gs::sampled_from([A, C, G, T].as_slice())
}

fn strands() -> impl gs::PrintableGenerator<Strand> {
    gs::sampled_from([Strand::OT, Strand::OB].as_slice())
}

fn either_haplotype() -> impl gs::PrintableGenerator<Haplotype> {
    gs::sampled_from([Haplotype::One, Haplotype::Two].as_slice())
}

#[hegel::composite]
fn sequences(tc: &TestCase, max_len: usize) -> Vec<u8> {
    let len = tc.draw(gs::integers::<usize>().min_value(0).max_value(max_len));
    (0..len).map(|_| tc.draw(gs::sampled_from(b"ACGT".as_slice()))).collect()
}

/// A molecule of either haplotype, methylated or not, read on either strand:
/// whatever base it shows, it never speaks for the other haplotype, whether or
/// not the read's strand is known.
#[hegel::test]
fn a_read_base_never_votes_for_the_other_haplotype(tc: TestCase) {
    let first = tc.draw(bases());
    let second = tc.draw(gs::sampled_from(
        [A, C, G, T].into_iter().filter(|base| *base != first).collect::<Vec<_>>(),
    ));
    let alleles = [first, second];
    let truth = tc.draw(either_haplotype());
    let strand = tc.draw(strands());
    let (prev, next) = (tc.draw(bases()), tc.draw(bases()));
    let methylated = tc.draw(gs::booleans());

    let allele = alleles[truth.index()];
    let observed = match (strand, allele) {
        (Strand::OT, C) if methylated && next == G => T,
        (Strand::OB, G) if methylated && prev == C => A,
        _ => allele,
    };

    for known_strand in [Some(strand), None] {
        let vote = snv_vote(alleles, known_strand, observed, Some(prev), Some(next));
        assert!(vote.is_none_or(|vote| vote == truth), "{vote:?} for a read of {truth:?}");
    }
    if !methylated {
        // An unconverted base is silent only where it could be the other
        // allele converted.
        let other = alleles[1 - truth.index()];
        let confounded = reads_as(other, observed, Some(strand), prev, next);
        let vote = snv_vote(alleles, Some(strand), observed, Some(prev), Some(next));
        assert_eq!(vote, (!confounded).then_some(truth));
    }
}

/// The cases the conversion rule exists for, at a `C/T` het.
#[test]
fn a_converted_base_votes_only_outside_a_cpg() {
    let alleles = [C, T];
    let vote = |strand, observed, next| snv_vote(alleles, Some(strand), observed, Some(A), next);

    // OT, followed by G: a `T` is the alt or the methylated reference.
    assert_eq!(vote(Strand::OT, T, Some(G)), None);
    // Not followed by G on this molecule: no CpG, so the `T` is the alt.
    assert_eq!(vote(Strand::OT, T, Some(A)), Some(Haplotype::Two));
    // The read ends here, so the next base may be a G.
    assert_eq!(vote(Strand::OT, T, None), None);
    assert_eq!(vote(Strand::OT, T, Some(Base::Unknown)), None);
    // Nor does a missing neighbour make a `T` speak for the `C` of an `A/C`
    // het: it is a converted `C` only if a G follows.
    assert_eq!(snv_vote([A, C], Some(Strand::OT), T, Some(A), None), None);
    assert_eq!(snv_vote([A, C], Some(Strand::OT), T, Some(A), Some(G)), Some(Haplotype::Two));
    assert_eq!(snv_vote([A, C], Some(Strand::OT), T, Some(A), Some(A)), None);
    // A base no conversion makes votes whatever follows.
    assert_eq!(vote(Strand::OT, C, None), Some(Haplotype::One));
    // A `C` is never a converted `T`.
    assert_eq!(vote(Strand::OT, C, Some(G)), Some(Haplotype::One));
    // OB converts no `C`.
    assert_eq!(vote(Strand::OB, T, Some(G)), Some(Haplotype::Two));
    // A third base fits neither.
    assert_eq!(vote(Strand::OT, A, Some(G)), None);

    // The mirror image at a `G/A` het: OB, preceded by C.
    let vote = |strand, prev| snv_vote([G, A], Some(strand), A, prev, Some(T));
    assert_eq!(vote(Strand::OB, Some(C)), None);
    assert_eq!(vote(Strand::OB, Some(T)), Some(Haplotype::Two));
    assert_eq!(vote(Strand::OT, Some(C)), Some(Haplotype::Two));
}

/// A read of one haplotype, of any length and converted anywhere its strand
/// converts, never matches the other haplotype alone. The flank repeats a
/// short unit, where an indel is hardest to tell from the reference.
#[hegel::test]
fn a_read_after_an_indel_never_votes_for_the_other_haplotype(tc: TestCase) {
    let unit = tc.draw(sequences(3));
    let copies = tc.draw(gs::integers::<usize>().min_value(0).max_value(8));
    let mut flank: Vec<u8> = unit.iter().copied().cycle().take(unit.len() * copies).collect();
    flank.extend(tc.draw(sequences(12)));
    let tails = [tc.draw(sequences(6)), tc.draw(sequences(6))];
    let haplotypes = tails.map(|tail| [tail, flank.clone()].concat());
    let truth = tc.draw(either_haplotype());
    let strand = tc.draw(strands());

    let molecule = &haplotypes[truth.index()];
    let len = tc.draw(gs::integers::<usize>().min_value(0).max_value(molecule.len()));
    let read: Vec<u8> = molecule
        .iter()
        .take(len)
        .map(|&base| match (strand, base) {
            (Strand::OT, b'C') if tc.draw(gs::booleans()) => b'T',
            (Strand::OB, b'G') if tc.draw(gs::booleans()) => b'A',
            _ => base,
        })
        .collect();

    for known_strand in [Some(strand), None] {
        let vote = indel_vote(&haplotypes, known_strand, &read);
        assert!(
            vote.is_none_or(|(vote, _)| vote == truth),
            "{vote:?} for a read of {truth:?}: {}",
            String::from_utf8_lossy(&read)
        );
        if let Some((_, decisive)) = vote {
            let [one, two] = &haplotypes;
            assert_ne!(one.get(decisive), two.get(decisive), "not where the haplotypes differ");
            assert!(decisive < read.len(), "the read does not reach the deciding base");
        }
    }
}

/// One `A` deleted from `AAAA`, anchored on the base before: a read has to
/// leave the run to say which haplotype it is.
#[test]
fn a_read_ending_inside_a_repeat_casts_no_indel_vote() {
    let haplotypes = [b"AAAACGTT".to_vec(), b"AAACGTT".to_vec()];
    let vote = |read: &[u8]| indel_vote(&haplotypes, Some(Strand::OB), read);

    assert_eq!(vote(b"AAA"), None, "three A fit both");
    assert_eq!(vote(b"AAAA"), Some((Haplotype::One, 3)));
    assert_eq!(vote(b"AAAC"), Some((Haplotype::Two, 3)));
    assert_eq!(vote(b"AAACGTT"), Some((Haplotype::Two, 3)));
    assert_eq!(vote(b"AAACGAT"), None, "a mismatch to both haplotypes");
    assert_eq!(vote(b""), None);
}

/// An inserted `C` against a reference `T` cannot be told from a converted
/// reference on OT, and can on OB.
#[test]
fn an_indel_that_differs_by_a_conversion_is_silent_on_the_converting_strand() {
    let haplotypes = [b"CAG".to_vec(), b"TAG".to_vec()];
    assert_eq!(indel_vote(&haplotypes, Some(Strand::OT), b"TAG"), None);
    assert_eq!(indel_vote(&haplotypes, Some(Strand::OB), b"TAG"), Some((Haplotype::Two, 0)));
    assert_eq!(indel_vote(&haplotypes, Some(Strand::OT), b"CAG"), Some((Haplotype::One, 0)));
    assert_eq!(indel_vote(&haplotypes, None, b"TAG"), None);
}

fn vote(set: i32, haplotype: Haplotype, weight: u8) -> Vote {
    Vote { set, haplotype, weight }
}

#[test]
fn the_phase_set_with_the_widest_margin_tags_the_fragment() {
    use Haplotype::{One, Two};
    let tag = |set, haplotype| Some(HaplotypeTag { set, haplotype });

    assert_eq!(assign([]), None);
    assert_eq!(assign([vote(7, Two, 30)]), tag(7, Two));
    assert_eq!(assign([vote(7, One, 30), vote(7, Two, 30)]), None, "votes that cancel");
    assert_eq!(assign([vote(7, One, 30), vote(7, Two, 31)]), tag(7, Two));
    assert_eq!(assign([vote(7, One, 0)]), None, "a base of quality 0 says nothing");
    // Set 9 is more one-sided than set 7.
    assert_eq!(
        assign([vote(7, One, 40), vote(7, Two, 30), vote(9, Two, 20), vote(9, Two, 20)]),
        tag(9, Two)
    );
    // A tie between sets goes to the lower one, whatever the order.
    assert_eq!(assign([vote(9, Two, 30), vote(7, One, 30)]), tag(7, One));
    assert_eq!(assign([vote(7, One, 30), vote(9, Two, 30)]), tag(7, One));
}

/// A forward first-in-pair read (OT) of `seq` at `pos`, every base at
/// quality 30.
fn read(name: &[u8], pos: i64, cigar: Vec<Cigar>, seq: &[u8]) -> Record {
    let mut record = Record::new();
    record.set(name, Some(&CigarString(cigar.into())), seq, &vec![30; seq.len()]);
    record.set_pos(pos);
    record.set_flags(99);
    record
}

fn snv(pos: u32, set: i32, alleles: [Base; 2]) -> PhasedSite {
    PhasedSite { pos, set, alleles: SiteAlleles::Snv(alleles) }
}

#[test]
fn a_read_votes_at_the_sites_it_aligns_a_base_to() -> Result<()> {
    let sites = [
        snv(99, 100, [A, G]),  // before the read
        snv(102, 100, [A, G]), // read shows G
        snv(105, 100, [A, G]), // deleted
        snv(108, 100, [G, A]), // read shows G, after the deletion
        snv(109, 200, [A, G]), // read shows T: neither
        snv(110, 200, [A, C]), // past the read's end
    ];
    //            100  101  102  103 | 104-106 deleted | 107 108 109
    let seq = b"AAGAAGT";
    let record = read(b"r", 100, vec![Cigar::Match(4), Cigar::Del(3), Cigar::Match(3)], seq);

    let votes = votes_of(&record, &sites)?;
    assert_eq!(
        &votes[..],
        &[vote(100, Haplotype::Two, 30), vote(100, Haplotype::One, 30)],
        "one vote per site the read shows an allele of"
    );
    Ok(())
}

/// The vote weighs what the deciding base weighs, also after an insertion
/// has moved the read's bases against the reference.
#[test]
fn a_vote_carries_the_quality_of_the_base_that_decides_it() -> Result<()> {
    let mut record =
        read(b"r", 100, vec![Cigar::Match(2), Cigar::Ins(2), Cigar::Match(2)], b"AATTGA");
    let quals = [10, 11, 12, 13, 14, 15];
    let seq = record.seq().as_bytes();
    let cigar = record.cigar().take();
    record.set(b"r", Some(&cigar), &seq, &quals);

    // Reference 102 is the read's fifth base.
    let votes = votes_of(&record, &[snv(102, 100, [C, G])])?;
    assert_eq!(&votes[..], &[vote(100, Haplotype::Two, 14)]);

    // An insertion of TT after reference 101, against a reference that goes
    // on GA: the first base after the anchor decides.
    let site = PhasedSite {
        pos: 101,
        set: 100,
        alleles: SiteAlleles::Indel([b"GAC".to_vec(), b"TTGAC".to_vec()]),
    };
    let votes = votes_of(&record, &[site])?;
    assert_eq!(&votes[..], &[vote(100, Haplotype::Two, 12)]);
    Ok(())
}

/// A soft clip is still the molecule's sequence: an insertion the aligner
/// clipped at the read's end is read from it.
#[test]
fn an_indel_is_read_from_a_soft_clip() -> Result<()> {
    let record = read(b"r", 100, vec![Cigar::Match(3), Cigar::SoftClip(4)], b"ACATTGA");
    let site = PhasedSite {
        pos: 102,
        set: 100,
        alleles: SiteAlleles::Indel([b"GACCA".to_vec(), b"TTGACCA".to_vec()]),
    };
    assert_eq!(&votes_of(&record, &[site])?[..], &[vote(100, Haplotype::Two, 30)]);
    Ok(())
}

/// Mates vote together, so the one over no site is tagged as well; an
/// alignment of the same name further off than [`MATE_REACH`] is not.
#[test]
fn mates_within_reach_share_one_tag() -> Result<()> {
    let sites = [snv(102, 100, [A, G])];
    let left = read(b"pair", 100, vec![Cigar::Match(5)], b"AAGAA");
    let near = read(b"pair", 100 + MATE_REACH.cast_signed(), vec![Cigar::Match(5)], b"AAAAA");
    let far = read(b"pair", 101 + MATE_REACH.cast_signed(), vec![Cigar::Match(5)], b"AAAAA");
    let stranger = read(b"other", 100, vec![Cigar::Match(5)], b"AATAA");

    let mut fragments = Fragments::default();
    for record in [&left, &near, &far, &stranger] {
        fragments.add(record, &sites)?;
    }

    let tag = Some(HaplotypeTag { set: 100, haplotype: Haplotype::Two });
    assert_eq!(fragments.tag_of(&left), tag);
    assert_eq!(fragments.tag_of(&near), tag);
    assert_eq!(fragments.tag_of(&far), None);
    assert_eq!(fragments.tag_of(&stranger), None);
    Ok(())
}

#[test]
fn tagging_replaces_the_tags_a_read_came_with() -> Result<()> {
    let sites = [snv(102, 100, [A, G])];
    let int = |record: &Record, tag: &[u8]| match record.aux(tag) {
        Ok(Aux::U8(value)) => Some(i32::from(value)),
        Ok(Aux::I32(value)) => Some(value),
        _ => None,
    };

    let mut tagged = read(b"tagged", 100, vec![Cigar::Match(5)], b"AAGAA");
    let mut untagged = read(b"untagged", 100, vec![Cigar::Match(5)], b"AATAA");
    for record in [&mut tagged, &mut untagged] {
        record.push_aux(b"HP", Aux::U8(1))?;
        record.push_aux(b"PS", Aux::I32(5))?;
        record.push_aux(b"NM", Aux::U8(0))?;
    }
    let mut fragments = Fragments::default();
    fragments.add(&tagged, &sites)?;
    fragments.add(&untagged, &sites)?;
    fragments.tag(&mut tagged)?;
    fragments.tag(&mut untagged)?;

    assert_eq!((int(&tagged, b"HP"), int(&tagged, b"PS")), (Some(2), Some(100)));
    assert_eq!((int(&untagged, b"HP"), int(&untagged, b"PS")), (None, None), "stale tags");
    assert_eq!(int(&untagged, b"NM"), Some(0), "other tags stay");
    Ok(())
}

/// The phased sites of VCF data lines on contig `chr1`, whose reference from
/// any position on is `flank`.
fn sites_of(lines: &[&str], flank: &[u8]) -> Result<Vec<Option<PhasedSite>>> {
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("calls.vcf");
    let mut vcf = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n",
    );
    for line in lines {
        vcf.push_str(&line.replace(' ', "\t"));
        vcf.push('\n');
    }
    std::fs::write(&path, vcf)?;

    let mut reader = bcf::Reader::from_path(&path)?;
    reader
        .records()
        .map(|record| PhasedSite::from_record(&record?, |_| Ok(flank.to_vec())))
        .collect()
}

#[test]
fn a_phased_het_with_a_phase_set_is_a_site() -> Result<()> {
    let sites = sites_of(
        &[
            "chr1 10 . A G . PASS . GT:PS 0|1:10",
            "chr1 20 . A G . PASS . GT:PS 1|0:10",
            "chr1 30 . A C,G . PASS . GT:PS 2|1:10",
            "chr1 40 . A G . PASS . GT:PS 0/1:.",
            "chr1 50 . A G . PASS . GT:PS 1|1:10",
            "chr1 60 . A G . PASS . GT:PS 0|1:.",
            "chr1 70 . A G . PASS . GT 0|1",
            "chr1 80 . A * . PASS . GT:PS 0|1:10",
            "chr1 90 . A . . PASS . GT:PS 0/0:.",
        ],
        b"",
    )?;
    assert_eq!(
        sites,
        [
            Some(snv(9, 10, [A, G])),
            Some(snv(19, 10, [G, A])),
            Some(snv(29, 10, [G, C])),
            None,
            None,
            None,
            None,
            None,
            None,
        ]
    );
    Ok(())
}

/// An indel site's haplotypes are spelled past the anchor: the allele's own
/// bases, then the reference after `REF`.
#[test]
fn an_indel_site_spells_both_haplotypes_into_the_reference() -> Result<()> {
    let sites = sites_of(
        &[
            "chr1 10 . A ATT . PASS . GT:PS 0|1:10",
            "chr1 20 . ACC A . PASS . GT:PS 1|0:10",
            "chr1 30 . ACC A,AGGCC . PASS . GT:PS 1|2:10",
            "chr1 40 . AC GT . PASS . GT:PS 0|1:10",
        ],
        b"gact",
    )?;
    let indel = |pos, haplotypes: [&[u8]; 2]| {
        Some(PhasedSite {
            pos,
            set: 10,
            alleles: SiteAlleles::Indel(haplotypes.map(<[u8]>::to_vec)),
        })
    };
    assert_eq!(
        sites,
        [
            indel(9, [b"GACT", b"TTGACT"]),
            indel(19, [b"GACT", b"CCGACT"]),
            indel(29, [b"GACT", b"GGCCGACT"]),
            // No shared anchor base to compare from.
            None,
        ]
    );
    Ok(())
}
