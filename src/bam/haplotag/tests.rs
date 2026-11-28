use super::*;
use Base::{A, C, G, T};
use color_eyre::eyre::ensure;
use hegel::{TestCase, generators as gs};
use seqair::{bam::aux::find_tag, reader::AlignmentStream};

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

/// The records of SAM lines on contig `chr1`, their fields separated by
/// spaces; read through the stream reader, since a record has no public
/// constructor.
fn records(lines: &[String]) -> Result<Vec<RawBamRecord>> {
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n");
    for line in lines {
        sam.push_str(&line.replace(' ', "\t"));
        sam.push('\n');
    }
    let mut stream = AlignmentStream::from_reader(std::io::Cursor::new(sam.into_bytes()))?;
    let mut records = Vec::new();
    let mut record = RawBamRecord::default();
    while stream.read_record(&mut record)? {
        records.push(record.clone());
    }
    ensure!(records.len() == lines.len(), "a SAM line did not parse");
    Ok(records)
}

/// The SAM line of a forward first-in-pair read (OT) of `seq` at the 0-based
/// `pos`, every base at quality 30.
fn sam(name: &str, pos: u64, cigar: &str, seq: &str) -> String {
    format!("{name} 99 chr1 {} 60 {cigar} = {} 0 {seq} {}", pos + 1, pos + 1, "?".repeat(seq.len()))
}

fn read(name: &str, pos: u64, cigar: &str, seq: &str) -> Result<RawBamRecord> {
    records(&[sam(name, pos, cigar, seq)])?.pop().wrap_err("no record")
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
    // 100 101 102 103 | 104-106 deleted | 107 108 109
    let record = read("r", 100, "4M3D3M", "AAGAAGT")?;

    let votes = votes_of(&record, &sites);
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
    // Qualities 10 to 15.
    let line = "r 99 chr1 101 60 2M2I2M = 101 0 AATTGA +,-./0".to_owned();
    let record = records(&[line])?.pop().wrap_err("no record")?;

    // Reference 102 is the read's fifth base.
    let votes = votes_of(&record, &[snv(102, 100, [C, G])]);
    assert_eq!(&votes[..], &[vote(100, Haplotype::Two, 14)]);

    // An insertion of TT after reference 101, against a reference that goes
    // on GA: the first base after the anchor decides.
    let site = PhasedSite {
        pos: 101,
        set: 100,
        alleles: SiteAlleles::Indel([b"GAC".to_vec(), b"TTGAC".to_vec()]),
    };
    let votes = votes_of(&record, &[site]);
    assert_eq!(&votes[..], &[vote(100, Haplotype::Two, 12)]);
    Ok(())
}

/// A read stored without qualities still votes, every base alike.
#[test]
fn a_read_without_qualities_votes() -> Result<()> {
    let line = "r 99 chr1 101 60 5M = 101 0 AAGAA *".to_owned();
    let record = records(&[line])?.pop().wrap_err("no record")?;
    let votes = votes_of(&record, &[snv(102, 100, [A, G])]);
    assert_eq!(&votes[..], &[vote(100, Haplotype::Two, NO_QUALITY_WEIGHT)]);
    Ok(())
}

/// A soft clip is still the molecule's sequence: an insertion the aligner
/// clipped at the read's end is read from it.
#[test]
fn an_indel_is_read_from_a_soft_clip() -> Result<()> {
    let record = read("r", 100, "3M4S", "ACATTGA")?;
    let site = PhasedSite {
        pos: 102,
        set: 100,
        alleles: SiteAlleles::Indel([b"GACCA".to_vec(), b"TTGACCA".to_vec()]),
    };
    assert_eq!(&votes_of(&record, &[site])[..], &[vote(100, Haplotype::Two, 30)]);
    Ok(())
}

fn int_tag(record: &RawBamRecord, tag: &[u8; 2]) -> Option<i64> {
    find_tag(record.aux().as_bytes(), *tag)?.as_i64()
}

/// A read's name and start, and its `(HP, PS)`.
type Tagged = (String, u64, Option<(i64, i64)>);

fn tags_of(batches: &[(Vec<RawBamRecord>, u8)]) -> Vec<Tagged> {
    batches
        .iter()
        .flat_map(|(records, _)| records)
        .map(|record| {
            (
                String::from_utf8_lossy(record.qname()).into_owned(),
                record.pos().map_or(0, |pos| pos.as_u64()),
                int_tag(record, b"HP").zip(int_tag(record, b"PS")),
            )
        })
        .collect()
}

/// Push `lines` as one batch, numbered `with`, voting at `sites`.
fn push(
    tagger: &mut Tagger<u8>,
    sites: &[PhasedSite],
    lines: &[String],
    with: u8,
) -> Result<Vec<(Vec<RawBamRecord>, u8)>> {
    let records = records(lines)?;
    let mut votes = BatchVotes::default();
    for record in &records {
        votes.add(record, sites);
    }
    tagger.push(records, votes, with)
}

/// Mates vote together, so the one over no site is tagged as well; an
/// alignment of the same name further off than [`MATE_REACH`] is not.
#[test]
fn mates_within_reach_share_one_tag() -> Result<()> {
    let sites = [snv(102, 100, [A, G])];
    let mut tagger = Tagger::default();
    let held = push(
        &mut tagger,
        &sites,
        &[
            sam("pair", 100, "5M", "AAGAA"),
            sam("other", 100, "5M", "AATAA"),
            sam("pair", 100 + MATE_REACH, "5M", "AAAAA"),
            sam("pair", 101 + MATE_REACH, "5M", "AAAAA"),
        ],
        1,
    )?;
    assert!(held.is_empty(), "a batch waits for the input to move past it");

    let tag = Some((2, 100));
    assert_eq!(
        tags_of(&tagger.finish()?),
        [
            ("pair".to_owned(), 100, tag),
            ("other".to_owned(), 100, None),
            ("pair".to_owned(), 100 + MATE_REACH, tag),
            ("pair".to_owned(), 101 + MATE_REACH, None),
        ]
    );
    Ok(())
}

/// The mate arrives a batch later, and the read over no site is the earlier
/// one: its batch is held until the input has moved out of reach, then both
/// come out tagged alike, in input order, each batch with what rode along.
#[test]
fn a_batch_waits_for_the_mates_in_the_next_one() -> Result<()> {
    let sites = [snv(402, 100, [A, G])];
    let mut tagger = Tagger::default();

    let first = push(&mut tagger, &sites, &[sam("pair", 100, "5M", "AAAAA")], 1)?;
    assert!(first.is_empty());
    let second = push(&mut tagger, &sites, &[sam("pair", 400, "5M", "AAGAA")], 2)?;
    assert!(second.is_empty(), "the first batch's last read still reaches the second's first");

    // Out of reach of the second batch's last read, so both are settled.
    let far = 401 + MATE_REACH;
    let third = push(&mut tagger, &sites, &[sam("far", far, "5M", "AAAAA")], 3)?;
    assert_eq!(third.iter().map(|(_, with)| *with).collect::<Vec<_>>(), [1, 2]);
    let tag = Some((2, 100));
    assert_eq!(tags_of(&third), [("pair".to_owned(), 100, tag), ("pair".to_owned(), 400, tag)]);
    assert!(tagger.votes.is_empty(), "votes out of reach are dropped: {:?}", tagger.votes);

    assert_eq!(tags_of(&tagger.finish()?), [("far".to_owned(), far, None)]);
    Ok(())
}

/// Votes are kept past their batch for the mates that come later: the batch
/// with the voting read is handed on while its mate's batch, which reaches
/// further, is still held.
#[test]
fn a_handed_on_batch_still_votes_for_later_mates() -> Result<()> {
    let sites = [snv(102, 100, [A, G])];
    let tag = Some((2, 100));
    let mut tagger = Tagger::default();
    push(&mut tagger, &sites, &[sam("pair", 100, "5M", "AAGAA")], 1)?;
    let mate = [sam("pair", 900, "5M", "AAAAA"), sam("far", 5000, "5M", "AAAAA")];
    assert!(push(&mut tagger, &sites, &mate, 2)?.is_empty());

    let handed_on = push(&mut tagger, &sites, &[sam("next", 5001, "5M", "AAAAA")], 3)?;
    assert_eq!(tags_of(&handed_on), [("pair".to_owned(), 100, tag)]);

    let later = push(&mut tagger, &sites, &[sam("last", 9000, "5M", "AAAAA")], 4)?;
    assert_eq!(
        tags_of(&later),
        [
            ("pair".to_owned(), 900, tag),
            ("far".to_owned(), 5000, None),
            ("next".to_owned(), 5001, None)
        ]
    );
    assert!(tagger.votes.is_empty(), "votes out of reach are dropped: {:?}", tagger.votes);
    Ok(())
}

/// Unplaced records end a sorted input, so the batches before them are
/// handed on, and they come out untagged themselves.
#[test]
fn unplaced_records_settle_the_batches_before_them() -> Result<()> {
    let sites = [snv(102, 100, [A, G])];
    let mut tagger = Tagger::default();
    push(&mut tagger, &sites, &[sam("read", 100, "5M", "AAGAA")], 1)?;
    let unplaced = "lost 4 * 0 0 * * 0 0 ACGT ???? HP:i:1 PS:i:7".to_owned();
    let handed_on = push(&mut tagger, &sites, &[unplaced], 2)?;
    assert_eq!(tags_of(&handed_on), [("read".to_owned(), 100, Some((2, 100)))]);
    assert_eq!(tags_of(&tagger.finish()?), [("lost".to_owned(), 0, None)]);
    Ok(())
}

#[test]
fn tagging_replaces_the_tags_a_read_came_with() -> Result<()> {
    let sites = [snv(102, 100, [A, G])];
    let stale = " HP:i:1 PS:i:5 NM:i:0";
    let mut tagger = Tagger::default();
    push(
        &mut tagger,
        &sites,
        &[
            sam("tagged", 100, "5M", "AAGAA") + stale,
            sam("untagged", 100, "5M", "AATAA") + stale,
            sam("untouched", 100, "5M", "AATAA") + " NM:i:0",
        ],
        1,
    )?;
    let batches = tagger.finish()?;
    let [tagged, untagged, untouched] = batches
        .first()
        .and_then(|(records, _)| <&[RawBamRecord; 3]>::try_from(&records[..]).ok())
        .wrap_err("three records")?;

    assert_eq!((int_tag(tagged, b"HP"), int_tag(tagged, b"PS")), (Some(2), Some(100)));
    assert_eq!((int_tag(untagged, b"HP"), int_tag(untagged, b"PS")), (None, None), "stale tags");
    for record in [tagged, untagged, untouched] {
        assert_eq!(int_tag(record, b"NM"), Some(0), "other tags stay");
    }
    Ok(())
}
