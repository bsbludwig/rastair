use bitflags::bitflags;
use clap::value_parser;
use clio::ClioPath;
use color_eyre::eyre::{Result, eyre};
use rust_htslib::bam::{self, Read as _, record::Cigar};
use smallvec::SmallVec;
use std::path::Path;
use tracing::{debug, info, instrument, trace, warn};
use tracing_subscriber::layer::SubscriberExt as _;

use crate::utils::{Base, TryAsBase as _, file_helpers::open_maybe_bgzip};

#[instrument(skip(path))]
pub fn read(path: &Path) -> Result<()> {
    let fetch_range = 6_105_700..6_105_800;

    let mut fasta = {
        let fasta_file = open_maybe_bgzip("test_data/test.fasta.gz")?;
        let fasta_index = bio::io::fasta::Index::from_file(&path.with_extension("fasta.fai"))
            .map_err(|err| eyre!(Box::new(err)))?;
        bio::io::fasta::IndexedReader::with_index(fasta_file, fasta_index)
    };
    // indexed_reader.fetch("chr19", fetch_range.start, fetch_range.end + 1)?;

    let mut bam = bam::IndexedReader::from_path(path)?;
    bam.set_threads(8)?;
    bam.fetch(("chr19", fetch_range.start, fetch_range.end + 1))?;

    let mut seq = Vec::new();

    for pile in bam
        .pileup()
        .filter_map(|p| p.ok())
        .filter(|p| fetch_range.contains(&(p.pos() as u64)))
        .take(100)
    {
        fasta.fetch("chr19", pile.pos() as u64, pile.pos() as u64 + 2)?;
        fasta.read(&mut seq)?;
        let bases = SeenBases(
            pile.alignments()
                .filter_map(|a| {
                    let pos = a.qpos()?;
                    let record = a.record();
                    if !record.is_proper_pair() {
                        // fixme: maybe be more lenient here
                        return None;
                    }
                    if record.is_quality_check_failed() {
                        return None;
                    }
                    // fixme: understand this better:
                    // if record.cigar().iter().any(|c| matches!(c, Cigar::SoftClip(_))) {
                    //     return None;
                    // }

                    Some(SeenBase {
                        base: record.seq()[pos].as_base().unwrap(),
                        qual: record.qual()[pos],
                        mapq: record.mapq(),
                        reverse: record.is_reverse(),
                        at_fringe: pos == 0 || pos == record.seq().len() - 1,
                    })
                })
                .collect(),
        );
        let reference_base = seq[0].as_base()?;
        let next_base = seq.get(1).and_then(|x| x.as_base().ok());
        if bases.interesting() {
            // info!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "found pile of interest");
            let pileup = InterestingPileup { pos: pile.pos(), bases, reference_base, next_base };
            info!(?pileup, "interesting pileup");
        } else if bases.matches(reference_base) {
            // Matches reference base
            // boring.
            // trace!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "pile matches reference");
        } else {
            warn!(
                ?bases,
                pos = pile.pos(),
                ?reference_base,
                ?next_base,
                "pile does not match reference but is also not interesting"
            );
        }
    }

    return Ok(());
}

#[derive(Debug)]
struct InterestingPileup {
    pos: u32,
    bases: SeenBases,
    reference_base: Base,
    next_base: Option<Base>,
}

struct SeenBases(SmallVec<SeenBase, 20>);

impl std::fmt::Debug for SeenBases {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.0.iter()).finish()
    }
}

struct SeenBase {
    base: Base,
    qual: u8,
    mapq: u8,
    reverse: bool,
    at_fringe: bool,
}

impl std::fmt::Debug for SeenBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if f.alternate() {
            // use the alternate format with `{:#?}` for no color
            write!(f, "{}", self.base)?;
        } else {
            write!(f, "{}", self.base.display_colored())?;
        }
        write!(
            f,
            " {}{} {}/{}",
            if self.reverse { "rev" } else { "fwd" },
            if self.at_fringe { " fr" } else { "" },
            self.qual,
            self.mapq,
        )
    }
}

impl SeenBases {
    fn matches(&self, base: Base) -> bool {
        self.0.iter().all(|b| b.base == base)
    }

    fn interesting(&self) -> bool {
        let mut counter = Counter::default();
        for b in &self.0 {
            match b.base {
                Base::A => counter.a += 1,
                Base::C => counter.c += 1,
                Base::G => counter.g += 1,
                Base::T => counter.t += 1,
            }
        }
        counter.interesting()
    }
}

#[derive(Debug, Default)]
struct Counter {
    c: usize,
    t: usize,
    a: usize,
    g: usize,
}

impl Counter {
    /// Interesting if there are multiple different bases seen
    fn interesting(&self) -> bool {
        let mut count = 0;
        if self.c > 0 {
            count += 1;
        }
        if self.t > 0 {
            count += 1;
        }
        if self.a > 0 {
            count += 1;
        }
        if self.g > 0 {
            count += 1;
        }
        count > 1
    }
}

impl FromIterator<u8> for Counter {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut counter = Counter { c: 0, t: 0, a: 0, g: 0 };
        for c in iter {
            match c {
                b'C' => counter.c += 1,
                b'T' => counter.t += 1,
                b'A' => counter.a += 1,
                b'G' => counter.g += 1,
                _ => {}
            }
        }
        counter
    }
}
