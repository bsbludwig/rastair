//! Generated scenarios for whole-run properties: a patch of tiled read pairs
//! with planted SNVs, CpGs, de-novo CpGs and indels, drawn by hegel so a
//! failure shrinks to the fewest, simplest events. The contig runs on past
//! the patch, so a property can put more reads there.

use super::scenario::{Hap, READ_LEN, Scenario};
use hegel::{TestCase, generators as gs};

pub const ORIGIN: usize = 600;
pub const PATCH: usize = 1_400;
pub const CONTIG_LEN: usize = ORIGIN + PATCH + 1_800;
/// Events sit where the tiling gives full depth (see [`Scenario::tiled`]).
pub const FIRST_EVENT: usize = ORIGIN + 350;
pub const LAST_EVENT: usize = ORIGIN + PATCH - 350;
/// Room for the widest event (a three-base deletion moved back by up to three
/// bases) before the next one.
const MIN_GAP: usize = 8;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Carrier {
    One,
    Two,
    Both,
}

impl Carrier {
    pub fn haps(self) -> &'static [Hap] {
        match self {
            Self::One => &[Hap::One],
            Self::Two => &[Hap::Two],
            Self::Both => &Hap::BOTH,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event {
    /// The reference base's `alt`-th other base, which may be a transition
    /// TAPS confounds with methylation.
    Snv {
        alt: usize,
        carrier: Carrier,
    },
    /// A reference `CG`.
    Cpg,
    /// A reference `TG` whose `T>C` makes a CpG.
    DenovoC {
        carrier: Carrier,
    },
    /// A reference `CA` whose `A>G` makes a CpG.
    DenovoG {
        carrier: Carrier,
    },
    Insertion {
        bases: Vec<u8>,
        carrier: Carrier,
    },
    Deletion {
        len: usize,
        carrier: Carrier,
    },
}

#[derive(Debug, Clone)]
pub struct Plan {
    pub seed: u64,
    /// Ascending positions, at least [`MIN_GAP`] apart.
    pub events: Vec<(usize, Event)>,
    /// Divides 50, so both reads of every pair start on the same grid.
    pub step: usize,
    pub copies: usize,
    /// Per haplotype, the chance a read pair is methylated.
    pub methylation: [f64; 2],
    pub indel_flanks: IndelFlanks,
}

/// Whether an indel may join a `C` and a `G` into a CpG.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndelFlanks {
    /// The reference and the inserted bases as drawn.
    Random,
    /// An `A` before the indel and a `T` after it, and no `G` inserted.
    Neutral,
}

/// A planted SNV, 0-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PlantedSnv {
    pub pos: usize,
    pub reference: u8,
    pub alt: u8,
    pub carrier: Carrier,
}

/// A planted indel after the base at `anchor`, 0-based: `inserted` bases, or
/// `deleted` reference bases.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlantedIndel {
    pub anchor: usize,
    pub inserted: Vec<u8>,
    pub deleted: usize,
    pub carrier: Carrier,
}

/// A planted CpG by its `C`, 0-based: in the reference (`carrier` `None`) or
/// made by an SNV on the carrier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PlantedCpg {
    pub c: usize,
    pub carrier: Option<Carrier>,
}

impl Plan {
    pub fn scenario(&self) -> Scenario {
        let mut scenario = self.reference();
        for PlantedSnv { pos, alt, carrier, .. } in self.snvs() {
            scenario = match carrier {
                Carrier::One => scenario.het(pos, alt, Hap::One),
                Carrier::Two => scenario.het(pos, alt, Hap::Two),
                Carrier::Both => scenario.hom(pos, alt),
            };
        }
        for PlantedIndel { anchor, inserted, deleted, carrier } in self.indels() {
            for &hap in carrier.haps() {
                scenario = if deleted > 0 {
                    scenario.deletion(anchor, deleted, hap)
                } else {
                    scenario.insertion(anchor, &inserted, hap)
                };
            }
        }
        let scenario = scenario.tiled(ORIGIN, ORIGIN + PATCH, self.step, self.copies);
        Hap::BOTH.iter().zip(self.methylation).fold(scenario, |scenario, (&hap, level)| {
            scenario.methylation_level(hap, level, self.seed ^ hap_seed(hap))
        })
    }

    /// The reference with the planted contexts written in, and no reads.
    pub fn reference(&self) -> Scenario {
        let scenario = self.events.iter().fold(
            Scenario::with_seed(CONTIG_LEN, self.seed),
            |scenario, (pos, event)| match event {
                Event::Cpg => scenario.reference_at(*pos, b"CG"),
                Event::DenovoC { .. } => scenario.reference_at(*pos, b"TG"),
                Event::DenovoG { .. } => scenario.reference_at(*pos, b"CA"),
                _ => scenario,
            },
        );
        match self.indel_flanks {
            IndelFlanks::Random => scenario,
            IndelFlanks::Neutral => self.indels().iter().fold(scenario, |scenario, indel| {
                scenario
                    .reference_at(indel.anchor, b"A")
                    .reference_at(indel.anchor + 1 + indel.deleted, b"T")
            }),
        }
    }

    pub fn snvs(&self) -> Vec<PlantedSnv> {
        let reference = self.reference();
        let base = |pos: usize| reference.reference().get(pos).copied().unwrap_or(b'N');
        self.events
            .iter()
            .filter_map(|(pos, event)| {
                let (pos, alt, carrier) = match event {
                    Event::Snv { alt, carrier } => {
                        let alt = b"ACGT".iter().copied().filter(|&b| b != base(*pos)).nth(*alt)?;
                        (*pos, alt, *carrier)
                    }
                    Event::DenovoC { carrier } => (*pos, b'C', *carrier),
                    Event::DenovoG { carrier } => (pos + 1, b'G', *carrier),
                    _ => return None,
                };
                Some(PlantedSnv { pos, reference: base(pos), alt, carrier })
            })
            .collect()
    }

    pub fn indels(&self) -> Vec<PlantedIndel> {
        self.events
            .iter()
            .filter_map(|(pos, event)| match event {
                Event::Insertion { bases, carrier } => Some(PlantedIndel {
                    anchor: *pos,
                    inserted: match self.indel_flanks {
                        IndelFlanks::Random => bases.clone(),
                        IndelFlanks::Neutral => bases
                            .iter()
                            .map(|&base| if base == b'G' { b'A' } else { base })
                            .collect(),
                    },
                    deleted: 0,
                    carrier: *carrier,
                }),
                Event::Deletion { len, carrier } => Some(PlantedIndel {
                    anchor: self.deletion_anchor(*pos, *len),
                    inserted: Vec::new(),
                    deleted: *len,
                    carrier: *carrier,
                }),
                _ => None,
            })
            .collect()
    }

    pub fn cpgs(&self) -> Vec<PlantedCpg> {
        self.events
            .iter()
            .filter_map(|(pos, event)| match event {
                Event::Cpg => Some(PlantedCpg { c: *pos, carrier: None }),
                Event::DenovoC { carrier } | Event::DenovoG { carrier } => {
                    Some(PlantedCpg { c: *pos, carrier: Some(*carrier) })
                }
                _ => None,
            })
            .collect()
    }

    /// The anchor of a deletion drawn at `pos`, moved onto the read start it
    /// would delete: no read can start inside a deletion.
    fn deletion_anchor(&self, pos: usize, len: usize) -> usize {
        (pos + 1..=pos + len)
            .find(|deleted| (deleted - ORIGIN).is_multiple_of(self.step))
            .unwrap_or(pos)
    }

    /// The pairs of each haplotype and strand over an event.
    pub fn depth_per_hap_and_strand(&self) -> usize {
        2 * READ_LEN / self.step * self.copies
    }
}

fn hap_seed(hap: Hap) -> u64 {
    match hap {
        Hap::One => 0x9e37_79b9_7f4a_7c15,
        Hap::Two => 0xbf58_476d_1ce4_e5b9,
    }
}

/// What a plan may plant.
#[derive(Debug, Clone, Copy)]
pub struct Kinds {
    pub indels: bool,
    pub max_events: usize,
    pub indel_flanks: IndelFlanks,
}

#[hegel::composite]
pub fn carriers(tc: &TestCase) -> Carrier {
    tc.draw(gs::sampled_from([Carrier::One, Carrier::Two, Carrier::Both].as_slice()))
}

#[hegel::composite]
pub fn events(tc: &TestCase, indels: bool) -> Event {
    let kinds: u8 = if indels { 6 } else { 4 };
    match tc.draw(gs::integers::<u8>().max_value(kinds - 1)) {
        0 => Event::Snv {
            alt: tc.draw(gs::integers::<usize>().max_value(2)),
            carrier: tc.draw(carriers()),
        },
        1 => Event::Cpg,
        2 => Event::DenovoC { carrier: tc.draw(carriers()) },
        3 => Event::DenovoG { carrier: tc.draw(carriers()) },
        4 => {
            let len = tc.draw(gs::integers::<usize>().min_value(1).max_value(3));
            let bases = (0..len).map(|_| tc.draw(gs::sampled_from(b"ACGT".as_slice()))).collect();
            Event::Insertion { bases, carrier: tc.draw(carriers()) }
        }
        _ => Event::Deletion {
            len: tc.draw(gs::integers::<usize>().min_value(1).max_value(3)),
            carrier: tc.draw(carriers()),
        },
    }
}

#[hegel::composite]
pub fn plans(tc: &TestCase, kinds: Kinds) -> Plan {
    let count = tc.draw(gs::integers::<usize>().max_value(kinds.max_events));
    let mut planted = Vec::with_capacity(count);
    let mut pos = FIRST_EVENT;
    for _ in 0..count {
        pos += tc.draw(gs::integers::<usize>().max_value(120));
        if pos > LAST_EVENT {
            break;
        }
        planted.push((pos, tc.draw(events(kinds.indels))));
        pos += MIN_GAP;
    }
    let level = || gs::sampled_from([0.0, 0.25, 0.5, 1.0].as_slice());
    Plan {
        seed: tc.draw(gs::integers::<u64>().min_value(1)),
        events: planted,
        step: tc.draw(gs::sampled_from([25, 50].as_slice())),
        copies: tc.draw(gs::integers::<usize>().min_value(1).max_value(2)),
        methylation: [tc.draw(level()), tc.draw(level())],
        indel_flanks: kinds.indel_flanks,
    }
}
