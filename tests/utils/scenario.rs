//! Two haplotypes, the TAPS reads a library of them would hold, and the phase
//! a correct caller has to report.
//!
//! A test states the biology — which haplotype carries which allele, which one
//! is methylated, where its read pairs sit — and the scenario writes the
//! reference and the BAM. Because the haplotypes are known, so is the answer:
//! [`Scenario::check_phasing`] reads the VCF back and fails on any record
//! phased against the truth.

use color_eyre::eyre::{Context as _, OptionExt as _, Result, bail, ensure, eyre};
use rust_htslib::bam::{self, Record, header::HeaderRecord, record::Cigar, record::CigarString};
use std::{collections::BTreeMap, path::Path};

pub const CONTIG: &str = "scenario";
pub const READ_LEN: usize = 100;
/// Where a template's second read starts, relative to its first, unless the
/// test says otherwise. Far enough apart that a pair is two separate reads.
pub const INSERT: usize = 250;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hap {
    One,
    Two,
}

impl Hap {
    const BOTH: [Self; 2] = [Self::One, Self::Two];

    const fn index(self) -> usize {
        match self {
            Self::One => 0,
            Self::Two => 1,
        }
    }
}

/// The TAPS strand a template comes from. A methylated cytosine reads as `T`
/// on the top strand, and its CpG partner `G` as `A` on the bottom strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Top,
    Bottom,
}

impl Strand {
    const BOTH: [Self; 2] = [Self::Top, Self::Bottom];
}

/// One read pair: the haplotype and strand it comes from, and where its two
/// reads start.
#[derive(Debug, Clone, Copy)]
struct Template {
    hap: Hap,
    strand: Strand,
    left: usize,
    right: usize,
}

#[derive(Debug, Clone)]
pub struct Scenario {
    reference: Vec<u8>,
    /// Per position, the base each haplotype carries where it differs from the
    /// reference.
    variants: BTreeMap<usize, [u8; 2]>,
    methylated: [bool; 2],
    templates: Vec<Template>,
}

impl Scenario {
    /// A pseudo-random reference of `len` bases with no CpG in it, so the only
    /// CpGs are the ones a test writes with [`Self::reference_at`] or creates
    /// with a variant.
    pub fn new(len: usize) -> Self {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut reference: Vec<u8> = Vec::with_capacity(len);
        while reference.len() < len {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let base = match state >> 62 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            if !(base == b'G' && reference.last() == Some(&b'C')) {
                reference.push(base);
            }
        }
        Self { reference, variants: BTreeMap::new(), methylated: [false; 2], templates: Vec::new() }
    }

    /// Overwrite the reference from `pos` on, e.g. with `b"CG"` for a CpG.
    pub fn reference_at(mut self, pos: usize, bases: &[u8]) -> Self {
        for (offset, &base) in bases.iter().enumerate() {
            if let Some(slot) = self.reference.get_mut(pos + offset) {
                *slot = base;
            }
        }
        self
    }

    /// A heterozygous SNV: `hap` carries `alt`, the other haplotype the
    /// reference.
    pub fn het(mut self, pos: usize, alt: u8, hap: Hap) -> Self {
        let reference = self.reference_base(pos);
        let mut bases = [reference; 2];
        if let Some(slot) = bases.get_mut(hap.index()) {
            *slot = alt;
        }
        self.variants.insert(pos, bases);
        self
    }

    /// A heterozygous transversion on `hap` that creates no CpG, for a site
    /// whose alleles should not matter to the test.
    pub fn snv(self, pos: usize, hap: Hap) -> Self {
        let before = pos.checked_sub(1).map(|p| self.reference_base(p));
        let after = self.reference_base(pos + 1);
        let candidates: &[u8] = match self.reference_base(pos) {
            b'A' | b'G' => b"CT",
            _ => b"AG",
        };
        let alt = candidates
            .iter()
            .copied()
            .find(|&alt| !(alt == b'C' && after == b'G') && !(alt == b'G' && before == Some(b'C')))
            .unwrap_or(b'T');
        self.het(pos, alt, hap)
    }

    /// Every CpG on `hap`, including any its variants create, is methylated.
    pub fn methylated(mut self, hap: Hap) -> Self {
        if let Some(slot) = self.methylated.get_mut(hap.index()) {
            *slot = true;
        }
        self
    }

    /// `copies` read pairs from `hap` on `strand`, the first read at `left` and
    /// its mate [`INSERT`] bases on.
    pub fn pairs(self, hap: Hap, strand: Strand, left: usize, copies: usize) -> Self {
        self.pairs_to(hap, strand, left, left + INSERT, copies)
    }

    pub fn pairs_to(
        mut self,
        hap: Hap,
        strand: Strand,
        left: usize,
        right: usize,
        copies: usize,
    ) -> Self {
        let template = Template { hap, strand, left, right };
        self.templates.extend(std::iter::repeat_n(template, copies));
        self
    }

    /// `copies` read pairs from each haplotype on each strand.
    pub fn balanced(self, left: usize, right: usize, copies: usize) -> Self {
        let mut scenario = self;
        for hap in Hap::BOTH {
            for strand in Strand::BOTH {
                scenario = scenario.pairs_to(hap, strand, left, right, copies);
            }
        }
        scenario
    }

    pub fn reference(&self) -> &[u8] {
        &self.reference
    }

    /// The `call --region` argument covering the whole contig.
    pub fn region(&self) -> String {
        format!("--region={CONTIG}:1-{}", self.reference.len())
    }

    fn reference_base(&self, pos: usize) -> u8 {
        self.reference.get(pos).copied().unwrap_or(b'N')
    }

    fn base(&self, hap: Hap, pos: usize) -> u8 {
        self.variants
            .get(&pos)
            .and_then(|bases| bases.get(hap.index()).copied())
            .unwrap_or_else(|| self.reference_base(pos))
    }

    /// One read as a BAM stores it: in reference orientation, with every
    /// methylated CpG converted on the read's own strand.
    fn read(&self, hap: Hap, strand: Strand, start: usize) -> Vec<u8> {
        let methylated = self.methylated.get(hap.index()).copied().unwrap_or(false);
        (start..start + READ_LEN)
            .map(|pos| {
                let here = self.base(hap, pos);
                let cpg_c = here == b'C' && self.base(hap, pos + 1) == b'G';
                let cpg_g =
                    here == b'G' && pos.checked_sub(1).map(|p| self.base(hap, p)) == Some(b'C');
                match strand {
                    Strand::Top if methylated && cpg_c => b'T',
                    Strand::Bottom if methylated && cpg_g => b'A',
                    _ => here,
                }
            })
            .collect()
    }

    /// Write `ref.fa` and its index into `dir`, returning the FASTA path.
    pub fn write_fasta(&self, dir: &Path) -> Result<std::path::PathBuf> {
        let fasta = dir.join("ref.fa");
        let header = format!(">{CONTIG}\n");
        let mut body = header.clone().into_bytes();
        body.extend_from_slice(&self.reference);
        body.push(b'\n');
        std::fs::write(&fasta, body).wrap_err("write fasta")?;
        let len = self.reference.len();
        let fai = format!("{CONTIG}\t{len}\t{}\t{len}\t{}\n", header.len(), len + 1);
        std::fs::write(dir.join("ref.fa.fai"), fai).wrap_err("write fai")?;
        Ok(fasta)
    }

    /// Write the reads as a sorted, indexed BAM.
    pub fn write_bam(&self, path: &Path) -> Result<()> {
        for template in &self.templates {
            ensure!(
                template.right + READ_LEN <= self.reference.len(),
                "a read at {} runs off the {}-base reference",
                template.right,
                self.reference.len()
            );
        }

        let mut header = bam::Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", CONTIG);
        sq.push_tag(b"LN", self.reference.len());
        header.push_record(&sq);

        let cigar = CigarString(vec![Cigar::Match(u32::try_from(READ_LEN)?)].into());
        let mut records = Vec::new();
        for (index, template) in self.templates.iter().enumerate() {
            let name = format!("pair_{index}");
            let span = i64::try_from(template.right + READ_LEN - template.left)?;
            // Top strand: first read forward on the left (99), mate reverse
            // (147). Bottom strand: the left read is the forward second read
            // (163), its mate the reverse first read (83).
            let (left_flags, right_flags) = match template.strand {
                Strand::Top => (99, 147),
                Strand::Bottom => (163, 83),
            };
            for (start, mate, flags, insert_size) in [
                (template.left, template.right, left_flags, span),
                (template.right, template.left, right_flags, -span),
            ] {
                let seq = self.read(template.hap, template.strand, start);
                let mut record = Record::new();
                record.set(name.as_bytes(), Some(&cigar), &seq, &[40u8; READ_LEN]);
                record.set_tid(0);
                record.set_pos(i64::try_from(start)?);
                record.set_mapq(60);
                record.set_flags(flags);
                record.set_mtid(0);
                record.set_mpos(i64::try_from(mate)?);
                record.set_insert_size(insert_size);
                records.push(record);
            }
        }
        records.sort_by_key(Record::pos);

        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)?;
        for record in &records {
            writer.write(record)?;
        }
        drop(writer);
        bam::index::build(path, None, bam::index::Type::Bai, 1)?;
        Ok(())
    }

    /// Check every record of a VCF against the haplotypes, and report how the
    /// heterozygous sites were phased.
    ///
    /// Fails when a phased record is not a heterozygote of the scenario, when a
    /// block puts one haplotype's alleles on both sides (a switch error), when a
    /// `PS` does not name the first record of its block, or when that record
    /// does not open with the lower allele (`0|1`, not `1|0`).
    pub fn check_phasing(&self, vcf: &str) -> Result<Phasing> {
        let mut blocks: BTreeMap<usize, Vec<(usize, Hap)>> = BTreeMap::new();
        let mut opens_low: BTreeMap<usize, bool> = BTreeMap::new();
        let mut unphased = Vec::new();

        for line in vcf.lines().filter(|line| !line.starts_with('#')) {
            let record = VcfLine::parse(line)?;
            let truth = self.variants.get(&record.pos);
            let Some(gt) = record.gt.split_once('|') else {
                if truth.is_some_and(|[one, two]| one != two) {
                    unphased.push(record.pos);
                }
                continue;
            };
            let [one, two] = truth.copied().ok_or_else(|| {
                eyre!("{} is phased but carries no variant in the scenario: {line}", record.pos)
            })?;
            ensure!(one != two, "{} is phased but homozygous in the scenario: {line}", record.pos);

            let allele = |index: &str| -> Result<u8> {
                let index: usize = index.parse().wrap_err_with(|| format!("GT in {line}"))?;
                record.alleles.get(index).copied().ok_or_eyre("GT names a missing allele")
            };
            opens_low.insert(record.pos, gt.0 < gt.1);
            let (first, second) = (allele(gt.0)?, allele(gt.1)?);
            // Which haplotype the record's haplotype 1 is.
            let hap = if (first, second) == (one, two) {
                Hap::One
            } else if (first, second) == (two, one) {
                Hap::Two
            } else {
                bail!(
                    "{} is phased as {first}|{second}, the haplotypes carry {one}/{two}",
                    record.pos
                )
            };
            let set = record.ps.ok_or_else(|| eyre!("phased without PS: {line}"))?;
            blocks.entry(set).or_default().push((record.pos, hap));
        }

        let mut phased = BTreeMap::new();
        for (set, sites) in blocks {
            let (first_pos, first_hap) = *sites.first().ok_or_eyre("a block is never empty")?;
            ensure!(
                first_pos + 1 == set,
                "PS={set} does not name its block's first site {first_pos}"
            );
            ensure!(
                opens_low.get(&first_pos) == Some(&true),
                "block PS={set} does not open with its lower allele first"
            );
            if let Some((pos, _)) = sites.iter().find(|(_, hap)| *hap != first_hap) {
                bail!("switch error at {pos} in block PS={set}: {sites:?}");
            }
            phased.insert(first_pos, sites.into_iter().map(|(pos, _)| pos).collect());
        }
        Ok(Phasing { blocks: phased, unphased })
    }
}

/// How a VCF phased a scenario's heterozygous sites, all consistent with the
/// truth by the time [`Scenario::check_phasing`] returns one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Phasing {
    /// Each block's sites (0-based), keyed on its first.
    pub blocks: BTreeMap<usize, Vec<usize>>,
    /// Heterozygous sites with a record but no phase.
    pub unphased: Vec<usize>,
}

/// The few columns of a VCF data line phasing is judged on.
struct VcfLine {
    /// 0-based.
    pos: usize,
    /// `REF` then each `ALT`, single bases only.
    alleles: Vec<u8>,
    gt: String,
    /// 1-based, as written.
    ps: Option<usize>,
}

impl VcfLine {
    fn parse(line: &str) -> Result<Self> {
        let columns: Vec<&str> = line.split('\t').collect();
        let column = |index: usize| columns.get(index).copied().ok_or_eyre("short VCF line");
        let pos = column(1)?.parse::<usize>()?.checked_sub(1).ok_or_eyre("POS 0")?;
        let mut alleles: Vec<u8> = column(3)?.bytes().take(1).collect();
        alleles.extend(column(4)?.split(',').filter_map(|alt| alt.bytes().next()));
        let format_value = |key: &str| -> Option<&str> {
            let index = column(8).ok()?.split(':').position(|candidate| candidate == key)?;
            column(9).ok()?.split(':').nth(index)
        };
        Ok(Self {
            pos,
            alleles,
            gt: format_value("GT").unwrap_or(".").to_owned(),
            ps: format_value("PS").and_then(|ps| ps.parse().ok()),
        })
    }
}
