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
    pub const BOTH: [Self; 2] = [Self::One, Self::Two];

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
    /// Whether this molecule's CpGs are methylated, where the test set it
    /// apart from its haplotype's.
    methylated: Option<bool>,
}

/// An indel after the base at its position.
#[derive(Debug, Clone)]
enum Indel {
    Insertion(Vec<u8>),
    /// This many reference bases.
    Deletion(usize),
}

#[derive(Debug, Clone)]
pub struct Scenario {
    reference: Vec<u8>,
    /// Per position, the base each haplotype carries where it differs from the
    /// reference.
    variants: BTreeMap<usize, [u8; 2]>,
    /// Per anchor position, the indel each haplotype carries there.
    indels: BTreeMap<usize, [Option<Indel>; 2]>,
    methylated: [bool; 2],
    templates: Vec<Template>,
    /// Read names are this followed by the template's index.
    read_name_prefix: String,
}

impl Scenario {
    /// A pseudo-random reference of `len` bases with no CpG in it, so the only
    /// CpGs are the ones a test writes with [`Self::reference_at`] or creates
    /// with a variant.
    pub fn new(len: usize) -> Self {
        Self::with_seed(len, 0x2545_f491_4f6c_dd1d)
    }

    /// [`Self::new`] with another reference, one per non-zero `seed`.
    pub fn with_seed(len: usize, seed: u64) -> Self {
        let mut state = seed.max(1);
        let mut reference: Vec<u8> = Vec::with_capacity(len);
        while reference.len() < len {
            state = xorshift(state);
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
        Self {
            reference,
            variants: BTreeMap::new(),
            indels: BTreeMap::new(),
            methylated: [false; 2],
            templates: Vec::new(),
            read_name_prefix: "pair_".to_owned(),
        }
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

    /// A homozygous SNV: both haplotypes carry `alt`.
    pub fn hom(mut self, pos: usize, alt: u8) -> Self {
        self.variants.insert(pos, [alt; 2]);
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

    /// An insertion of `bases` after `pos`, on `hap`; call it for each
    /// haplotype that carries one.
    pub fn insertion(self, pos: usize, bases: &[u8], hap: Hap) -> Self {
        self.indel(pos, Indel::Insertion(bases.to_vec()), hap)
    }

    /// A deletion of the `len` bases after `pos`, on `hap`; call it for each
    /// haplotype that carries one.
    pub fn deletion(self, pos: usize, len: usize, hap: Hap) -> Self {
        self.indel(pos, Indel::Deletion(len), hap)
    }

    fn indel(mut self, pos: usize, indel: Indel, hap: Hap) -> Self {
        if let Some(slot) = self.indels.entry(pos).or_default().get_mut(hap.index()) {
            *slot = Some(indel);
        }
        self
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
        let template = Template { hap, strand, left, right, methylated: None };
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

    /// `copies` read pairs from each haplotype on each strand every `step`
    /// bases, the first starting at `from`, the last ending by `to`. From
    /// `from + INSERT + READ_LEN` to `to - INSERT - READ_LEN` a position lies
    /// under `2 * READ_LEN / step * copies` reads of each haplotype and strand
    /// when `step` divides [`READ_LEN`].
    pub fn tiled(self, from: usize, to: usize, step: usize, copies: usize) -> Self {
        let mut scenario = self;
        let mut left = from;
        while left + INSERT + READ_LEN <= to {
            scenario = scenario.balanced(left, left + INSERT, copies);
            left += step.max(1);
        }
        scenario
    }

    /// Each read pair of `hap` added so far is methylated, at every CpG it
    /// covers, with probability `level`, drawn from `seed`.
    pub fn methylation_level(mut self, hap: Hap, level: f64, seed: u64) -> Self {
        let mut state = seed.max(1);
        for template in self.templates.iter_mut().filter(|template| template.hap == hap) {
            state = xorshift(state);
            template.methylated = Some(unit_interval(state) < level);
        }
        self
    }

    /// The same read pairs in another order, drawn from `seed`: the BAM holds
    /// the reads of one position in another order, under other names.
    pub fn shuffled(mut self, seed: u64) -> Self {
        let mut state = seed.max(1);
        for index in (1..self.templates.len()).rev() {
            state = xorshift(state);
            let other = usize::try_from(state % (index as u64 + 1)).unwrap_or(index);
            self.templates.swap(index, other);
        }
        self
    }

    /// Read names start with `prefix` instead of `pair_`.
    pub fn read_names(mut self, prefix: &str) -> Self {
        prefix.clone_into(&mut self.read_name_prefix);
        self
    }

    /// The scenario seen from the other strand: the reverse complement of the
    /// reference and of both haplotypes, every read pair on the opposite TAPS
    /// strand at the mirrored place. Position `p` becomes `len - 1 - p`.
    pub fn mirrored(&self) -> Result<Self> {
        ensure!(self.indels.is_empty(), "a mirrored indel would need a new anchor");
        let len = self.reference.len();
        let mirror = |pos: usize| len - 1 - pos;
        let templates = self
            .templates
            .iter()
            .map(|template| -> Result<Template> {
                let mirror_start = |start: usize| {
                    len.checked_sub(READ_LEN + start).ok_or_eyre("a read runs off the reference")
                };
                Ok(Template {
                    strand: match template.strand {
                        Strand::Top => Strand::Bottom,
                        Strand::Bottom => Strand::Top,
                    },
                    left: mirror_start(template.right)?,
                    right: mirror_start(template.left)?,
                    ..*template
                })
            })
            .collect::<Result<_>>()?;
        Ok(Self {
            reference: self.reference.iter().rev().map(|&base| complement(base)).collect(),
            variants: self
                .variants
                .iter()
                .map(|(&pos, bases)| (mirror(pos), bases.map(complement)))
                .collect(),
            indels: BTreeMap::new(),
            methylated: self.methylated,
            templates,
            read_name_prefix: self.read_name_prefix.clone(),
        })
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

    /// `hap` as reads see it: each base with the reference position it aligns
    /// to, `None` inside an insertion.
    fn haplotype(&self, hap: Hap) -> Vec<(Option<usize>, u8)> {
        let mut bases = Vec::with_capacity(self.reference.len());
        let mut pos = 0;
        while pos < self.reference.len() {
            bases.push((Some(pos), self.base(hap, pos)));
            match self.indel_at(pos, hap) {
                Some(Indel::Insertion(inserted)) => {
                    bases.extend(inserted.iter().map(|&base| (None, base)));
                }
                Some(Indel::Deletion(len)) => pos += len,
                None => {}
            }
            pos += 1;
        }
        bases
    }

    /// One read as a BAM stores it: in reference orientation, with every
    /// methylated CpG converted on the read's own strand, and the CIGAR that
    /// aligns it.
    fn read(&self, template: &Template, start: usize) -> Result<(Vec<u8>, CigarString)> {
        let Template { hap, strand, .. } = *template;
        let methylated = template
            .methylated
            .unwrap_or_else(|| self.methylated.get(hap.index()).copied().unwrap_or(false));
        let haplotype = self.haplotype(hap);
        let first = haplotype
            .iter()
            .position(|&(pos, _)| pos == Some(start))
            .ok_or_else(|| eyre!("a read cannot start at {start}, which {hap:?} deletes"))?;
        let base_at = |index: Option<usize>| index.and_then(|i| haplotype.get(i)).map(|&(_, b)| b);

        let mut seq = Vec::with_capacity(READ_LEN);
        let mut cigar: Vec<Cigar> = Vec::new();
        let mut last_pos: Option<usize> = None;
        for (index, &(pos, here)) in haplotype.iter().enumerate().skip(first).take(READ_LEN) {
            let cpg_c = here == b'C' && base_at(index.checked_add(1)) == Some(b'G');
            let cpg_g = here == b'G' && base_at(index.checked_sub(1)) == Some(b'C');
            seq.push(match strand {
                Strand::Top if methylated && cpg_c => b'T',
                Strand::Bottom if methylated && cpg_g => b'A',
                _ => here,
            });

            let op = match (pos, last_pos) {
                (None, _) => Cigar::Ins(1),
                (Some(pos), Some(last)) if pos > last + 1 => {
                    cigar.push(Cigar::Del(u32::try_from(pos - last - 1)?));
                    Cigar::Match(1)
                }
                (Some(_), _) => Cigar::Match(1),
            };
            last_pos = pos.or(last_pos);
            match (cigar.last_mut(), op) {
                (Some(Cigar::Match(n)), Cigar::Match(1)) | (Some(Cigar::Ins(n)), Cigar::Ins(1)) => {
                    *n += 1;
                }
                _ => cigar.push(op),
            }
        }
        Ok((seq, CigarString(cigar.into())))
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

        let mut records = Vec::new();
        for (index, template) in self.templates.iter().enumerate() {
            let name = format!("{}{index}", self.read_name_prefix);
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
                let (seq, cigar) = self.read(template, start)?;
                let mut record = Record::new();
                record.set(name.as_bytes(), Some(&cigar), &seq, &vec![40u8; seq.len()]);
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
    /// `PS` does not name a record of its block, or when the block's first
    /// record does not open with the lower allele (`0|1`, not `1|0`).
    pub fn check_phasing(&self, vcf: &str) -> Result<Phasing> {
        // Each record's position, haplotype, and whether it opens low.
        let mut blocks: BTreeMap<usize, Vec<(usize, Hap, bool)>> = BTreeMap::new();
        let mut unphased = Vec::new();

        for line in vcf.lines().filter(|line| !line.starts_with('#')) {
            let record = VcfLine::parse(line)?;
            let truth = self.truth(&record);
            let Some(gt) = record.gt.split_once('|') else {
                if truth.as_ref().is_some_and(|[one, two]| one != two) {
                    unphased.push(record.pos);
                }
                continue;
            };
            let [one, two] = truth.ok_or_else(|| {
                eyre!("{} is phased but carries no variant in the scenario: {line}", record.pos)
            })?;
            ensure!(one != two, "{} is phased but homozygous in the scenario: {line}", record.pos);

            let allele = |index: &str| -> Result<&str> {
                let index: usize = index.parse().wrap_err_with(|| format!("GT in {line}"))?;
                record
                    .alleles
                    .get(index)
                    .map(String::as_str)
                    .ok_or_eyre("GT names a missing allele")
            };
            let (first, second) = (allele(gt.0)?, allele(gt.1)?);
            // Which haplotype the record's haplotype 1 is.
            let hap = if (first, second) == (one.as_str(), two.as_str()) {
                Hap::One
            } else if (first, second) == (two.as_str(), one.as_str()) {
                Hap::Two
            } else {
                bail!(
                    "{} is phased as {first}|{second}, the haplotypes carry {one}/{two}",
                    record.pos
                )
            };
            let set = record.ps.ok_or_else(|| eyre!("phased without PS: {line}"))?;
            blocks.entry(set).or_default().push((record.pos, hap, gt.0 < gt.1));
        }

        let mut phased = BTreeMap::new();
        for (set, sites) in blocks {
            let (first_pos, first_hap, opens_low) =
                *sites.first().ok_or_eyre("a block is never empty")?;
            ensure!(
                sites.iter().any(|(pos, ..)| pos + 1 == set),
                "PS={set} does not name a site of its block"
            );
            ensure!(opens_low, "block PS={set} does not open with its lower allele first");
            if let Some((pos, ..)) = sites.iter().find(|(_, hap, _)| *hap != first_hap) {
                bail!("switch error at {pos} in block PS={set}: {sites:?}");
            }
            phased.insert(first_pos, sites.into_iter().map(|(pos, ..)| pos).collect());
        }
        Ok(Phasing { blocks: phased, unphased })
    }

    /// What each haplotype carries at a record, spelled the way the record
    /// spells its alleles: an indel record is anchored on the base before.
    fn truth(&self, record: &VcfLine) -> Option<[String; 2]> {
        if !record.is_indel() {
            return Some(self.variants.get(&record.pos)?.map(|base| char::from(base).to_string()));
        }
        let pos = record.pos;
        self.indels.get(&pos)?;
        // REF spans the longest deletion at the anchor, and every allele is
        // spelled against it.
        let deleted = Hap::BOTH
            .iter()
            .filter_map(|&hap| match self.indel_at(pos, hap)? {
                Indel::Deletion(len) => Some(*len),
                Indel::Insertion(_) => None,
            })
            .max()
            .unwrap_or(0);
        let spell = |hap: Hap| -> Option<String> {
            let (inserted, kept): (&[u8], _) = match self.indel_at(pos, hap) {
                Some(Indel::Insertion(bases)) => (bases, pos + 1),
                Some(Indel::Deletion(len)) => (&[], pos + 1 + len),
                None => (&[], pos + 1),
            };
            let tail = self.reference.get(kept..=pos + deleted)?;
            Some(
                std::iter::once(self.reference_base(pos))
                    .chain(inserted.iter().copied())
                    .chain(tail.iter().copied())
                    .map(char::from)
                    .collect(),
            )
        };
        Some([spell(Hap::One)?, spell(Hap::Two)?])
    }

    fn indel_at(&self, pos: usize, hap: Hap) -> Option<&Indel> {
        self.indels.get(&pos)?.get(hap.index())?.as_ref()
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
    /// `REF` then each `ALT`.
    alleles: Vec<String>,
    gt: String,
    /// 1-based, as written.
    ps: Option<usize>,
}

impl VcfLine {
    fn is_indel(&self) -> bool {
        self.alleles.iter().any(|allele| allele.len() != 1)
    }

    fn parse(line: &str) -> Result<Self> {
        let columns: Vec<&str> = line.split('\t').collect();
        let column = |index: usize| columns.get(index).copied().ok_or_eyre("short VCF line");
        let pos = column(1)?.parse::<usize>()?.checked_sub(1).ok_or_eyre("POS 0")?;
        let mut alleles = vec![column(3)?.to_owned()];
        alleles.extend(column(4)?.split(',').map(str::to_owned));
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

fn xorshift(mut state: u64) -> u64 {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    state
}

/// The top 53 bits of `state` as a number in `[0, 1)`.
fn unit_interval(state: u64) -> f64 {
    (state >> 11) as f64 / (1_u64 << 53) as f64
}

pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}
