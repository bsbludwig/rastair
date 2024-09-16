pub mod cpg_buffer;

use bio::bio_types::strand::Strand;
use bio::io::bed::Reader;
use bio::io::fasta::{Index, IndexedReader};
use fxhash::FxBuildHasher;
use std::path::PathBuf;
use std::str::FromStr;
//use log::{trace, debug, error};
use anyhow::{anyhow, Result};
use log::{debug, trace, warn};

use std::{
    collections::{BTreeMap, HashMap},
    fmt::Debug,
    io::{stdout, Read, Seek, Write},
    path::Path,
};

use crate::utils::file_helpers::open_file;
use cpg_buffer::{CpgBuffer, CpgInfo};

use super::{ReadMask, ReadMaskSetting};

const INITIAL_HASH_SIZE: usize = 1_000; // how many reads in between the average two matching pairs?
const FLUSH_BUFFER_THRESHOLD: usize = 5_000;
const CPG_SEARCH_RANGE: u64 = 2_000;
const MAX_INSERT_SIZE: u64 = 10_000;

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
enum MethylationState
{
    Methylated,
    Unmethylated,
    Unknown,
}
use MethylationState::*;

struct ReadInfo
{
    contig: Vec<u8>,
    start: u64,
    flag: u16,
    cpg_info: Vec<(usize, MethylationState)>,
}

impl ReadInfo
{
    fn new(contig: Vec<u8>, start: u64, flag: u16, cpg_info: Vec<(usize, MethylationState)>)
           -> Self
    {
        ReadInfo { contig,
                   start,
                   flag,
                   cpg_info }
    }

    fn strand(&self) -> Strand
    {
        return Self::strand_from_flag(self.flag);
    }

    fn strand_from_flag(flag: u16) -> Strand
    {
        // F1
        if flag & 64 == 64 && flag & 16 == 0 {
            Strand::Forward
        }
        // R2
        else if flag & 144 == 144 {
            Strand::Forward
        }
        // R1
        else if flag & 80 == 80 {
            Strand::Reverse
        }
        // F2
        else if flag & 128 == 128 && flag & 16 == 0 {
            Strand::Reverse
        } else {
            Strand::Unknown
        }
    }
}

pub struct PatGenerator<R: Read + Seek>
{
    // map from read name to chr/(start,end)/flag/[cpg_ind/mod_status, ...]
    read_hash: HashMap<String, ReadInfo, FxBuildHasher>,
    read_reader: Reader<R>,
    cpg_buffer: CpgBuffer<R>,
    current_chromosome: Vec<u8>,
    output_buffer: BTreeMap<usize, Vec<(String, usize)>>,
    n_ot: ReadMaskSetting,
    n_ob: ReadMaskSetting,
}

impl<R: Read + Seek> PatGenerator<R>
{
    pub fn with_read_and_fasta(bed_reader: Reader<R>,
                               fasta_reader: IndexedReader<R>,
                               n_ot: ReadMaskSetting,
                               n_ob: ReadMaskSetting,
                               step_size_option: &Option<u32>)
                               -> Result<Self>
    {
        let read_hash: HashMap<String, ReadInfo, FxBuildHasher> =
            HashMap::with_capacity_and_hasher(INITIAL_HASH_SIZE, FxBuildHasher::default());

        let cpg_buffer = if let Some(step_size) = step_size_option {
            CpgBuffer::with_reader_and_stepsize(fasta_reader, *step_size as usize)?
        } else {
            CpgBuffer::with_reader(fasta_reader)?
        };
        Ok(Self { read_hash,
                  read_reader: bed_reader,
                  cpg_buffer,
                  current_chromosome: Vec::new(),
                  output_buffer: BTreeMap::new(),
                  n_ot,
                  n_ob })
    }

    pub fn process(&mut self) -> Result<()>
    {
        // Iterate over records
        for r in self.read_reader.records() {
            if let Ok(record) = r {
                let (chr, start, end, qname) = (record.chrom(),
                                                record.start(),
                                                record.end(),
                                                record.name()
                                                      .ok_or(anyhow!("No read name in record"))?);

                let this_chromosome = Vec::from(chr.as_bytes());
                if self.current_chromosome.len() == 0 {
                    self.current_chromosome = this_chromosome.clone();
                } else if self.current_chromosome != this_chromosome {
                    // Make sure we fast-forward in the bed file to the next chromosome
                    self.cpg_buffer.progress_to_contig(&this_chromosome);
                    // flush all buffers
                    flush_singletons(&mut self.read_hash, &mut self.output_buffer);
                    self.read_hash.clear();
                    flush_write_buffer_until(&mut self.output_buffer,
                                             self.current_chromosome.as_slice(),
                                             std::usize::MAX)?; // flush the cache
                    self.current_chromosome = this_chromosome.clone();
                }
                let read_length: u32 =
                    record.aux(7)
                          .ok_or(anyhow!("No read length in record {}:{}-{}", chr, start, end))
                          .unwrap()
                          .parse()?;

                let flag: u16 = record.aux(8)
                                      .ok_or(anyhow!("No flag in record {}:{}-{}", chr, start, end))
                                      .unwrap()
                                      .parse()?;
                // create a single representation of CpGs in this read
                let mods = parse_mod_str(record.aux(11).unwrap_or(""));
                let unmods = parse_mod_str(record.aux(12).unwrap_or(""));
                let snps = parse_mod_str(record.aux(13).unwrap_or(""));
                // Convert in-read coordinates to absolute CpG indices

                let strand = ReadInfo::strand_from_flag(flag);
                let read_mask = match strand {
                    Strand::Forward => {
                        if flag & 64 == 64
                        // first in pair
                        {
                            self.n_ot.r1
                        } else {
                            self.n_ot.r2
                        }
                    }
                    Strand::Reverse => {
                        if flag & 64 == 64
                        // first in pair
                        {
                            self.n_ob.r1
                        } else {
                            self.n_ob.r2
                        }
                    }
                    Strand::Unknown => ReadMask(0, 0),
                };
                let all_mods: Vec<(usize, MethylationState)> = if start < end {
                    if let Some(cpg_slice) = self.cpg_buffer
                                                 .cpgs_in_range(chr.as_bytes().as_ref(), start, end)
                    {
                        match zip_mods(&mods,
                                       &unmods,
                                       &snps,
                                       strand,
                                       cpg_slice,
                                       read_mask,
                                       read_length as usize)
                        {
                            None => {
                                warn!("Skipping read {} due to parser error", qname);
                                continue;
                            }
                            Some(mods) => mods,
                        }
                    } else {
                        warn!("Failed to fetch sequence for record {}", qname);
                        Vec::new()
                    }
                } else {
                    warn!("Empty sequence for record {}", qname);
                    Vec::new()
                };

                let r2 = ReadInfo::new(this_chromosome, start, flag, all_mods);

                // Check if the read pair is in the cache already, otherwise just put data in cache and continue
                if let Some(r1) = self.read_hash.remove(qname) {
                    debug!("Processing read pair {} from {}:{} to {}:{}",
                           qname, chr, r1.start, chr, end);

                    // then combine the two pairs into an output string and put that into the output cache.
                    process_read_pair(&r1, &r2, &mut self.output_buffer);

                    /* The read cache has changed. Flush the output cache to the lowest start pos
                     * remaining in the read buffer.
                     * This is an expensive operation, so we will only do this once the output buffer
                     * has exceeded some size.
                     */
                    if self.output_buffer.len() > FLUSH_BUFFER_THRESHOLD
                       || self.read_hash.len() > FLUSH_BUFFER_THRESHOLD
                    {
                        // remove orphaned single reads that are unreasonably far away from a pair
                        // TODO make this configurable
                        let keys_to_remove: Vec<String> = self.read_hash
                            .iter()
                            .filter(|(_, v)| {
                                if (v.start + MAX_INSERT_SIZE) < start
                                {
                                    trace!("Read at {}:{} likely orphaned, more than {} from {}:{}", std::str::from_utf8(&v.contig).unwrap_or_default(), v.start, MAX_INSERT_SIZE, chr, start);
                                    true
                                }
                                else {
                                    false
                                }
                            })
                            .map(|(k, _)| String::from(k))
                            .collect();
                        debug!("Will treat {} reads in cache that are too far before current position as singletons", keys_to_remove.len());
                        for k in keys_to_remove {
                            // insert as singleton
                            if let Some(singleton) = self.read_hash.remove(&k) {
                                if let Some(new_flag) = get_mate_flag(singleton.flag) {
                                    let fake_pair = ReadInfo::new(singleton.contig.clone(),
                                                                  singleton.start,
                                                                  new_flag,
                                                                  Vec::new());
                                    process_read_pair(&singleton,
                                                      &fake_pair,
                                                      &mut self.output_buffer);
                                }
                            }
                        }

                        let right_margin = if let Some(leftmost_remaining_read) =
                            self.read_hash.values().min_by(|a, b| a.start.cmp(&b.start))
                        {
                            debug!("Left-most read in read buffer: {}:{}",
                                   chr, leftmost_remaining_read.start);
                            leftmost_remaining_read.start
                        } else {
                            debug!("No reads left in hash, flush to {}", r2.start);
                            r2.start
                        };
                        if let Some(cpgs) = self.cpg_buffer
                                                .cpgs_in_range(&self.current_chromosome,
                                                               right_margin
                                                               - std::cmp::min(right_margin,
                                                                               CPG_SEARCH_RANGE),
                                                               right_margin)
                        {
                            if let Some(last_cpg_before) = cpgs.last() {
                                let len_before = self.output_buffer.len();
                                flush_write_buffer_until(&mut self.output_buffer,
                                                         &self.current_chromosome,
                                                         last_cpg_before.index / 2 + 1)?;
                                debug!("Flushed {} positions from output buffer before {}:{}",
                                       len_before - self.output_buffer.len(),
                                       chr,
                                       last_cpg_before.index / 2 + 1);
                            } else {
                                warn!("Could not find max index in previous {} bases before {}",
                                      CPG_SEARCH_RANGE, right_margin);
                            }
                        } else {
                            warn!("No CpGs in the {} bases before {}:{}, can't flush buffer this time", CPG_SEARCH_RANGE, chr, right_margin);
                        }
                    }
                } else {
                    self.read_hash.insert(qname.to_string(), r2);
                }
            }
        }
        // Final flush. First dump all remaining unpaired reads as singletons
        flush_singletons(&mut self.read_hash, &mut self.output_buffer);
        flush_write_buffer_until(&mut self.output_buffer,
                                 self.current_chromosome.as_slice(),
                                 std::usize::MAX)?; // flush the cache
        Ok(())
    }
}

fn flush_singletons(read_hash: &mut HashMap<String, ReadInfo, FxBuildHasher>,
                    output_buffer: &mut BTreeMap<usize, Vec<(String, usize)>>)
                    -> ()
{
    // Final flush. First dump all remaining unpaired reads as singletons
    debug!("\n{}",
           read_hash.keys()
                    .map(|f| format!("Dumping singleton {}", f))
                    .collect::<Vec<String>>()
                    .join("\n"));
    for singleton in read_hash.values() {
        if let Some(new_flag) = get_mate_flag(singleton.flag) {
            let fake_pair = ReadInfo::new(singleton.contig.clone(),
                                          singleton.start,
                                          new_flag,
                                          Vec::new());
            process_read_pair(singleton, &fake_pair, output_buffer);
        }
    }
    ()
}

fn get_mate_flag(flag: u16) -> Option<u16>
{
    if flag & 64 == 64 {
        if flag & 16 == 16 {
            // F2
            Some((flag ^ 80) | 160)
        } else {
            // F1
            Some((flag ^ 96) | 144)
        }
    } else if flag & 128 == 128 {
        if flag & 16 == 16 {
            // R2
            Some((flag ^ 144) | 96)
        } else {
            // R1
            Some((flag ^ 160) | 80)
        }
    } else {
        None
    }
}
fn process_read_pair(r1: &ReadInfo,
                     r2: &ReadInfo,
                     output_buffer: &mut BTreeMap<usize, Vec<(String, usize)>>)
                     -> ()
{
    if r1.cpg_info.len() == 0 && r2.cpg_info.len() == 0 {
        debug!("Skipping fragment with no information at {}:{}",
               std::str::from_utf8(&r1.contig).unwrap_or_default(),
               r1.start);
        return;
    }

    if let Some((pos, meth_string)) = read_to_output_tuple(r1, r2) {
        debug!("Combined into pattern {} {}", pos, meth_string);
        add_to_output_buffer(output_buffer, (pos, meth_string));
    }
}
fn flush_write_buffer_until(output_buffer: &mut BTreeMap<usize, Vec<(String, usize)>>,
                            current_chromosome: &[u8],
                            end: usize)
                            -> Result<()>
{
    let mut lock = stdout().lock();
    let chr_string = std::str::from_utf8(current_chromosome).unwrap_or_default();

    loop {
        if let Some((index, mut value)) = output_buffer.pop_first() {
            if index >= end {
                // re-insert
                output_buffer.insert(index, value);
                break;
            }
            // pat documentation suggests they sort by pattern, which
            // probably improves gzip compressed of the output file
            value.sort_by(|a, b| a.0.cmp(&b.0));
            for (meth_string, count) in value.iter() {
                writeln!(lock,
                         "{}\t{}\t{}\t{}",
                         chr_string, index, meth_string, count)?;
            }
        } else {
            break;
        }
    }

    Ok(())
}

fn add_to_output_buffer(output_buffer: &mut BTreeMap<usize, Vec<(String, usize)>>,
                        search_pattern: (usize, String))
                        -> ()
{
    if let Some(strings) = output_buffer.get_mut(&search_pattern.0) {
        if let Some(index) = strings.iter()
                                    .position(|(pattern, _)| pattern == &search_pattern.1)
        {
            strings[index].1 += 1;
        } else {
            // new pattern at this position, add to list
            strings.push((search_pattern.1, 1));
        }
    } else {
        output_buffer.insert(search_pattern.0, vec![(search_pattern.1, 1)]);
    }
}

fn min_max_indices(read1: &ReadInfo, read2: &ReadInfo) -> (usize, usize)
{
    let default_info_min = (std::usize::MAX, Unknown);
    let default_info_max = (0, Unknown);
    (std::cmp::min(read1.cpg_info.first().unwrap_or(&default_info_min).0,
                   read2.cpg_info.first().unwrap_or(&default_info_min).0),
     std::cmp::max(read1.cpg_info.last().unwrap_or(&default_info_max).0,
                   read2.cpg_info.last().unwrap_or(&default_info_max).0))
}

fn read_to_output_tuple(read1: &ReadInfo, read2: &ReadInfo) -> Option<(usize, String)>
{
    assert_eq!(read1.strand(),
               read2.strand(),
               "Trying to merge reads from different orientations");
    assert_eq!(read1.contig, read2.contig, "Read contigs do not match");

    if read1.cpg_info.is_empty() && read2.cpg_info.is_empty() {
        return None;
    }

    let (mut min_index, max_index) = min_max_indices(read1, read2);

    let site_count = max_index - min_index + 1;

    if site_count == 0 {
        return None;
    }

    // Create an empty string with no methlation info yet
    let mut output = vec!['.' as u8; site_count];
    // the first x1 CpGs are from r1, the last x2 CpGs are from r2
    for cpg in read1.cpg_info.iter().chain(read2.cpg_info.iter()) {
        let methyl_char: u8 = match cpg.1 {
            Methylated => b'C',
            Unmethylated => b'T',
            Unknown => b'.',
        };
        let cur_char = output[cpg.0 - min_index];
        if cur_char != methyl_char && cur_char != b'.' {
            // TODO: This should be configurable in some way.
            // In cfDNA, this would lead to loss of methylation info due to the R2 being poorly converted
            warn!("Non-matching methylation state at index {}, will set to Unknown",
                  cpg.0);
            output[cpg.0 - min_index] = b'.';
        } else {
            output[cpg.0 - min_index] = match cpg.1 {
                Methylated => b'C',
                Unmethylated => b'T',
                Unknown => b'.',
            };
        }
    }

    // trim leading and trailing .'s
    loop {
        if output.len() > 0 && output[0] == b'.' {
            min_index += 1;
            output.remove(0);
        } else {
            break;
        }
    }
    if output.len() == 0 {
        return None;
    }
    loop {
        if output.len() == 0 {
            break;
        }
        let c = output.len() - 1;
        if output[c] == b'.' {
            output.remove(c);
        } else {
            break;
        }
    }
    if output.len() == 0 {
        return None;
    }

    Some((min_index, String::from_utf8(output).unwrap_or_default()))
}
fn parse_mod_str(mod_str: &str) -> Vec<u8>
{
    // expect format [1,6,10]

    // empty or malformed, return empty
    if mod_str.trim_start().len() == 0 {
        return Vec::new();
    }

    // remove brackets if they exist, split by comma, turn into u8
    mod_str.split(&['[', ']', ',', ' '])
           .filter(|f| f.len() > 0)
           .map(|f| f.parse::<u8>().unwrap_or_default()) // this is a bit risky, no error checking
           .collect()
}

fn zip_mods<'a>(mods: &Vec<u8>,
                unmod: &Vec<u8>,
                snps: &Vec<u8>,
                strand: Strand,
                cpg_info: Vec<&'a CpgInfo>,
                read_mask: ReadMask,
                read_length: usize)
                -> Option<Vec<(usize, MethylationState)>>
{
    let filtered_cpgs: Vec<&CpgInfo> = cpg_info.into_iter()
                                               .filter(|f| f.strand() == strand)
                                               .collect();

    // TODO I will have to relax this, as currently I can't deal with deletions in the read, which will lead to out-of-sync
    // errors in rare instances!
    if filtered_cpgs.len() > mods.len() + unmod.len() + snps.len() {
        warn!("Likely indel covering a CpG, CpG count in read does not match CpG count in reference");
        return None;
    } else if filtered_cpgs.len() < mods.len() + unmod.len() + snps.len() {
        warn!("Fewer CpGs in slice than in read: {} vs {}",
              filtered_cpgs.len(),
              mods.len() + unmod.len() + snps.len());
        return None;
    }

    let mut mod_pairs: Vec<(u8, MethylationState)> =
        mods.iter()
            .map(|f| {
                let pos_in_read = *f;
                assert!((pos_in_read as usize) < read_length,
                        "Position {} outside read length {}",
                        pos_in_read,
                        read_length);
                if (pos_in_read as usize) < read_mask.0
                   || (read_length - pos_in_read as usize - 1) < read_mask.1
                {
                    (pos_in_read, Unknown)
                } else {
                    (pos_in_read, Methylated)
                }
            })
            .collect();
    let mut unmod_pairs: Vec<(u8, MethylationState)> =
        unmod.iter()
             .map(|f| {
                 let pos_in_read = *f;
                 assert!((pos_in_read as usize) < read_length,
                         "Position {} outside read length {}",
                         pos_in_read,
                         read_length);
                 if (pos_in_read as usize) < read_mask.0
                    || (read_length - pos_in_read as usize - 1) < read_mask.1
                 {
                     (pos_in_read, Unknown)
                 } else {
                     (pos_in_read, Unmethylated)
                 }
             })
             .collect();
    let mut snp_pairs: Vec<(u8, MethylationState)> = snps.iter().map(|f| (*f, Unknown)).collect();
    mod_pairs.append(&mut unmod_pairs);
    mod_pairs.append(&mut snp_pairs);
    // Sort methylated sutes by position
    mod_pairs.sort_by(|a, b| a.0.cmp(&b.0));
    // remove leading/trailing Unknowns
    let mut found_value = false;
    mod_pairs = mod_pairs.into_iter()
                         .filter(|entry| {
                             if found_value {
                                 true
                             } else {
                                 if entry.1 == Unknown {
                                     false
                                 } else {
                                     found_value = true;
                                     true
                                 }
                             }
                         })
                         .collect();
    found_value = false;
    mod_pairs = mod_pairs.into_iter()
                         .rev()
                         .filter(|entry| {
                             if found_value {
                                 true
                             } else {
                                 if entry.1 == Unknown {
                                     false
                                 } else {
                                     found_value = true;
                                     true
                                 }
                             }
                         })
                         .rev()
                         .collect();
    // PAT format counts CpG positions as 1, rather than counting C and G separately, so I need to fix the coordinate system
    Some(mod_pairs.iter()
                  .enumerate()
                  .map(|(i, m)| (filtered_cpgs[i].index / 2 + 1, m.1))
                  .collect())
}

#[allow(non_snake_case)]
pub fn run_bed2pat<P: AsRef<Path> + std::fmt::Debug>(fasta_path: P,
                                                     read_bed: P,
                                                     nOT_option: &Option<String>,
                                                     nOB_option: &Option<String>,
                                                     chunk_size_option: &Option<u32>)
                                                     -> Result<()>
{
    let read_file = open_file(&read_bed)?;
    let read_reader = Reader::new(read_file);
    let fasta_file = open_file(&fasta_path)?;
    let index_path =
        PathBuf::from(format!("{}.fai", fasta_path.as_ref().to_str().unwrap_or_default()).as_str());
    let fasta_index = Index::from_file(&index_path)?;
    let indexed_reader = IndexedReader::with_index(fasta_file, fasta_index);

    #[allow(non_snake_case)]
    let mut n_ot =
        ReadMaskSetting::from_str(nOT_option.as_ref().unwrap_or(&"0,0,0,0".to_string())).unwrap();
    let mut n_ob =
        ReadMaskSetting::from_str(nOB_option.as_ref().unwrap_or(&"0,0,0,0".to_string())).unwrap();

    n_ot.r2 = ReadMask(n_ot.r2.1, n_ot.r2.0); // R2 is mapped in reverse
    n_ob.r1 = ReadMask(n_ob.r1.1, n_ob.r1.0); // F2 is mapped in reverse

    let mut generator = PatGenerator::with_read_and_fasta(read_reader,
                                                          indexed_reader,
                                                          n_ot,
                                                          n_ob,
                                                          chunk_size_option)?;
    generator.process()
}
/*====================================================
 = Unit Tests
====================================================*/
#[cfg(test)]
mod tests
{
    //use std::io::Cursor;

    use super::*;
    use anyhow::{Ok, Result};

    // Need to write some tests that test the actual reading, but it's a hassle
    const _COORD_BED: &[u8] = b"chr1\t0\t1\t0\t+
chr1\t1\t2\t1\t-
chr1\t6\t7\t2\t+
chr1\t7\t8\t3\t-
chr1\t10\t11\t4\t+
chr1\t11\t12\t5\t-
chr1\t20\t21\t6\t+
chr1\t21\t22\t7\t-
chr1\t100\t101\t8\t+
chr1\t101\t102\t9\t-
chr1\t110\t111\t10\t+
chr1\t111\t112\t11\t-
chr2\t10\t11\t0\t+
chr2\t11\t12\t1\t-
chr2\t16\t17\t2\t+
chr2\t17\t18\t3\t-
chr2\t18\t19\t4\t+
chr2\t19\t20\t5\t-
";
    #[test]
    fn can_parse_mod_str() -> Result<()>
    {
        let mut test_str = "[1,2,3]";
        let mut parsed = parse_mod_str(test_str);

        assert_eq!(parsed.len(), 3);
        assert_eq!(parsed, [1, 2, 3]);

        test_str = "[]";
        parsed = parse_mod_str(test_str);
        assert_eq!(parsed.len(), 0);

        test_str = "";
        parsed = parse_mod_str(test_str);
        assert_eq!(parsed.len(), 0);

        test_str = "[1, 2,3]";
        parsed = parse_mod_str(test_str);
        assert_eq!(parsed.len(), 3);
        assert_eq!(parsed, [1, 2, 3]);

        Ok(())
    }

    #[test]
    fn can_zip_mods() -> Result<()>
    {
        let cpg_info = vec![CpgInfo::new(Vec::from("chr1".as_bytes()), 1, 0),
                            CpgInfo::new(Vec::from("chr1".as_bytes()), 2, 1),
                            CpgInfo::new(Vec::from("chr1".as_bytes()), 4, 2),
                            CpgInfo::new(Vec::from("chr1".as_bytes()), 5, 3),
                            CpgInfo::new(Vec::from("chr1".as_bytes()), 10, 4),
                            CpgInfo::new(Vec::from("chr1".as_bytes()), 11, 5),
                            CpgInfo::new(Vec::from("chr1".as_bytes()), 20, 6),
                            CpgInfo::new(Vec::from("chr1".as_bytes()), 21, 7)];
        let mods: Vec<u8> = [1, 4].to_vec();
        let unmods: Vec<u8> = [2, 3].to_vec();
        let snps: Vec<u8> = Vec::new();
        let read_mask = ReadMask(0, 0);
        assert_eq!(zip_mods(&mods,
                            &unmods,
                            &snps,
                            Strand::Forward,
                            cpg_info.iter().collect(),
                            read_mask,
                            6).expect("Failed to zip"),
                   [(1, Methylated),
                    (2, Unmethylated),
                    (3, Unmethylated),
                    (4, Methylated)]);

        let mods: Vec<u8> = [1].to_vec();
        let unmods: Vec<u8> = [2, 3].to_vec();
        let snps: Vec<u8> = [4].to_vec();
        assert_eq!(zip_mods(&mods,
                            &unmods,
                            &snps,
                            Strand::Forward,
                            cpg_info.iter().collect(),
                            read_mask,
                            6).expect("Failed to zip"),
                   [(1, Methylated),
                    (2, Unmethylated),
                    (3, Unmethylated),
                    (4, Unknown)]);
        Ok(())
    }

    #[test]
    fn test_read_strand() -> Result<()>
    {
        let mut read = ReadInfo::new(Vec::from("chr1".as_bytes()), 0, 99, Vec::new());
        assert_eq!(read.strand(), Strand::Forward);
        read.flag = 147;
        assert_eq!(read.strand(), Strand::Forward);
        read.flag = 83;
        assert_eq!(read.strand(), Strand::Reverse);
        read.flag = 163;
        assert_eq!(read.strand(), Strand::Reverse);
        Ok(())
    }
}
