pub mod cpg_buffer;

use bgzip::{index::BGZFIndex, read::IndexedBGZFReader, BGZFReader};
use bio::io::fasta::IndexedReader;
use bio::io::{bed::Reader, fasta::Index};
use bio::bio_types::strand::Strand;
use std::path::PathBuf;
use std::str::FromStr;
use fxhash::FxBuildHasher;
//use log::{trace, debug, error};
use anyhow::{anyhow, bail, Result};
use log::{debug, warn};

use std::{collections::{HashMap, VecDeque}, fmt::Debug, fs::File, io::{stdout, Read, Seek, Write}, path::Path};

use crate::operations::bed2pat::cpg_buffer::{CpgInfo, CpgBuffer};

const INITIAL_HASH_SIZE: usize = 1000; // how many reads in between the average two matching pairs?

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
enum MethylationState {
    Methylated,
    Unmethylated,
    Unknown
}
use MethylationState::*;

use super::{ReadMask, ReadMaskSetting};

struct ReadInfo
{
    contig  : Vec<u8>,
    start   : u64,
    flag    : u16,
    cpg_info: Vec<(usize, MethylationState)>
}

impl ReadInfo
{
    fn new(contig: Vec<u8>, start: u64, flag: u16, cpg_info: Vec<(usize, MethylationState)>) -> Self
    {
        ReadInfo
        {
            contig,
            start,
            flag,
            cpg_info
        }
    }

    fn strand(&self) -> Strand
    {
        return Self::strand_from_flag(self.flag);
    }

    fn strand_from_flag(flag: u16) -> Strand
    {
        // F1
        if flag & 64 == 64 && flag & 16 == 0
        {
            Strand::Forward
        }
        // R2
        else if flag & 144 == 144
        {
            Strand::Forward
        }
        // R1
        else if flag & 80 == 80
        {
            Strand::Reverse
        }
        // F2
        else if flag & 128 == 128 && flag & 16 == 0
        {
            Strand::Reverse
        }
        else
        {
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
    output_buffer: VecDeque<(usize, Vec<(String, usize)>)>,
    n_ot: ReadMaskSetting,
    n_ob: ReadMaskSetting
}

impl <R: Read+Seek> PatGenerator<R> {
    pub fn with_read_and_fasta(bam_reader: Reader<R>, fasta_reader: IndexedReader<R>, n_ot: ReadMaskSetting, n_ob: ReadMaskSetting, step_size_option: &Option<usize>) -> Result<Self>
    {
        let read_hash: HashMap<String, ReadInfo, FxBuildHasher> = HashMap::with_capacity_and_hasher(INITIAL_HASH_SIZE, FxBuildHasher::default());

        let cpg_buffer =
            if let Some(step_size) = step_size_option
            {
                CpgBuffer::with_reader_and_stepsize(fasta_reader, *step_size)?
            }
            else
            {
                CpgBuffer::with_reader(fasta_reader)?
            };
        Ok(
            Self
            {
                read_hash,
                read_reader: bam_reader,
                cpg_buffer,
                current_chromosome: Vec::new(),
                output_buffer: VecDeque::new(),
                n_ot,
                n_ob
            }
        )
    }

    pub fn process(&mut self) -> Result<()>
    {
        /* Assume both bed files are coordinate sorted in the same way.
        * In that case, we need to do the following:
        * 1. Create a hash that contains the position, orientation and
        *    CpG identities of each fragment
        * 2. Go through the per-read bed file, and get the next read
        * 3. Move forward in the CpG bed file until I reach the end
        *    position of the current read, and save all CpG info in a
        *    map (coord->index)
        * 4. Drop coordinates before the current read start
        * 5. Check if the current read is a pair of a previous read. If
        *    so, combine them into a pair and store as a pattern.
        * 6. If the pattern is different from the previous pattern, print the
        *    previous pattern and its count. If it's the same, increment count
        */
        // This maps from read_name -> [chr, first_cpg_index, cpg_string]

        // Iterate over records
        for r in self.read_reader.records()
        {
            if let Ok(record) = r
            {
                let (chr, start, end, qname) = (record.chrom(), record.start(), record.end(), record.name().ok_or(anyhow!("No read name in record"))?);

                let this_chromosome = Vec::from(chr.as_bytes());
                if self.current_chromosome.len() == 0
                {
                    self.current_chromosome = this_chromosome.clone();
                }
                else if self.current_chromosome != this_chromosome
                {
                    // Make sure we fast-forward in the bed file to the next chromosome
                    self.cpg_buffer.progress_to_contig(&this_chromosome);
                    // flush all buffers
                    self.read_hash.clear();
                    flush_write_buffer_until(&mut self.output_buffer, self.current_chromosome.as_slice(), std::usize::MAX)?; // flush the cache
                    self.current_chromosome = this_chromosome.clone();
                }

                let flag: u16 = record
                                    .aux(7)
                                    .ok_or(anyhow!("No flag in record {}:{}-{}", chr, start, end))
                                    .unwrap()
                                    .parse()?;

                // create a single representation of CpGs in this read
                let mods = parse_mod_str(record.aux(10).unwrap_or(""));
                let unmods = parse_mod_str(record.aux(11).unwrap_or(""));
                let snps = parse_mod_str(record.aux(12).unwrap_or(""));
                // Convert in-read coordinates to absolute CpG indices

                let strand = ReadInfo::strand_from_flag(flag);
                let read_mask = match strand
                {
                    Strand::Forward => {
                        if flag & 64 == 64 // first in pair
                        {
                            self.n_ot.r1
                        }
                        else
                        {
                            self.n_ot.r2
                        }
                    },
                    Strand::Reverse => {
                        if flag & 64 == 64 // first in pair
                        {
                            self.n_ob.r1
                        }
                        else
                        {
                            self.n_ob.r2
                        }
                    },
                    Strand::Unknown => {
                        ReadMask(0, 0)
                    },
                };
                let all_mods: Vec<(usize, MethylationState)> =
                    if let Some(cpg_slice) = self.cpg_buffer.cpgs_in_range(chr.as_bytes().as_ref(), start, end)
                    {
                        zip_mods(&mods, &unmods, &snps, strand, cpg_slice, read_mask, (end-start) as usize)
                    }
                    else
                    {
                        Vec::new()
                    };

                let r2 = ReadInfo::new(this_chromosome, start, flag, all_mods);

                // Check if the read pair is in the cache already, otherwise just put data in cache and continue
                if let Some(r1) = self.read_hash.remove(qname)
                {
                    debug!("Processing read pair {} from {}:{} to {}:{}", qname, chr, r1.start, chr, end);

                    // then combine the two pairs into an output string and put that into the output cache.
                    if let Some((pos, meth_string)) = read_to_output_tuple(&r1, &r2)
                    {
                        debug!("Combined into pattern {} {}", pos, meth_string);
                        add_to_output_buffer(&mut self.output_buffer, (pos, meth_string));
                    }

                    // TODO I can't seem to get this to work, let's just cache everything for now...
                    // // The read cache has changed. Flush the output cache to the lowest CpG index in the
                    // // remaining read cache
                    // if let Some(min_index_in_cache) = min_index_in_read_buffer(&self.read_hash)
                    // {
                    //     debug!("Min index in cache: {}", min_index_in_cache);
                    //     if last_index > min_index_in_cache
                    //     {
                    //         error!("Found read in cache that's before the earliest read supposedly already processed - how can this be?");
                    //     }
                    //     last_index = min_index_in_cache;
                    //     flush_write_buffer_until(&mut self.output_buffer, &self.current_chromosome, min_index_in_cache)?;
                    // }
                    // else
                    // {
                    //     // Read cache empty. Flush all output to the last CpG in the current read
                    //     let (_, max_index) = min_max_indices(&r1, &r2);
                    //     flush_write_buffer_until(&mut self.output_buffer, &self.current_chromosome, max_index-1)?;
                    // }
                    // drop all CpGs before the beginning of the leftmost read still in the cache

                    /* flush CpG buffer before leftmost read still in buffer. Inefficient to do a full search
                    *  every time but it's not straightforward to keep track of the number of reads starting
                    *  at specific positions within the buffer. I'd have to create a dedicated data structure
                    *  for the read_buffer that has a separate table of start positions with counts, and do
                    *  "double bookkeeping"
                    */
                }
                else
                {
                    self.read_hash.insert(qname.to_string(), r2);
                }
            }
        }
        // Final flush
        flush_write_buffer_until(&mut self.output_buffer, self.current_chromosome.as_slice(), std::usize::MAX)?; // flush the cache
        Ok(())
    }
}

/// Go through all reads in the buffer and find the smallest CpG index.
/// Not efficient, but since the read buffer shouldn't be huge, this is
/// probably acceptable.
fn _min_index_in_read_buffer(read_buffer: &HashMap<String, ReadInfo, FxBuildHasher>) -> Option<usize>
{
    if read_buffer.len() > 0
    {
        let mut min_index = std::usize::MAX;
        for read in read_buffer.values()
        {
            if let Some((min_in_read,_)) = read.cpg_info.first()
            {
                if *min_in_read < min_index
                {
                    min_index = *min_in_read;
                }
            }
        }
        Some(min_index)
    }
    else
    {
        None
    }
}

fn flush_write_buffer_until(output_buffer: &mut VecDeque<(usize, Vec<(String, usize)>)>, current_chromosome: &[u8], end: usize) -> Result<()>
{
    let mut lock = stdout().lock();
    let mut write_buffer: Vec<(usize, Vec<(String, usize)>)> = Vec::new();
    output_buffer.make_contiguous();
    loop
    {
        if let Some((index, meth_strings)) = output_buffer.pop_front()
        {
            if index >= end
            {
                // put it back, end loop
                output_buffer.push_front((index, meth_strings));
                break;
            }

            write_buffer.push((index, meth_strings));
        }
        else
        {
            // Empty, stop
            break;
        }
    }
    if write_buffer.len() > 0
    {
        write_buffer.sort_by(|a, b| a.0.cmp(&b.0));

        for (index, meth_strings) in write_buffer.iter()
        {
            for (meth_string, count) in meth_strings
            {
                writeln!(lock, "{}\t{}\t{}\t{}", std::str::from_utf8(current_chromosome).unwrap_or_default(), index, meth_string, count)?;
            }
        }

    }
    Ok(())
}

fn add_to_output_buffer(output_buffer: &mut VecDeque<(usize, Vec<(String, usize)>)>, search_pattern: (usize, String)) -> ()
{
    if let Some(index_ref) = output_buffer
    .iter_mut()
    .find(|f| f.0 == search_pattern.0 )
    {
        // found an entry with this first CpG index
        // look for the same meth_string
        if let Some((index_of_pattern, _)) = index_ref.1
        .iter()
        .enumerate()
        .find(|(_i, f)| f.0 == search_pattern.1)
        {
            index_ref.1[index_of_pattern].1 += 1;
        }
        else
        {
            index_ref.1.push((search_pattern.1, 1));
        }
    }
    else
    {
        output_buffer.push_back((search_pattern.0, vec![(search_pattern.1, 1)]));
    }

    ()
}

fn min_max_indices (read1: &ReadInfo, read2: &ReadInfo) -> (usize, usize)
{
    let default_info_min = (std::usize::MAX, Unknown);
    let default_info_max = (0, Unknown);
    (std::cmp::min(read1.cpg_info.first().unwrap_or(&default_info_min).0, read2.cpg_info.first().unwrap_or(&default_info_min).0),
     std::cmp::max(read1.cpg_info.last().unwrap_or(&default_info_max).0, read2.cpg_info.last().unwrap_or(&default_info_max).0))
}

fn read_to_output_tuple(read1: &ReadInfo, read2: &ReadInfo) -> Option<(usize, String)>
{
    assert_eq!(read1.strand(), read2.strand(), "Trying to merge reads from different orientations");
    assert_eq!(read1.contig, read2.contig, "Read contigs do not match");

    if read1.cpg_info.is_empty() && read2.cpg_info.is_empty()
    {
        return None;
    }

    let (min_index, max_index) = min_max_indices(read1, read2);

    let site_count = max_index-min_index+1;
    // Create an empty string with no methlation info yet
    let mut output = vec!['.' as u8; site_count];
    // the first x1 CpGs are from r1, the last x2 CpGs are from r2
    for cpg in read1.cpg_info.iter().chain(read2.cpg_info.iter())
    {
        let methyl_char = match cpg.1 {
            Methylated   => 'C' as u8,
            Unmethylated => 'T' as u8,
            Unknown      => '.' as u8
        };
        let cur_char = output[cpg.0-min_index];
        if  cur_char != methyl_char && cur_char != '.' as u8
        {
            warn!("Non-matching methylation state at index {}", cpg.0);
        }
        output[cpg.0-min_index] = match cpg.1 {
            Methylated   => 'C' as u8,
            Unmethylated => 'T' as u8,
            Unknown      => '.' as u8
        }
    }

    Some((min_index, String::from_utf8(output).unwrap_or_default()))
}
fn parse_mod_str(mod_str: &str) -> Vec<u8>
{
    // expect format [1,6,10]

    // empty or malformed, return empty
    if mod_str.trim_start().len() == 0
    {
        return Vec::new();
    }

    // remove brackets if they exist, split by comma, turn into u8
    mod_str
        .split(&['[',']',',',' '])
        .filter(|f| f.len()>0)
        .map(|f| f.parse::<u8>().unwrap_or_default()) // this is a bit risky, no error checking
        .collect()
}

fn zip_mods<'a>(mods: &Vec<u8>, unmod: &Vec<u8>, snps: &Vec<u8>, strand: Strand, cpg_info: Vec<&'a CpgInfo>, read_mask: ReadMask, read_length: usize) -> Vec<(usize, MethylationState)>
{
    let filtered_cpgs: Vec<&CpgInfo> = cpg_info.into_iter().filter(|f| f.strand() == strand).collect();

    // TODO I will have to relax this, as currently I can't deal with deletions in the read, which will lead to out-of-sync
    // errors in rare instances!
    if filtered_cpgs.len() > mods.len()+unmod.len()+snps.len()
    {
        warn!("Likely indel covering a CpG, CpG count in read does not match CpG count in reference");
    }
    else if filtered_cpgs.len() < mods.len()+unmod.len()+snps.len()
    {
        assert!(false, "Fewer CpGs in slice than in reads");
    }

    let mut mod_pairs : Vec<(u8, MethylationState)> = mods.iter()
        .map(|f|
            {
                let pos_in_read = *f;
                if pos_in_read as usize <= read_mask.0 || pos_in_read as usize > read_length - read_mask.1
                {
                    (pos_in_read, Unknown)
                }
                else
                {
                    (pos_in_read, Methylated)
                }
            })
        .collect();
    let mut unmod_pairs : Vec<(u8, MethylationState)> = unmod.iter()
        .map(|f|
            {
                let pos_in_read = *f;
                if pos_in_read as usize <= read_mask.0 || pos_in_read as usize > read_length - read_mask.1
                {
                    (pos_in_read, Unknown)
                }
                else
                {
                    (pos_in_read, Unmethylated)
                }
            })
        .collect();
    let mut snp_pairs : Vec<(u8, MethylationState)> = snps.iter().map(|f| (*f, Unknown)).collect();
    mod_pairs.append( &mut unmod_pairs );
    mod_pairs.append( &mut snp_pairs );
    mod_pairs.sort_by(|a, b| a.0.cmp(&b.0));

    mod_pairs.iter().enumerate().map(|(i, m)| (filtered_cpgs[i].index/2+1, m.1)).collect()
}

trait ReadAndSeek: Read+Seek {}
impl <R: Read+Seek> ReadAndSeek for R {}

fn open_file<P: AsRef<Path> + Debug>(path: P) -> Result<Box<dyn ReadAndSeek>>
{
    if path.as_ref().extension().unwrap_or_default() == "gz"
    {
        let mut index_path = Path::new(path.as_ref()).to_owned();
        index_path.set_extension("gz.gzi");
        if !index_path.exists()
        {
            bail!("{} does not exist. bgzip input must be indexed (try bgzip -r {})", index_path.to_str().unwrap_or_default(), path.as_ref().to_str().unwrap_or_default());
        }
        let index = BGZFIndex::from_reader(File::open(index_path)?)?;
        let gzreader = BGZFReader::new(File::open(path)?)?;
        let in_file = IndexedBGZFReader::new(gzreader, index)?;
        Ok(Box::new(in_file))
    }
    else
    {
        let in_file = File::open(path)?;
        Ok(Box::new(in_file))
    }
}

#[allow(non_snake_case)]
pub fn run_bed2pat<P: AsRef<Path> + std::fmt::Debug>(
    fasta_path: P,
    read_bed: P,
    nOT_option: &Option<String>,
    nOB_option: &Option<String>,
    chunk_size_option: &Option<usize>) -> Result<()>
{
    let read_file = open_file(&read_bed)?;
    let read_reader = Reader::new(read_file);
    let fasta_file = open_file(&fasta_path)?;
    let index_path = PathBuf::from(format!("{}.fai", fasta_path.as_ref().to_str().unwrap_or_default()).as_str());
    let fasta_index = Index::from_file(&index_path)?;
    let indexed_reader = IndexedReader::with_index(fasta_file, fasta_index);

    #[allow(non_snake_case)]
    let mut n_ot = ReadMaskSetting::from_str(nOT_option.as_ref().unwrap_or(&"0,0,0,0".to_string())).unwrap();
    let mut n_ob = ReadMaskSetting::from_str(nOB_option.as_ref().unwrap_or(&"0,0,0,0".to_string())).unwrap();

    n_ot.r2 = ReadMask(n_ot.r2.1, n_ot.r2.0); // R2 is mapped in reverse
    n_ob.r1 = ReadMask(n_ob.r1.1, n_ob.r1.0); // F2 is mapped in reverse

    let mut generator = PatGenerator::with_read_and_fasta(read_reader, indexed_reader, n_ot, n_ob, chunk_size_option)?;
    generator.process()
}
/*====================================================
 = Unit Tests
====================================================*/
#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use anyhow::{Ok, Result};

    const COORD_BED: &[u8] = b"chr1\t0\t1\t0\t+
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
        assert_eq!(parsed, [1,2,3]);

        test_str = "[]";
        parsed = parse_mod_str(test_str);
        assert_eq!(parsed.len(), 0);

        test_str = "";
        parsed = parse_mod_str(test_str);
        assert_eq!(parsed.len(), 0);

        test_str = "[1, 2,3]";
        parsed = parse_mod_str(test_str);
        assert_eq!(parsed.len(), 3);
        assert_eq!(parsed, [1,2,3]);

        Ok(())
    }

    // #[test]
    // fn can_zip_mods() -> Result<()>
    // {
    //     let mut buffer = CpgBuffer::with_file(Cursor::new(COORD_BED))?;
    //     let cpg_info = buffer.cpgs_in_range("chr1".as_ref(), 0, 90).expect("Cannot load cpg info");
    //     let mods: Vec<u8> = [1, 4].to_vec();
    //     let unmods: Vec<u8> = [2, 3].to_vec();
    //     let snps: Vec<u8> = Vec::new();
    //     let read_mask = ReadMask(0,0);
    //     assert_eq!(zip_mods(&mods, &unmods, &snps, Strand::Forward, cpg_info, read_mask, 6), [(1, Methylated), (2, Unmethylated), (3,Unmethylated), (4, Methylated)]);

    //     let mods: Vec<u8> = [1].to_vec();
    //     let unmods: Vec<u8> = [2, 3].to_vec();
    //     let snps: Vec<u8> = [4].to_vec();
    //     assert_eq!(zip_mods(&mods, &unmods, &snps, Strand::Forward, cpg_info, read_mask, 6), [(1, Methylated), (2, Unmethylated), (3,Unmethylated), (4, Unknown)]);
    //     Ok(())
    // }


    #[test]
    fn test_read_strand() ->Result<()>
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