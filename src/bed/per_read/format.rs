use crate::{bed::BedRecord, sequence::Region};
use color_eyre::eyre::Result;
use rastair_types::SmallVec;
use std::io::Write;

/// Store methylation information for a single read
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PerRead {
    /// ID of the sequence
    pub region: Region,
    /// Flag of read
    pub flag: u16,
    /// Mapq of read
    pub mapq: u8,
    /// Absolute fragment length (non-directional)
    pub frag_length: u64,
    /// Read length
    pub read_length: usize,
    /// Name of read
    pub read_id: String,
    /// Number of CpGs in a read
    pub cpg_count: u16,
    /// Number of modified CpGs
    pub mod_count: usize,
    /// Positions in read of modified CpGs
    pub mod_cpgs: SmallVec<usize, 24>,
    /// Positions in read of unmodified CpGs
    pub unmod_cpgs: SmallVec<usize, 24>,
    /// Positions in read of CpGs that are mutated
    pub snp_cpgs: SmallVec<usize, 24>,
    /// Positions in read of de-novo CpGs that are mutated
    pub mod_denovos: SmallVec<usize, 10>,
    /// Positions in read of de-novo CpGs that are unmodified
    pub unmod_denovos: SmallVec<usize, 10>,
}

impl BedRecord for PerRead {
    const HEADER: &'static str = "#chr\tstart\tend\tread_id\tmapq\torientation\tinsert_size\tread_length\tflag\tnum_cpg\tnum_mod\tmod_cpgs\tunmod_cpgs\tsnp_cpgs\tmod_denovos\tunmod_denovos";

    fn write<W: Write>(&self, f: &mut W) -> Result<()> {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.region.contig,
            self.region.start,
            self.region.end,
            self.read_id,
            self.mapq,
            if self.flag & 16 == 16 { "-" } else { "+" },
            self.frag_length,
            self.read_length,
            self.flag,
            self.cpg_count,
            self.mod_count,
        )?;

        write!(f, "\t")?;
        write_list(f, &self.mod_cpgs)?;
        write!(f, "\t")?;
        write_list(f, &self.unmod_cpgs)?;

        write!(f, "\t")?;
        write_list(f, &self.snp_cpgs)?;

        write!(f, "\t")?;
        write_list(f, &self.mod_denovos)?;
        write!(f, "\t")?;
        write_list(f, &self.unmod_denovos)?;
        writeln!(f)?;

        Ok(())
    }

    fn chr(&self) -> &str {
        &self.region.contig
    }

    fn start(&self) -> usize {
        usize::try_from(self.region.start).expect("region start should fit into usize")
    }

    fn end(&self) -> usize {
        usize::try_from(self.region.end).expect("region end should fit into usize")
    }
}

fn write_list<W: Write>(f: &mut W, list: &[usize]) -> Result<()> {
    for (i, pos) in list.iter().enumerate() {
        if i > 0 {
            write!(f, ",")?;
        }
        write!(f, "{pos}")?;
    }
    Ok(())
}
