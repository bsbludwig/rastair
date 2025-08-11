use crate::{
    bed::{BedFormat, BedRecord, writer::BedWriter},
    io::formats::FromFileExtension as _,
    sequence::Region,
};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, Result};
use std::io::Write;
use tracing::debug;

#[derive(Debug, Clone, clap::Args)]
pub struct BedReadsParams {
    /// Output BED file with all reads
    #[arg(long, required = false, default_missing_value = "-", num_args = 0..=1)]
    pub bed_reads: Option<ClioPath>,

    /// Format of the output BED reads file
    ///
    /// If not specified, the format is guessed based on the file extension.
    #[arg(long, requires = "bed_reads")]
    pub bed_reads_format: Option<BedFormat>,
}

impl BedReadsParams {
    pub fn bed_format(&self) -> BedFormat {
        if let Some(format) = self.bed_reads_format {
            format
        } else if let Some(path) = &self.bed_reads
            && let Some(path) = path.path().to_str()
            && let Some(format) = BedFormat::from_file_extension(path)
        {
            format
        } else {
            debug!(
                "Could not determine BED output format from file extension, defaulting to uncompressed"
            );
            BedFormat::Bed
        }
    }

    pub fn writer(&self) -> Result<Option<BedWriter<PerRead>>> {
        let Some(path) = &self.bed_reads else {
            return Ok(None);
        };

        let format = self.bed_format();
        let writer = BedWriter::new(path, format)
            .wrap_err_with(|| format!("Failed to create BED writer for {path}"))?;
        Ok(Some(writer))
    }
}

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
    pub frag_length: u32,
    /// Read length
    pub read_length: u32,
    /// Name of read
    pub read_id: String,
    /// Number of CpGs in a read
    pub cpg_count: u16,
    /// Number of modified CpGs
    pub mod_count: u16,
    /// Positions in read  of modified CpGs
    pub mod_cpgs: Vec<usize>,
    /// Positions in read of unmodified CpGs
    pub unmod_cpgs: Vec<usize>,
    /// Positions in read of CpGs that are mutated
    pub snp_cpgs: Vec<usize>,
}

impl BedRecord for PerRead {
    const HEADER: &'static str = "#chr\tstart\tend\tread_id\tmapq\torientation\tinsert_size\tread_length\tflag\tnum_cpg\tnum_mod\tmod_cpgs\tunmod_cpgs\tsnp_cpgs";

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
        for (i, pos) in self.mod_cpgs.iter().enumerate() {
            if i > 0 {
                write!(f, ",")?;
            }
            write!(f, "{pos}")?;
        }

        write!(f, "\t")?;
        for (i, pos) in self.unmod_cpgs.iter().enumerate() {
            if i > 0 {
                write!(f, ",")?;
            }
            write!(f, "{pos}")?;
        }

        write!(f, "\t")?;
        for (i, pos) in self.snp_cpgs.iter().enumerate() {
            if i > 0 {
                write!(f, ",")?;
            }
            write!(f, "{pos}")?;
        }

        writeln!(f)?;

        Ok(())
    }
}
