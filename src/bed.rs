use crate::{
    io::formats::FromFileExtension as _,
    utils::Phred,
    vcf::{GenotypeConfidence, GenotypeLikelihood},
};
use clio::ClioPath;
use color_eyre::{Result, eyre::Context as _};
use rastair2_vcf::standard_fields::{Genotype, GenotypeAllele};
use smol_str::SmolStr;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};
use tracing::{debug, instrument};

#[derive(Debug, Clone, clap::Args)]
pub struct BedParams {
    /// Output BED file with the called methylated positions
    #[arg(long = "bed", required = false, default_missing_value = "-", num_args = 0..=1)]
    pub bed: Option<ClioPath>,

    /// Format of the output BED file
    ///
    /// If not specified, the format is guessed based on the file extension.
    #[arg(long, requires = "bed")]
    pub bed_format: Option<BedFormat>,
}

impl BedParams {
    pub fn bed_format(&self) -> BedFormat {
        if let Some(format) = self.bed_format {
            format
        } else if let Some(path) = &self.bed
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

    pub fn writer(&self) -> Result<Option<BedWriter>> {
        let Some(path) = &self.bed else {
            return Ok(None);
        };

        let format = self.bed_format();
        let writer = BedWriter::new(path, format)
            .wrap_err_with(|| format!("Failed to create BED writer for {path}"))?;
        Ok(Some(writer))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum BedFormat {
    /// BGZIP compressed file, usually `.bed.gz`
    BedGz,
    /// Regular BED file, usually `.bed`
    Bed,
}

pub struct BedWriter {
    pub path: ClioPath,
    pub format: BedFormat,
    writer: Writer,
}

enum Writer {
    Bed(Box<dyn Write + Send + Sync>),
    BedGz(bgzip::BGZFWriter<Box<dyn Write + Send + Sync>>),
}

impl Write for Writer {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Writer::Bed(writer) => writer.write(buf),
            Writer::BedGz(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Writer::Bed(writer) => writer.flush(),
            Writer::BedGz(writer) => writer.flush(),
        }
    }
}

impl BedWriter {
    #[instrument(level = "debug")]
    pub fn new(path: &ClioPath, format: BedFormat) -> Result<Self> {
        let writer = path.clone().create().wrap_err("Failed to create output")?;
        let writer: Box<dyn Write + Send + Sync> = Box::new(BufWriter::new(writer));
        let mut writer: Writer = match format {
            BedFormat::BedGz => {
                let writer = bgzip::BGZFWriter::with_compress_unit_size(
                    writer,
                    bgzip::Compression::fast(),
                    bgzip::write::DEFAULT_COMPRESS_UNIT_SIZE,
                    // Write index if the path is a local file
                    path.is_local(),
                )
                .wrap_err("Failed to create BGZF writer")?;
                Writer::BedGz(writer)
            }
            BedFormat::Bed => Writer::Bed(writer),
        };
        writeln!(&mut writer, "{}", Rastair1BedFormat::HEADER)
            .wrap_err("Failed to write header")?;
        Ok(Self { path: path.clone(), format, writer })
    }

    pub fn write_record(&mut self, record: &Rastair1BedFormat) -> Result<()> {
        record
            .write(&mut self.writer)
            .wrap_err_with(|| format!("Failed to write record {}:{}", record.contig, record.pos))?;
        writeln!(self.writer)?;
        Ok(())
    }

    pub fn close(mut self) -> Result<()> {
        self.writer.flush().wrap_err("Failed to flush writer")?;
        if let Writer::BedGz(bgzfwriter) = self.writer
            && let Some(index) = bgzfwriter.close().wrap_err("Failed to close BGZF writer")?
            && self.path.is_local()
        {
            write_index(self.path.path(), index)
                .wrap_err_with(|| format!("Failed to write index for `{}`", self.path.display()))?;
        }
        Ok(())
    }
}

fn write_index(original_path: &Path, index: bgzip::index::BGZFIndex) -> Result<()> {
    let index_path = original_path.with_extension("gz.gzi");
    let mut index_file = File::create(&index_path)
        .wrap_err_with(|| format!("Failed to create index file `{}`", index_path.display()))?;
    index
        .write(&mut index_file)
        .wrap_err_with(|| format!("Failed to write index to `{}`", index_path.display()))?;
    Ok(())
}

#[derive(Debug, Clone)]
pub struct Rastair1BedFormat {
    pub contig: SmolStr,
    pub pos: usize,
    pub r#ref: SmolStr,
    pub beta: f32,
    pub unmod: u32,
    pub r#mod: u32,
    pub no_snp: u32,
    pub snp: u32,
    pub coverage: usize,
    pub genotype: Genotype,
    pub genotype_likelihood: GenotypeLikelihood,
    pub genotype_confidence: GenotypeConfidence,
}

impl Rastair1BedFormat {
    pub const HEADER: &'static str = "#chr\tstart\tend\tname\tbeta_est\tstrand\tunmod\tmod\tno_snp\tsnp\tcoverage\tgenotype\tgt_p_score\tgt_conf_score";
}

mod internal_to_bed;
mod vcf_to_bed;

impl Rastair1BedFormat {
    fn write(&self, mut f: impl Write) -> Result<()> {
        let Rastair1BedFormat {
            contig,
            pos: start,
            r#ref,
            beta,
            unmod,
            r#mod,
            no_snp,
            snp,
            coverage,
            genotype,
            genotype_likelihood,
            genotype_confidence,
        } = self;
        let end = start + 1;
        let name = ".";
        let strand = match r#ref.as_str() {
            "C" => "+",
            "G" => "-",
            _ => ".",
        };

        write!(
            f,
            "{contig}\t{start}\t{end}\t{name}\t{beta:.2}\t{strand}\t{unmod}\t{mod}\t{no_snp}\t{snp}\t{coverage}\t"
        )?;

        let genotype = genotype_to_rastair1_string(genotype, r#ref);
        let likelihood = genotype_likelihood
            .first()
            .and_then(|x| *x)
            .and_then(|x| Phred::new(x).ok())
            .map(|x| *x)
            .unwrap_or_default();
        let confidence = genotype_confidence
            .first()
            .and_then(|x| *x)
            .and_then(|x| Phred::new(x).ok())
            .map(|x| *x)
            .unwrap_or_default();
        write!(f, "{genotype}\t{likelihood:.2}\t{confidence:.2}")?;

        Ok(())
    }
}

pub fn genotype_to_rastair1_string(genotype: &Genotype, ref_base: &str) -> SmolStr {
    match genotype.0.as_slice() {
        [GenotypeAllele::Phased(0)] | [GenotypeAllele::Unphased(0)] => {
            if ref_base == "C" {
                SmolStr::new_static("C/C")
            } else {
                SmolStr::new_static("G/G")
            }
        }
        [GenotypeAllele::Phased(0), GenotypeAllele::Phased(1)]
        | [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)] => {
            if ref_base == "C" {
                SmolStr::new_static("C/T")
            } else {
                SmolStr::new_static("G/A")
            }
        }
        [GenotypeAllele::Phased(1)] | [GenotypeAllele::Unphased(1)] => {
            if ref_base == "C" {
                SmolStr::new_static("T/T")
            } else {
                SmolStr::new_static("A/A")
            }
        }
        _ => SmolStr::new_static("."),
    }
}
