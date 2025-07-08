use crate::{
    utils::{Base::*, Phred},
    vcf::{
        GenotypeConfidence, GenotypeLikelihood, Methylated, Record as Rastair2Record,
        utils::NoStrandBiasForBaseErrorExt as _,
    },
};
use clio::ClioPath;
use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat as _},
};
use rastair2_vcf::{
    VcfField,
    standard_fields::{Genotype, GenotypeAllele},
};
use rust_htslib::bcf::Record as HtslibRecord;
use smol_str::SmolStr;
use std::io::Write;
use tracing::instrument;

#[derive(Debug)]
pub struct BedWriter {
    pub path: ClioPath,
    writer: clio::Output,
}

impl BedWriter {
    #[instrument(level = "debug")]
    pub fn new(path: &ClioPath) -> Result<Self> {
        let mut writer =
            path.clone().create().wrap_err_with(|| format!("Failed to create output {path}"))?;
        writeln!(&mut writer, "{}", Rastair1BedFormat::HEADER)
            .wrap_err_with(|| format!("Failed to write header to {path}"))?;
        Ok(Self { path: path.clone(), writer })
    }

    pub fn write_record(&mut self, record: &Rastair1BedFormat) -> Result<()> {
        record.write(&mut self.writer).wrap_err_with(|| {
            format!("Failed to write record for {}:{}", record.contig, record.pos)
        })?;
        writeln!(self.writer)?;
        Ok(())
    }
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

impl TryFrom<&Rastair2Record> for Rastair1BedFormat {
    type Error = color_eyre::eyre::Report;

    #[allow(clippy::cast_possible_truncation)]
    fn try_from(record: &Rastair2Record) -> Result<Self, Self::Error> {
        let r#ref = record.main.r#ref.clone();

        let (unmod, r#mod, no_snp, snp) = if r#ref == "C" {
            (
                record.strand_count(C).or_empty().ot,
                record.strand_count(T).or_empty().ot,
                record.strand_count(C).or_empty().ob,
                record.strand_count(T).or_empty().ob,
            )
        } else if r#ref == "G" {
            (
                record.strand_count(G).or_empty().ob,
                record.strand_count(A).or_empty().ob,
                record.strand_count(G).or_empty().ot,
                record.strand_count(A).or_empty().ot,
            )
        } else {
            (0, 0, 0, 0)
        };

        Ok(Rastair1BedFormat {
            contig: record.main.chrom.clone(),
            pos: record.main.pos as usize,
            r#ref,
            beta: record.samples[0].methylated.beta().unwrap_or_default() as f32,
            unmod,
            r#mod,
            no_snp,
            snp,
            coverage: *record.info.read_depth,
            genotype: record.samples[0].genotype.clone(),
            genotype_likelihood: record.samples[0].genotype_likelihood.clone(),
            genotype_confidence: record.samples[0].genotype_confidence.clone(),
        })
    }
}

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
            "{contig}\t{start}\t{end}\t{name}\t{beta}\t{strand}\t{unmod}\t{mod}\t{no_snp}\t{snp}\t{coverage}\t"
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
        write!(f, "{genotype}\t{likelihood}\t{confidence}")?;

        Ok(())
    }
}

fn genotype_to_rastair1_string(genotype: &Genotype, ref_base: &str) -> SmolStr {
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
