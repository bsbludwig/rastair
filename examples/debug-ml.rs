use std::io::{BufWriter, Write};

use clap::Parser as _;
use clio::ClioPath;
use color_eyre::{Result, eyre::Context as _};
use ndarray::ArrayView1;
use rastair2::{
    call::{
        CallParams,
        ml::{self, MlModel},
        process,
    },
    sequence::ChunkRegion,
    utils::surrounding_records,
};

#[derive(Debug, clap::Parser)]
struct Cli {
    #[clap(flatten)]
    call: CallParams,

    /// Enable more logging
    #[arg(short, long, global = true)]
    verbose: bool,
}

fn main() -> Result<()> {
    color_eyre::install()?;

    let args = Cli::parse();
    rastair2::utils::setup_tracing(args.verbose);
    let mut params = args.call;

    if params.ml.ml.is_none() {
        // Always do ML, that's the point of this command
        params.ml.ml = Some(0.8);
    }

    let mut readers = params.segments.readers().wrap_err("Failed to read BAM/FASTA files")?;
    let regions: Vec<ChunkRegion> =
        readers.segments().wrap_err("Could not fetch segments from BAM file")?.collect();
    let ml = params.ml.init().wrap_err("Failed to initialize machine learning model")?;

    let mut cpg = MlDebugOutput::new(ClioPath::new("CpG_features.tsv.gz")?)?;
    let mut denovo = MlDebugOutput::new(ClioPath::new("denovo_features.tsv.gz")?)?;
    let mut other = MlDebugOutput::new(ClioPath::new("other_features.tsv.gz")?)?;

    for region in regions {
        let pileup_mapping_params = process::PileupMappingParams {
            include_cpgs: params.methylation.should_include_all_cpgs(),
            keep_overlapping_reads: params.variant_calling.keep_overlapping_reads,
            read_masking: params.variant_calling.read_masking.clone(),
            read_flags: params.variant_calling.read_flags.clone(),
        };

        let piles = region.process(&mut readers, &pileup_mapping_params)?;

        let mut records = piles
            .iter()
            .map(|pile| pile.variant_metrics(&params.variant_calling))
            .collect::<Result<Vec<_>>>()
            .wrap_err("Failed to collect metrics")?;

        // Call methylation events if requested
        let record_len = records.len();
        for i in 0..record_len {
            let (before, current, after) = surrounding_records(&mut records, i);

            params
                .methylation
                .call(current, before, after)
                .wrap_err("Failed to call methylation")?;
        }

        // Filter out piles that are not CpG if requested
        if params.variant_calling.cpgs_only {
            records.retain(|record| *record.info.in_cp_g || *record.info.de_novo_cp_g_candidate);
        }

        let record_len = records.len();
        for i in 0..record_len {
            let (before, current, after) = surrounding_records(&mut records, i);

            if let ml::MlResult::Prediction { model, prediction, features, .. } =
                ml.predict(current, before, after)
            {
                match model {
                    MlModel::Cpg => &mut cpg,
                    MlModel::DenovoCpg => &mut denovo,
                    MlModel::Others => &mut other,
                }
                .write(
                    &current.main.chrom,
                    current.main.pos,
                    prediction,
                    features.view(),
                )?;
            }
        }
    }

    Ok(())
}

pub struct MlDebugOutput {
    writer: Box<dyn Write>,
}

impl MlDebugOutput {
    pub fn new(path: ClioPath) -> Result<Self> {
        let file = path.create().wrap_err("Failed to create debug output file")?;
        let compressor = bgzip::BGZFWriter::new(file, bgzip::Compression::fast());
        let writer: Box<dyn Write> = Box::new(BufWriter::new(compressor));
        Ok(MlDebugOutput { writer })
    }

    pub fn write(
        &mut self,
        chr: &str,
        pos: u32,
        result: f64,
        features: ArrayView1<f64>,
    ) -> Result<()> {
        write!(self.writer, "{chr}\t{pos}\t{result}\t")?;
        for (i, value) in features.iter().enumerate() {
            if i > 0 {
                write!(self.writer, "\t")?;
            }
            write!(self.writer, "{value}")?;
        }
        writeln!(self.writer)?;
        Ok(())
    }
}
