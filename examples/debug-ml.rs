use clap::Parser as _;
use clio::ClioPath;
use color_eyre::{Result, eyre::Context as _};
use ndarray::ArrayView1;
use rastair::{
    call::{
        CallParams,
        ml::{self, MlModel, Prediction},
        process,
    },
    sequence::ChunkRegion,
    utils::{logging::setup_logging, surrounding_records},
};
use std::io::{BufWriter, Write};

#[derive(Debug, clap::Parser)]
struct Cli {
    #[clap(flatten)]
    call: CallParams,

    /// Enable more logging
    #[arg(short, long, global = true)]
    verbose: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    setup_logging(args.verbose);
    let params = args.call;

    let mut readers = params.segments.readers().wrap_err("Failed to read BAM/FASTA files")?;
    let regions: Vec<ChunkRegion> = readers
        .segments(100_000, 100)
        .wrap_err("Could not fetch segments from BAM file")?
        .collect();
    let ml = params.ml.init().wrap_err("Failed to initialize machine learning model")?;

    let mut cpg = MlDebugOutput::new(ClioPath::new("CpG_features.tsv.gz")?)?;
    writeln!(&mut cpg.writer, "{CPG_HEADER}")?;
    let mut denovo = MlDebugOutput::new(ClioPath::new("denovo_features.tsv.gz")?)?;
    writeln!(&mut denovo.writer, "{DENOVO_HEADER}")?;
    let mut other = MlDebugOutput::new(ClioPath::new("other_features.tsv.gz")?)?;
    writeln!(&mut other.writer, "{OTHER_HEADER}")?;

    for region in regions {
        let pileup_mapping_params = process::PileupMappingParams {
            include_cpgs: params.methylation.should_include_all_cpgs(),
            variant_calling: params.variant_calling.clone(),
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

            if let ml::MlResult::Predictions(predictions) = ml.predict(current, before, after) {
                for Prediction { model, prediction, features, .. } in predictions {
                    match model {
                        MlModel::Cpg => &mut cpg,
                        MlModel::DenovoCpg => &mut denovo,
                        MlModel::Others => &mut other,
                    }
                    .write(
                        &current.main.chrom,
                        current.main.pos,
                        *prediction,
                        features.view(),
                    )?;
                }
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

// Keep these in sync with the `params_from_record` functions!
const CPG_HEADER: &str = "chr\tpos\tresult\tad_alt_adj\talt_score_adj\tref_a\tref_c\tref_g\tref_t\talt_a\talt_c\talt_g\talt_t\tmapq\tnum_mapq0\tp1a\tp1c\tp1g\tp1t\tp2a\tp2c\tp2g\tp2t\tp4a\tp4c\tp4g\tp4t\tp5a\tp5c\tp5g\tp5t\tregion_entropy\tad_ref\tad_alt\tsb_ot_ref\tsb_ob_ref\tsb_ot_alt\tsb_ob_alt\talt_score\tbq_ref\tbq_alt\tbq_ot_ref\tbq_ob_ref\tbq_ot_alt\tbq_ob_alt\tmq_ot_ref\tmq_ob_ref\tmq_ot_alt\tmq_ob_alt\tmq_ref\tmq_alt\tpos_in_read_ref\tpos_in_read_alt\tnum_aligned_bases_ref\tnum_aligned_bases_alt\tnum_indels_ref\tnum_indels_alt\tbeta_ratio";
const DENOVO_HEADER: &str = "chr\tpos\tresult\tad_alt_adj\talt_score_adj\tsb_adj\tref_a\tref_c\tref_g\tref_t\talt_a\talt_c\talt_g\talt_t\tmapq\tnum_mapq0\tp1a\tp1c\tp1g\tp1t\tp2a\tp2c\tp2g\tp2t\tp4a\tp4c\tp4g\tp4t\tp5a\tp5c\tp5g\tp5t\tregion_entropy\tad_ref\tad_alt\tsb_ot_ref\tsb_ob_ref\tsb_ot_alt\tsb_ob_alt\talt_score\tbq_ref\tbq_alt\tbq_ot_ref\tbq_ob_ref\tbq_ot_alt\tbq_ob_alt\tmq_ot_ref\tmq_ob_ref\tmq_ot_alt\tmq_ob_alt\tmq_ref\tmq_alt\tpos_in_read_ref\tpos_in_read_alt\tnum_aligned_bases_ref\tnum_aligned_bases_alt\tnum_indels_ref\tnum_indels_alt\tbeta_ratio";
const OTHER_HEADER: &str = "chr\tpos\tresult\tref_a\tref_c\tref_g\tref_t\talt_a\talt_c\talt_g\talt_t\t*mapq\tnum_mapq0\tp1a\tp1c\tp1g\tp1t\tp2a\tp2c\tp2g\tp2t\tp4a\tp4c\tp4g\tp4t\tp5a\tp5c\tp5g\tp5t\tregion_entropy\tad_ref\tad_alt\tsb_ot_ref\tsb_ob_ref\tsb_ot_alt\tsb_ob_alt\tsb_alt\tsb_ref\talt_score\tbq_ref\tbq_alt\tbq_ot_ref\tbq_ob_ref\tbq_ot_alt\tbq_ob_alt\tmq_ot_ref\tmq_ob_ref\tmq_ot_alt\tmq_ob_alt\tmq_ref\tmq_alt\tpos_in_read_ref\tpos_in_read_alt\tnum_aligned_bases_ref\tnum_aligned_bases_alt\tnum_indels_ref\tnum_indels_alt";
