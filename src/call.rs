use crate::utils::{RegionString, TryAsBase as _, file_helpers::open_fasta};
use clap::value_parser;
use clio::ClioPath;
use color_eyre::eyre::Result;
use rust_htslib::bam::{self, FetchDefinition, Read as _};
use tracing::{info, instrument, warn};

mod scores;
mod variants;
use variants::{SeenBases, VariantCandidatePileup, pileup_mapper};

#[derive(Debug, clap::Args)]
pub struct CallParams {
    /// A sorted and indexed bam file
    #[arg(value_name="BAM_FILE", value_parser=value_parser!(ClioPath).exists().is_file())]
    bam_file: ClioPath,

    /// A sorted and indexed (via samtools faidx) fasta file. Can be bgzip
    /// compressed, but requires both a gzi index and a fai index
    #[arg(short='r', long, value_name="FASTA_FILE", required=true, value_parser=value_parser!(ClioPath).exists().is_file())]
    fasta_file: ClioPath,

    /// Restrict to a specific chromosome or region of a chromosome. Format is
    /// "chr", "chr:start" or "chr:start-end", where start is 1-based and end is
    /// inclusive.
    #[arg(short = 'l', long)]
    region: Option<RegionString>,
}

#[instrument(skip(params))]
pub fn read(params: &CallParams) -> Result<()> {
    let mut fasta = open_fasta(&params.fasta_file)?;
    // indexed_reader.fetch("chr19", fetch_range.start, fetch_range.end + 1)?;

    let mut bam = bam::IndexedReader::from_path(params.bam_file.path())?;
    bam.set_threads(8)?;
    if let Some(region) = &params.region {
        bam.fetch(region)?;
    } else {
        bam.fetch(FetchDefinition::All)?;
    }
    // bam.fetch(("chr19", fetch_range.start, fetch_range.end + 1))?;

    let mut seq = Vec::new();

    bam.pileup()
        .filter_map(|p| p.ok())
        // .filter(|p| fetch_range.contains(&(p.pos() as u64)))
        .take(100)
        .map(|pile| -> Result<Option<VariantCandidatePileup>> {
            fasta.fetch("chr19", pile.pos() as u64, pile.pos() as u64 + 2)?;
            fasta.read(&mut seq)?;
            let bases = SeenBases(pile.alignments().filter_map(pileup_mapper).collect());
            let reference_base = seq[0].as_base()?;
            let next_base = seq.get(1).and_then(|x| x.as_base().ok());
            if bases.is_variant_candidate() {
                // info!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "found pile of interest");
                Ok(Some(VariantCandidatePileup {
                    pos: pile.pos(),
                    bases,
                    reference_base,
                    next_base,
                }))
                // info!(?pileup, metrics=?pileup.metrics(), "variant candidate");
            } else if bases.matches(reference_base) {
                // Matches reference base
                // boring.
                // trace!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "pile matches reference");
                Ok(None)
            } else {
                warn!(
                    ?bases,
                    pos = pile.pos(),
                    ?reference_base,
                    ?next_base,
                    "pile does not match reference but is also not interesting"
                );
                Ok(None)
            }
        })
        .flat_map(
            |x| -> Option<Result<(VariantCandidatePileup, scores::VariantCandidatePileupMetrics)>> {
                let x = x.transpose()?;
                Some(x.map(|pile| {
                    let metrics = pile.metrics();
                    (pile, metrics)
                }))
            },
        )
        .for_each(|x| {
            if let Ok((pile, metrics)) = x {
                info!(?pile, metrics=?metrics, "variant candidate");
            } else {
                warn!("failed to get pileup");
            }
        });

    return Ok(());
}
