use crate::{
    sequence::SegmentsParams,
    utils::TryAsBase as _,
};
use color_eyre::eyre::{Context, Result};
use rust_htslib::bam::{Read as _};
use tracing::{info, instrument, warn};

mod scores;
mod variants;
use variants::{SeenBases, VariantCandidatePileup, pileup_mapper};

#[derive(Debug, clap::Args)]
pub struct CallParams {
    #[command(flatten)]
    segments: SegmentsParams,
}

#[instrument(skip(params))]
pub fn call(params: &CallParams) -> Result<()> {
    let mut segments = params.segments.segments().wrap_err("failed to fetch segments")?;
    let mut seq = Vec::new();

    segments
        .bam
        .pileup()
        .filter_map(|p| p.ok())
        // .filter(|p| fetch_range.contains(&(p.pos() as u64)))
        .take(100)
        .map(|pile| -> Result<Option<VariantCandidatePileup>> {
            segments.fasta.fetch("chr19", u64::from(pile.pos()), u64::from(pile.pos()) + 2)?;
            segments.fasta.read(&mut seq)?;
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
