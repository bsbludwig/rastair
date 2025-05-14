use crate::{
    sequence::{ChunkRegion, Segment, SegmentsParams},
    utils::TryAsBase as _,
};
use clio::ClioPath;
use color_eyre::eyre::{Context, ContextCompat, Result};
use rust_htslib::bam::{FetchDefinition, Read as _, pileup::Pileup};
use scores::VariantCandidatePileupMetrics;
use std::io::{self, BufWriter};
use tracing::{instrument, warn};

mod methylation;
mod scores;
mod variants;
mod filtering {
    pub mod threshold;
}
use variants::{SeenBases, VariantCandidatePileup, pileup_mapper};

#[derive(Debug, clap::Args)]
pub struct CallParams {
    #[command(flatten)]
    segments: SegmentsParams,

    #[command(flatten)]
    thresholds: filtering::threshold::ThresholdConfig,

    #[arg(short = 'o', long)]
    vcf_output: ClioPath,
}

#[instrument(level = "debug", skip(params))]
pub fn call(params: &CallParams) -> Result<()> {
    let mut readers = params.segments.readers().wrap_err("failed to fetch segments")?;

    let mut output = BufWriter::new(
        params.vcf_output.clone().create().wrap_err("failed to create output file/stream")?,
    );
    MethylationEventWriter::write_header(&mut output)?;

    let regions = readers.segments().wrap_err("could not fetch segments")?;
    for region in regions {
        let segment = readers.segment(&region)?;

        readers.bam.fetch(
            FetchDefinition::try_from(&segment)
                .wrap_err("convert region string to fetch definition")?,
        )?;
        readers
            .bam
            .pileup()
            .filter_map(|p| p.ok())
            .filter(|p| {
                // Filter out pileups that are not in the region of interest
                region.contains(u64::from(p.pos()))
            })
            .flat_map(|pile| collect_candidate(pile, &segment, &region).transpose())
            .map(|pile| -> Result<_> {
                let pile = pile?;
                let metrics = pile.metrics().wrap_err("calculate metrics")?;
                Ok((pile, metrics))
            })
            .peekable()
            // TODO: also peek next pile,metrics if available
            .try_for_each(|x| -> Result<()> {
                match x {
                    Ok((pile, metrics)) => {
                        if pile.likely_methylation_event(&metrics, &params.thresholds) {
                            MethylationEventWriter(&pile, &metrics).write(&mut output)?;
                        }
                    }
                    Err(error) => {
                        warn!(%error, "failed to get pileup");
                    }
                }
                Ok(())
            })?;
    }

    return Ok(());
}

#[instrument(level = "debug", skip_all)]
fn collect_candidate(
    pile: Pileup,
    segment: &Segment,
    region: &ChunkRegion,
) -> Result<Option<VariantCandidatePileup>> {
    let idx = pile.pos() as usize
        - usize::try_from(segment.range.start).wrap_err("segment range fits in usize")?;
    let bases = SeenBases(pile.alignments().filter_map(pileup_mapper).collect());
    let reference_base =
        segment.sequence.get(idx).wrap_err("failed to get reference base")?.as_base()?;
    let next_base = segment.sequence.get(idx + 1).and_then(|x| x.as_base().ok());
    if bases.is_variant_candidate() {
        // info!(?bases, pos = pile.pos(), ?reference_base, ?next_base, "found pile of interest");
        Ok(Some(VariantCandidatePileup {
            chrom: region.chromosome.clone(),
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
}

struct MethylationEventWriter<'p, 'm>(
    &'p VariantCandidatePileup,
    &'m VariantCandidatePileupMetrics,
);

impl MethylationEventWriter<'_, '_> {
    fn write_header(mut w: impl io::Write) -> Result<()> {
        write!(w, "#")?;
        ["CHROM", "POS", "REF", "ALT", "VAF", "BINOM", "BETA"]
            .iter()
            .try_for_each(|x| write!(w, "{}\t", x))?;
        writeln!(w)?;
        Ok(())
    }

    fn write(&self, mut w: impl io::Write) -> Result<()> {
        let chrom = self.0.chrom.as_str();
        let pos = self.0.pos;
        let r#ref = self.1.reference_count;
        let alt = self.1.alt_count;
        let vaf = self.1.vaf;
        let binom = self.1.binomial;
        let beta = self.0.beta();

        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}", chrom, pos, r#ref, alt, vaf, binom, beta)?;

        Ok(())
    }
}
