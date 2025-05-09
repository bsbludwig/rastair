use crate::{sequence::SegmentsParams, utils::TryAsBase as _};
use clio::ClioPath;
use color_eyre::eyre;
use color_eyre::eyre::{Context, Result};
use rust_htslib::bam::Read as _;
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
    let chr = "chr19";
    let mut readers = params.segments.readers().wrap_err("failed to fetch segments")?;

    let mut output = BufWriter::new(
        params.vcf_output.clone().create().wrap_err("failed to create output file/stream")?,
    );
    MethylationEventWriter::write_header(&mut output)?;

    let regions = readers.segments()?;
    for region in regions {
        let segment = readers.segment(&region)?;

        readers.bam.fetch(&segment)?;
        readers
            .bam
            .pileup()
            .filter_map(|p| p.ok())
            .filter(|p| {
                // Filter out pileups that are not in the region of interest
                region.contains(u64::from(p.pos()))
            })
            .map(|pile| -> Result<Option<VariantCandidatePileup>> {
                let idx = pile.pos() as usize
                    - usize::try_from(segment.range.start).expect("segment range fits in usize");
                let bases = SeenBases(pile.alignments().filter_map(pileup_mapper).collect());
                let reference_base = segment.sequence[idx].as_base()?;
                let next_base = segment.sequence.get(idx + 1).and_then(|x| x.as_base().ok());
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
            .flat_map(Result::transpose)
            .map(|x| {
                x.map(|pile| {
                    let metrics = pile.metrics();
                    (pile, metrics)
                })
            })
            // TODO: also peek next pile,metrics if available
            .try_for_each(|x| -> Result<()> {
                match x {
                    Ok((pile, metrics)) => {
                        // trace!(?pile, ?metrics, "found variant candidate");
                        if pile.likely_methylation_event(&metrics, &params.thresholds) {
                            // let bases = pile.bases.iter().fold(String::new(), |mut acc, b| {
                            //     write!(&mut acc, "{}", b.base.display_colored()).unwrap();
                            //     acc
                            // });
                            // info!(
                            //     pos=pile.pos,
                            //     ref=%pile.reference_base.display_colored(),
                            //     %bases,
                            //     ?metrics,
                            //     "found methylation event"
                            // );
                            MethylationEventWriter(pile, metrics).write(chr, &mut output)?;
                        }
                    }
                    Err(_error) => {
                        // warn!(%error, "failed to get pileup");
                        // Err(error)
                    }
                }
                Ok(())
            })?;
    }

    return Ok(());
}

struct MethylationEventWriter(VariantCandidatePileup, VariantCandidatePileupMetrics);

impl MethylationEventWriter {
    fn write_header(mut w: impl io::Write) -> Result<()> {
        write!(w, "#")?;
        ["CHROM", "POS", "REF", "ALT", "VAF", "BINOM", "BETA"]
            .iter()
            .try_for_each(|x| write!(w, "{}\t", x))?;
        writeln!(w)?;
        Ok(())
    }

    fn write(&self, chrom: &str, mut w: impl io::Write) -> Result<()> {
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
