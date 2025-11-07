use crate::{
    call::{
        methylation::params::MethylationCallingParams, pileup::Pileup,
        variant_calling::VariantCallingParams,
    },
    metrics::{self, PileupMetrics, PositionMetricsExt},
    sequence::{ChunkRegion, Readers, Segment},
    utils::logging::ThisIsABug,
};
use color_eyre::eyre::{Result, WrapErr};
use rust_htslib::bam::{FetchDefinition, Read as _};
use std::{ops::Deref, rc::Rc};
use tracing::{Level, debug, instrument, trace, warn};

pub struct PileupMappingParams {
    pub variant_calling: VariantCallingParams,
}

impl Deref for PileupMappingParams {
    type Target = VariantCallingParams;

    fn deref(&self) -> &Self::Target {
        &self.variant_calling
    }
}

#[instrument(level = "info", skip_all)]
pub fn get_pileups(
    readers: &mut Readers,
    region: &ChunkRegion,
    params: &PileupMappingParams,
) -> Result<(Rc<Segment>, impl Iterator<Item = Pileup>)> {
    let segment = readers.segment(region, 2).wrap_err("Failed to fetch segment")?;
    debug!(len = segment.sequence.len(), "Processing region");

    // Fetch the pileups for the segment
    FetchDefinition::try_from(&segment.region)
        .wrap_err("Could not convert region string")
        .this_is_a_bug()
        .and_then(|r| readers.bam.fetch(r).wrap_err("Could not fetch segment from BAM file"))
        .wrap_err_with(|| format!("Could not fetch region `{}` from BAM file", region.region))?;

    let segment = Rc::new(segment);
    let segment_clone = segment.clone();

    // Go over each column in the pileup from htslib and build our own pileup
    let piles = readers
        .bam
        .pileup()
        .filter_map(|p| match p {
            Ok(p) => Some(p),
            Err(e) => {
                if tracing::enabled!(Level::TRACE) {
                    trace!(%e, "Failed to read pileup, skipping");
                }
                None
            }
        })
        .filter(|p| {
            // Sometimes we get pileups that are not in the region of interest.
            // TODO: figure out why and when this happens
            // TODO: Write test for this
            region.contains(u64::from(p.pos()))
        })
        .map(move |pile| {
            Pileup::from_hts(&pile, segment.clone(), params).wrap_err_with(|| {
                format!("Failed to get candidate from pileup at position {}", pile.pos())
            })
        })
        .filter_map(|res| match res {
            Ok(x) => Some(x),
            Err(error) => {
                warn!(error = format!("{error:#}"), "Failed to get pileup, skipping");
                None
            }
        });
    Ok((segment_clone, piles))
}

pub struct PileupMetricsParams {
    pub variant_calling: VariantCallingParams,
    pub methylation: MethylationCallingParams,
}

#[instrument(level = "info", skip_all)]
pub fn calculate_pileup_metrics(
    pileups: impl Iterator<Item = Pileup>,
    segment: &Segment,
    params: &PileupMetricsParams,
) -> impl Iterator<Item = Result<PileupMetrics>> {
    pileups.into_iter().map(PileupMetrics::new).map(move |metrics| {
        // Set "extended" metrics that depend on the pileup and some external params
        let mut current = metrics?;

        let genotype = current.pileup.estimate_genotype(params.variant_calling.error_model);
        let methylated = metrics::methylation::call(&params.methylation.thresholds, &current)?
            .unwrap_or_default();

        let region_entropy = segment
            .entropy_around::<100>(current.pileup.idx())
            .wrap_err("Failed to calculate region entropy")?;

        let ext = PositionMetricsExt {
            genotype,
            methylated,
            region_entropy,
            denovo_adj: metrics::DenovoAdjecent::No,
        };
        current.set_extended_metrics(ext);

        Ok(current)
    })
}
