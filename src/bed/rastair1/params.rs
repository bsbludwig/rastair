use crate::{
    bed::{BedFormat, writer::BedWriter},
    io::formats::FromFileExtension as _,
    utils::cli,
};
use better_default::Default;
use clio::ClioPath;
use color_eyre::{Result, eyre::Context as _};
use seqair_types::Probability;
use tracing::{debug, instrument};

#[derive(Debug, Clone, clap::Args, Default)]
pub struct BedParams {
    /// Output BED file with the called methylated positions
    #[arg(long = "bed", required = false, default_missing_value = "-", num_args = 0..=1)]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub bed: Option<ClioPath>,

    /// Format of the output BED file
    ///
    /// If not specified, the format is guessed based on the file extension.
    #[arg(long, requires = "bed")]
    #[arg(help_heading = cli::sections::OUTPUT)]
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

    #[instrument(level = "debug", skip(self))]
    pub fn writer(&self) -> Result<Option<BedWriter<super::Rastair1BedFormat>>> {
        let Some(path) = &self.bed else {
            return Ok(None);
        };

        let format = self.bed_format();
        let writer = BedWriter::new(path, format)
            .wrap_err_with(|| format!("Failed to create BED writer for {path}"))?;
        Ok(Some(writer))
    }
}

/// Parameters for filtering BED records produced by `convert`
#[derive(Debug, Clone, clap::Args, Default)]
pub struct BedRecordsFilterParams {
    /// Include CpG positions with zero coverage
    ///
    /// This can be useful to get a complete list of CpG positions in the output BED file.
    /// Note that this requires the input data to contain a complete list of CpG positions,
    /// e.g. by using the `--all --cpgs-only` options when calling methylation.
    #[arg(long = "bed-include-empty")]
    #[arg(help_heading = cli::sections::FILTER)]
    pub include_empty: bool,
}

#[derive(Debug, Clone)]
pub struct BedRecordsConvertParams {
    /// Minimum ML score to consider a position as variant
    pub ml_threshold: Probability,
    pub filters: BedRecordsFilterParams,
}
