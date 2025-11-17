use crate::{
    bed::{BedFormat, writer::BedWriter},
    io::formats::FromFileExtension as _,
    utils::cli,
};
use clio::ClioPath;
use color_eyre::eyre::{Context as _, Result};
use tracing::{debug, instrument};

#[derive(Debug, Clone, clap::Args)]
pub struct BedReadsParams {
    /// Output BED file with all reads
    #[arg(long, required = false, default_value = "-", num_args = 0..=1)]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub bed: ClioPath,

    /// Format of the output BED reads file
    ///
    /// If not specified, the format is guessed based on the file extension.
    #[arg(long, requires = "bed")]
    #[arg(help_heading = cli::sections::OUTPUT)]
    pub bed_format: Option<BedFormat>,
}

impl BedReadsParams {
    pub fn bed_format(&self) -> BedFormat {
        if let Some(format) = self.bed_format {
            format
        } else if let Some(path) = self.bed.path().to_str()
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
    pub fn writer(&self) -> Result<BedWriter<super::format::PerRead>> {
        let path = &self.bed;

        let format = self.bed_format();
        let writer = BedWriter::new(path, format)
            .wrap_err_with(|| format!("Failed to create BED writer for {path}"))?;
        Ok(writer)
    }
}
