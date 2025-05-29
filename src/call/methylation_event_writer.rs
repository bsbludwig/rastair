use super::{scores::VariantCandidatePileupMetrics, variants::VariantCandidatePileup};
use color_eyre::eyre::Result;
use std::io::{self};

pub struct MethylationEventWriter<'p, 'm>(
    pub &'p VariantCandidatePileup,
    pub &'m VariantCandidatePileupMetrics,
);

impl MethylationEventWriter<'_, '_> {
    pub fn write_header(mut w: impl io::Write) -> Result<()> {
        write!(w, "#")?;
        ["CHROM", "POS", "REF", "ALT", "VAF", "BINOM", "BETA"]
            .iter()
            .try_for_each(|x| write!(w, "{}\t", x))?;
        writeln!(w)?;
        Ok(())
    }

    pub fn write(&self, mut w: impl io::Write) -> Result<()> {
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
