use std::fmt::{self};
use tracing::instrument;

impl super::Rastair1BedFormat {
    #[instrument(level = "debug", skip(self), fields(contig=%self.contig, pos=self.pos))]
    pub fn sanity_check(&self) -> Option<Errors> {
        let mut errors = vec![];

        if self.beta.is_none() {
            errors.push("Missing beta value".to_string());
        }

        if let Some(beta) = self.beta
            && *beta > 0.
            && self.r#mod == 0
        {
            errors.push(format!(
                "Inconsistent beta and mod counts: beta={}, mod={}",
                beta, self.r#mod
            ));
        }

        if self.coverage < (self.unmod + self.r#mod + self.no_snp + self.snp) as usize {
            errors.push(format!(
                "Coverage mismatch: coverage={}, unmod+mod+no_snp+snp={}",
                self.coverage,
                self.unmod + self.r#mod + self.no_snp + self.snp
            ));
        }

        if errors.is_empty() { None } else { Some(Errors(errors)) }
    }
}

pub struct Errors(Vec<String>);

impl fmt::Display for Errors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for error in &self.0 {
            writeln!(f, "- {}", error)?;
        }
        Ok(())
    }
}
