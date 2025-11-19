use std::fmt;

/// The error rates for different Illumina sequencing platforms
#[derive(
    Debug,
    Copy,
    Clone,
    Default,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    clap::ValueEnum,
    serde::Serialize,
    serde::Deserialize,
)]
#[allow(clippy::doc_markdown, reason = "custom names for sequencing platforms")]
pub enum ErrorModel {
    /// MiSeq <https://support.illumina.com/sequencing/sequencing_instruments/miseq.html>
    Miseq,
    /// MiniSeq <https://support.illumina.com/sequencing/sequencing_instruments/miniseq.html>
    Miniseq,
    /// NextSeq 500 <https://support.illumina.com/sequencing/sequencing_instruments/nextseq-500.html>
    Nextseq500,
    /// NextSeq 550 <https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550.html>
    Nextseq550,
    /// HiSeq 2500 <https://support.illumina.com/sequencing/sequencing_instruments/hiseq_2500.html>
    Hiseq2500,
    /// NovaSeq 6000 <https://support.illumina.com/sequencing/sequencing_instruments/novaseq-6000.html>
    #[default] // as in rastair 1
    Novaseq6000,
    /// HiSeq X Ten <https://support.illumina.com/sequencing/sequencing_instruments/hiseq-x.html>
    HiseqXTen,
}

impl fmt::Display for ErrorModel {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ErrorModel::Miseq => write!(f, "MiSeq"),
            ErrorModel::Miniseq => write!(f, "MiniSeq"),
            ErrorModel::Nextseq500 => write!(f, "NextSeq500"),
            ErrorModel::Nextseq550 => write!(f, "NextSeq550"),
            ErrorModel::Hiseq2500 => write!(f, "HiSeq2500"),
            ErrorModel::Novaseq6000 => write!(f, "NovaSeq6000"),
            ErrorModel::HiseqXTen => write!(f, "HiSeq X Ten"),
        }
    }
}

// impl FromStr for ErrorModel {
//     type Err = String;

//     fn from_str(s: &str) -> Result<Self, Self::Err> {
//         match s.to_lowercase().as_str() {
//             "miseq" => Ok(ErrorModel::Miseq),
//             "miniseq" => Ok(ErrorModel::Miniseq),
//             "nextseq500" => Ok(ErrorModel::Nextseq500),
//             "nextseq550" => Ok(ErrorModel::Nextseq550),
//             "hiseq2500" => Ok(ErrorModel::Hiseq2500),
//             "novaseq6000" => Ok(ErrorModel::Novaseq6000),
//             "hiseq-x-ten" => Ok(ErrorModel::HiseqXTen),
//             _ => Err(format!("Unknown error model: {}", s)),
//         }
//     }
// }

impl ErrorModel {
    /// The error rate for the given error model
    ///
    /// Cf. Nicholas Stoler, Anton Nekrutenko, Sequencing error profiles of
    /// Illumina sequencing instruments, NAR Genomics and Bioinformatics, Volume
    /// 3, Issue 1, March 2021, lqab019, <https://doi.org/10.1093/nargab/lqab019>
    pub fn error_rate(&self) -> f64 {
        match self {
            ErrorModel::Miseq => 0.00473,
            ErrorModel::Miniseq => 0.00613,
            ErrorModel::Nextseq500 => 0.00429,
            ErrorModel::Nextseq550 => 0.00593,
            ErrorModel::Hiseq2500 => 0.00112,
            ErrorModel::Novaseq6000 => 0.00109,
            ErrorModel::HiseqXTen => 0.00087,
        }
    }
}
