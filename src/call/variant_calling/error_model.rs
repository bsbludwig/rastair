use seqair_types::Probability;
use std::fmt;
use std::str::FromStr;

macro_rules! define_error_models {
    (
        $(
            $variant:ident($cli:expr, $display:expr, $rate:expr, $url:expr)
        ),+ $(,)?
    ) => {
        /// The error rates for different Illumina sequencing platforms
        #[derive(Debug, Copy, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
        #[allow(clippy::doc_markdown, reason = "custom names for sequencing platforms")]
        pub enum ErrorModel {
            $(
                #[doc = concat!($display, " <", $url, ">")]
                $variant,
            )+
            /// Custom error rate
            Custom(Probability),
        }

        impl fmt::Display for ErrorModel {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                match self {
                    $(
                        ErrorModel::$variant => write!(f, $display),
                    )+
                    ErrorModel::Custom(p) => write!(f, "{}", **p),
                }
            }
        }

        impl FromStr for ErrorModel {
            type Err = color_eyre::Report;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                let lowercase = s.to_lowercase();
                match lowercase.as_str() {
                    $(
                        $cli => Ok(ErrorModel::$variant),
                    )+
                    _ => {
                        use color_eyre::eyre::WrapErr;

                        let value: f64 = s.parse().map_err(|_| {
                            color_eyre::eyre::eyre!(
                                "Invalid error model: '{}'. Expected a platform name ({}) or a numeric error rate",
                                s,
                                concat!($($cli, ", "),+).trim_end_matches(", ")
                            )
                        })?;

                        let probability = Probability::new(value).wrap_err("Error rate must be between 0.0 and 1.0")?;

                        Ok(ErrorModel::Custom(probability))
                    }
                }
            }
        }



        impl ErrorModel {
            /// The error rate for the given error model
            ///
            /// Cf. Nicholas Stoler, Anton Nekrutenko, Sequencing error profiles of
            /// Illumina sequencing instruments, NAR Genomics and Bioinformatics, Volume
            /// 3, Issue 1, March 2021, lqab019, <https://doi.org/10.1093/nargab/lqab019>
            pub fn error_rate(&self) -> Probability {
                match self {
                    $(
                        ErrorModel::$variant => Probability::new_panicky($rate),
                    )+
                    ErrorModel::Custom(p) => *p,
                }
            }

            /// Creates a custom value parser for clap that shows possible values
            /// but also accepts custom numeric error rates
            pub fn value_parser() -> ErrorModelValueParser {
                ErrorModelValueParser
            }
        }

        #[derive(Clone)]
        pub struct ErrorModelValueParser;

        impl clap::builder::TypedValueParser for ErrorModelValueParser {
            type Value = ErrorModel;

            fn parse_ref(
                &self,
                _cmd: &clap::Command,
                _arg: Option<&clap::Arg>,
                value: &std::ffi::OsStr,
            ) -> Result<Self::Value, clap::Error> {
                let s = value.to_str().ok_or_else(|| {
                    clap::Error::raw(
                        clap::error::ErrorKind::InvalidUtf8,
                        "Invalid UTF-8"
                    )
                })?;

                s.parse::<ErrorModel>().map_err(|e| {
                    clap::Error::raw(
                        clap::error::ErrorKind::InvalidValue,
                        format!("{}", e)
                    )
                })
            }

            fn possible_values(&self) -> Option<Box<dyn Iterator<Item = clap::builder::PossibleValue> + '_>> {
                Some(Box::new([
                    $(
                        clap::builder::PossibleValue::new($cli).help(format!("{} <{}>", $display, $url)),
                    )+
                ].into_iter()))
            }
        }
    };
}

define_error_models! {
    Miseq("miseq", "MiSeq", 0.00473, "https://support.illumina.com/sequencing/sequencing_instruments/miseq.html"),
    Miniseq("miniseq", "MiniSeq", 0.00613, "https://support.illumina.com/sequencing/sequencing_instruments/miniseq.html"),
    Nextseq500("nextseq500", "NextSeq500", 0.00429, "https://support.illumina.com/sequencing/sequencing_instruments/nextseq-500.html"),
    Nextseq550("nextseq550", "NextSeq550", 0.00593, "https://support.illumina.com/sequencing/sequencing_instruments/nextseq-550.html"),
    Hiseq2500("hiseq2500", "HiSeq2500", 0.00112, "https://support.illumina.com/sequencing/sequencing_instruments/hiseq_2500.html"),
    Novaseq6000("novaseq6000", "NovaSeq6000", 0.00109, "https://support.illumina.com/sequencing/sequencing_instruments/novaseq-6000.html"),
    Hiseqxten("hiseqxten", "HiSeq X Ten", 0.00087, "https://support.illumina.com/sequencing/sequencing_instruments/hiseq-x.html"),
}

#[expect(clippy::derivable_impls, reason = "generated by macro")]
impl Default for ErrorModel {
    fn default() -> Self {
        ErrorModel::Novaseq6000
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_platform_names() {
        assert!(matches!("miseq".parse::<ErrorModel>().unwrap(), ErrorModel::Miseq));
        assert!(matches!("MiSeq".parse::<ErrorModel>().unwrap(), ErrorModel::Miseq));
        assert!(matches!("MISEQ".parse::<ErrorModel>().unwrap(), ErrorModel::Miseq));
        assert!(matches!("novaseq6000".parse::<ErrorModel>().unwrap(), ErrorModel::Novaseq6000));
    }

    #[test]
    fn test_parse_custom_error_rate() {
        let model = "0.005".parse::<ErrorModel>().unwrap();
        assert!(matches!(model, ErrorModel::Custom(_)));
        assert_eq!(*model.error_rate(), 0.005);
    }

    #[test]
    fn test_parse_custom_error_rate_boundaries() {
        let model_zero = "0.0".parse::<ErrorModel>().unwrap();
        assert!(matches!(model_zero, ErrorModel::Custom(_)));
        assert_eq!(*model_zero.error_rate(), 0.0);

        let model_one = "1.0".parse::<ErrorModel>().unwrap();
        assert!(matches!(model_one, ErrorModel::Custom(_)));
        assert_eq!(*model_one.error_rate(), 1.0);
    }

    #[test]
    fn test_parse_invalid_error_rate_out_of_range() {
        assert!("1.5".parse::<ErrorModel>().is_err());
        assert!("-0.1".parse::<ErrorModel>().is_err());
    }

    #[test]
    fn test_parse_invalid_input() {
        assert!("invalid_platform".parse::<ErrorModel>().is_err());
        assert!("not_a_number".parse::<ErrorModel>().is_err());
    }

    #[test]
    fn test_error_rate_for_platforms() {
        assert_eq!(*ErrorModel::Miseq.error_rate(), 0.00473);
        assert_eq!(*ErrorModel::Novaseq6000.error_rate(), 0.00109);
    }

    #[test]
    fn test_error_rate_for_custom() {
        let model = ErrorModel::Custom(Probability::new(0.01).unwrap());
        assert_eq!(*model.error_rate(), 0.01);
    }

    #[test]
    fn test_display_platforms() {
        assert_eq!(ErrorModel::Miseq.to_string(), "MiSeq");
        assert_eq!(ErrorModel::Novaseq6000.to_string(), "NovaSeq6000");
    }

    #[test]
    fn test_display_custom() {
        let model = ErrorModel::Custom(Probability::new(0.005).unwrap());
        assert_eq!(model.to_string(), "0.005");
    }

    #[test]
    fn test_default() {
        assert!(matches!(ErrorModel::default(), ErrorModel::Novaseq6000));
    }
}
