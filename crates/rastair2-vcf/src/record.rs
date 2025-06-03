/// Write VCF records and headers.
pub trait WriteToVcf {
    /// Add headers for this type
    fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()>;

    /// Write all data to the VCF record
    fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()>;
}

/// Macro to define a VCF record with filters, info, format, and samples
#[macro_export]
macro_rules! vcf_record {
    (
        filters: [$($filter:ident),* $(,)?],
        info: [$($info:ident),* $(,)?],
        format: [$($format:ident),* $(,)?],
        min_samples: $min_samples:expr
    ) => {
        #[allow(non_camel_case_types, non_snake_case, clippy::upper_case_acronyms)]
        pub struct Filters(smallvec::SmallVec<&'static str, 5>);

        impl Filters {
            pub fn new() -> Self {
                Self(smallvec::SmallVec::new())
            }

            pub fn add(mut self, filter: impl $crate::VcfFilter) -> Self {
                self.0.push(filter.filter());
                self
            }
        }

        impl $crate::WriteToVcf for Filters {
            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    header.push_record($filter::header().as_bytes());
                )*
                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                use color_eyre::eyre::WrapErr;

                self.0.iter().try_for_each(|filter| {
                    record.push_filter(filter.as_bytes()).wrap_err_with(|| {
                        format!("Failed to push filter {filter}")
                    })
                })?;

                Ok(())
            }
        }

        #[allow(non_camel_case_types, non_snake_case, clippy::upper_case_acronyms)]
        pub struct Info {
            $(
                pub $info: $info,
            )*
        }

        impl $crate::WriteToVcf for Info {
            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    header.push_record($info::header().as_bytes());
                )*
                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                $(
                    self.$info.write(record)?;
                )*
                Ok(())
            }
        }

        #[allow(non_camel_case_types, non_snake_case, clippy::upper_case_acronyms)]
        pub struct Format {
            $(
                pub $format: $format,
            )*
        }

        impl $crate::WriteToVcf for Format {
            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    header.push_record($format::header().as_bytes());
                )*
                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                $(
                    self.$format.write(record)?;
                )*
                Ok(())
            }
        }

        #[allow(non_camel_case_types, non_snake_case, clippy::upper_case_acronyms)]
        pub struct Record {
            pub fixed_fields: VcfFixedFields,
            pub filters: Filters,
            pub info: Info,
            pub samples: smallvec::SmallVec<Format, $min_samples>,
        }

        impl $crate::WriteToVcf for Record {
            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                use color_eyre::eyre::WrapErr;

                Filters::write_header(header).wrap_err("Failed to write filters header")?;
                Info::write_header(header).wrap_err("Failed to write info header")?;
                Format::write_header(header).wrap_err("Failed to write format header")?;

                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                use color_eyre::eyre::WrapErr;

                self.fixed_fields.write(record)?;
                self.filters.write(record)
                    .wrap_err("Failed to write filters")?;
                self.info.write(record)
                    .wrap_err("Failed to write info")?;
                self.samples.iter().try_for_each(|sample| {
                    sample.write(record)
                        .wrap_err("Failed to write format")
                }).wrap_err("Failed to write samples")?;

                Ok(())
            }
        }
    };
}
