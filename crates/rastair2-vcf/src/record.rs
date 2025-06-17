/// Write VCF records and headers.
pub trait WriteToVcf {
    /// Add headers for this type
    fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()>;

    /// Write all data to the VCF record
    fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()>;
}

/// Macro to define a VCF `Record` struct alongside `Filters`, `Info` and `Format` structs.
///
/// Best call this in a module so you can refer to the types without this macro
/// polluting the namespace.
///
/// # Inline samples
///
/// Similar to how info and format fields are built by [`crate::info_field!`]
/// and [`crate::format_field!`], the `Record` struct created by this macro uses
/// a [`smallvec::SmallVec`] for the `samples` field to inline up to
/// `$min_samples` samples directly in the struct.
#[macro_export]
macro_rules! vcf_record {
    (
        filters: [$($filter:ident),* $(,)?],
        info: [$($info:ident),* $(,)?],
        format: [$($format:ident),* $(,)?],
        min_samples: $min_samples:expr
    ) => {pastey::paste!{
        /// Filters that can be applied to a VCF record
        #[derive(Debug, Clone)]
        pub struct Filters(smallvec::SmallVec<&'static str, 2>);

        impl Filters {
            pub fn new() -> Self {
                Self(smallvec::SmallVec::new())
            }

            pub fn add(&mut self, filter: impl $crate::VcfFilter) {
                self.0.push(filter.filter());
            }

            pub fn len(&self) -> usize {
                self.0.len()
            }

            pub fn is_empty(&self) -> bool {
                self.0.is_empty()
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

                if self.0.is_empty() {
                    record.set_filters::<[u8]>(&[]).wrap_err("Failed to clear filters")?;
                } else {
                    self.0.iter().try_for_each(|filter| {
                        record.push_filter(filter.as_bytes()).wrap_err_with(|| {
                            format!("Failed to push filter {filter}")
                        })
                    })?;
                }

                Ok(())
            }
        }

        /// Info fields for a VCF record containing various metadata
        #[derive(Debug, Clone)]
        pub struct Info {
            $(
                pub [<$info:snake>] : $info,
            )*
        }

        impl $crate::WriteToVcf for Info {
            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    header.push_record(<$info as $crate::InfoField>::header().as_bytes());
                )*
                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                $(
                    $crate::InfoField::write(&self.[<$info:snake>], record)?;
                )*
                Ok(())
            }
        }

        /// Format fields for a VCF record containing sample-specific data
        ///
        /// Used to add data that was "called"
        #[derive(Debug, Clone)]
        pub struct Format {
            $(
                pub [<$format:snake>]: $format,
            )*
        }

        impl $crate::WriteToVcf for Format {
            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    header.push_record(<$format as $crate::FormatField>::header().as_bytes());
                )*
                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record) -> color_eyre::Result<()> {
                $(
                    $crate::FormatField::write(&self.[<$format:snake>], record)?;
                )*
                Ok(())
            }
        }

        /// A VCF record containing fixed fields, filters, info, and format data
        #[derive(Debug, Clone)]
        #[allow(non_camel_case_types, non_snake_case, clippy::upper_case_acronyms)]
        pub struct Record {
            /// Fixed fields for the VCF record, such as chromosome, position, ID, reference, and alternate alleles
            pub fixed_fields: VcfFixedFields,
            /// Filters applied to the VCF record
            pub filters: Filters,
            /// Metrics and data about the variant
            pub info: Info,
            /// Sample-specific data for the VCF record
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
    }
}}
