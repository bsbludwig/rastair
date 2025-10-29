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
        #[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
        pub struct Filters($crate::FilterContainer);

        impl Filters {
            /// Create a new, empty set of filters
            pub fn new() -> Self {
                Self($crate::FilterContainer::new())
            }
        }

        impl std::ops::Deref for Filters {
            type Target = $crate::FilterContainer;

            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl std::ops::DerefMut for Filters {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
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
                self.0.write_to_record(record)
            }
        }

        /// Info fields for a VCF record containing various metadata
        #[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
        pub struct Info {
            $(
                pub [<$info:snake>] : $info,
            )*
        }

        impl $crate::WriteToVcf for Info {
            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    <$info as $crate::InfoField>::write_header(header)?;
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
        #[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
        pub struct Format {
            $(
                pub [<$format:snake>]: $format,
            )*
        }

        impl $crate::WriteToVcf for Format {
            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    <$format as $crate::FormatField>::write_header(header)?;
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
        #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
        #[allow(non_camel_case_types, non_snake_case, clippy::upper_case_acronyms)]
        pub struct Record {
            /// Fixed fields for the VCF record, such as chromosome, position, ID, reference, and alternate alleles
            pub main:  VcfFixedFields,
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

                self.main.write(record)?;
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

        #[allow(unused)]
        impl Record {
            /// Get description of the VCF record
            pub fn description() -> $crate::reflect::VcfDescription {
                $crate::reflect::VcfDescription {
                    filters: vec![$( $filter::description() ),*],
                    infos: {
                        let mut res = vec![];
                        $( res.extend($info::description()); )*
                        res
                    },
                    formats: {
                        let mut res = vec![];
                        $( res.extend($format::description()); )*
                        res
                    },
                }
            }
        }
    }
}}
