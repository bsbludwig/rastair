/// Write VCF records and headers.
pub trait WriteToVcf {
    /// Configuration type for controlling which fields to write
    type Config;

    /// Add headers for this type
    fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()>;

    /// Write all data to the VCF record
    ///
    /// The `config` parameter controls which fields are actually written to the VCF.
    fn write(
        &self,
        record: &mut rust_htslib::bcf::Record,
        config: &Self::Config,
    ) -> color_eyre::Result<()>;
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
        info: [$($info:tt)*],
        format: [$($format:tt)*],
        min_samples: $min_samples:expr
    ) => {
        $crate::vcf_record! {
            @with_defaults
            filters: [$($filter),*],
            info: [],
            info_defaults: [],
            info_remaining: [$($info)*],
            format: [],
            format_defaults: [],
            format_remaining: [$($format)*],
            min_samples: $min_samples
        }
    };

    // Parse info field with "default"
    (
        @with_defaults
        filters: [$($filter:ident),*],
        info: [$($info:ident),*],
        info_defaults: [$($info_default:ident),*],
        info_remaining: [$first:ident default $(, $($rest:tt)*)?],
        format: [$($format:ident),*],
        format_defaults: [$($format_default:ident),*],
        format_remaining: [$($format_rest:tt)*],
        min_samples: $min_samples:expr
    ) => {
        $crate::vcf_record! {
            @with_defaults
            filters: [$($filter),*],
            info: [$($info,)* $first],
            info_defaults: [$($info_default,)* $first],
            info_remaining: [$($($rest)*)?],
            format: [$($format),*],
            format_defaults: [$($format_default),*],
            format_remaining: [$($format_rest)*],
            min_samples: $min_samples
        }
    };

    // Parse info field without "default"
    (
        @with_defaults
        filters: [$($filter:ident),*],
        info: [$($info:ident),*],
        info_defaults: [$($info_default:ident),*],
        info_remaining: [$first:ident $(, $($rest:tt)*)?],
        format: [$($format:ident),*],
        format_defaults: [$($format_default:ident),*],
        format_remaining: [$($format_rest:tt)*],
        min_samples: $min_samples:expr
    ) => {
        $crate::vcf_record! {
            @with_defaults
            filters: [$($filter),*],
            info: [$($info,)* $first],
            info_defaults: [$($info_default),*],
            info_remaining: [$($($rest)*)?],
            format: [$($format),*],
            format_defaults: [$($format_default),*],
            format_remaining: [$($format_rest)*],
            min_samples: $min_samples
        }
    };

    // Parse format field with "default"
    (
        @with_defaults
        filters: [$($filter:ident),*],
        info: [$($info:ident),*],
        info_defaults: [$($info_default:ident),*],
        info_remaining: [],
        format: [$($format:ident),*],
        format_defaults: [$($format_default:ident),*],
        format_remaining: [$first:ident default $(, $($rest:tt)*)?],
        min_samples: $min_samples:expr
    ) => {
        $crate::vcf_record! {
            @with_defaults
            filters: [$($filter),*],
            info: [$($info),*],
            info_defaults: [$($info_default),*],
            info_remaining: [],
            format: [$($format,)* $first],
            format_defaults: [$($format_default,)* $first],
            format_remaining: [$($($rest)*)?],
            min_samples: $min_samples
        }
    };

    // Parse format field without "default"
    (
        @with_defaults
        filters: [$($filter:ident),*],
        info: [$($info:ident),*],
        info_defaults: [$($info_default:ident),*],
        info_remaining: [],
        format: [$($format:ident),*],
        format_defaults: [$($format_default:ident),*],
        format_remaining: [$first:ident $(, $($rest:tt)*)?],
        min_samples: $min_samples:expr
    ) => {
        $crate::vcf_record! {
            @with_defaults
            filters: [$($filter),*],
            info: [$($info),*],
            info_defaults: [$($info_default),*],
            info_remaining: [],
            format: [$($format,)* $first],
            format_defaults: [$($format_default),*],
            format_remaining: [$($($rest)*)?],
            min_samples: $min_samples
        }
    };

    // All fields parsed, generate the actual code
    (
        @with_defaults
        filters: [$($filter:ident),*],
        info: [$($info:ident),*],
        info_defaults: [$($info_default:ident),*],
        info_remaining: [],
        format: [$($format:ident),*],
        format_defaults: [$($format_default:ident),*],
        format_remaining: [],
        min_samples: $min_samples:expr
    ) => {pastey::paste!{
        /// Filters that can be applied to a VCF record
        #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
        pub struct Filters(rastair_types::SmallVec<rastair_types::SmolStr, 8>);

        #[allow(unused)]
        impl Filters {
            pub fn new() -> Self {
                Self(rastair_types::SmallVec::new())
            }

            pub fn add(&mut self, filter: rastair_types::SmolStr) {
                self.0.push(filter);
            }

            pub fn len(&self) -> usize {
                self.0.len()
            }

            pub fn is_empty(&self) -> bool {
                self.0.is_empty()
            }

            pub fn clear(&mut self) {
                self.0.clear();
            }

            pub fn pass(&self) -> bool {
                self.0.is_empty() || (self.0.len() == 1 && self.0[0] == "PASS")
            }
        }

        impl Default for Filters {
            fn default() -> Self {
                Self::new()
            }
        }

        impl $crate::WriteToVcf for Filters {
            type Config = FieldConfig;

            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    header.push_record($filter::header().as_bytes());
                )*
                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record, _config: &Self::Config) -> color_eyre::Result<()> {
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
        #[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
        pub struct Info {
            $(
                pub [<$info:snake>] : $info,
            )*
        }

        impl $crate::WriteToVcf for Info {
            type Config = FieldConfig;

            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    <$info as $crate::InfoField>::write_header(header)?;
                )*
                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record, config: &Self::Config) -> color_eyre::Result<()> {
                $(
                    if config.info.[<$info:snake>] {
                        $crate::InfoField::write(&self.[<$info:snake>], record)?;
                    }
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
            type Config = FieldConfig;

            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                $(
                    <$format as $crate::FormatField>::write_header(header)?;
                )*
                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record, config: &Self::Config) -> color_eyre::Result<()> {
                $(
                    if config.format.[<$format:snake>] {
                        $crate::FormatField::write(&self.[<$format:snake>], record)?;
                    }
                )*
                Ok(())
            }
        }

        /// Configuration for which INFO fields to write to VCF
        #[derive(Debug, Clone)]
        pub struct InfoFieldConfig {
            $(
                pub [<$info:snake>]: bool,
            )*
        }

        impl Default for InfoFieldConfig {
            fn default() -> Self {
                let mut config = Self {
                    $(
                        [<$info:snake>]: false,
                    )*
                };
                $(
                    config.[<$info_default:snake>] = true;
                )*
                config
            }
        }

        /// Configuration for which FORMAT fields to write to VCF
        #[derive(Debug, Clone)]
        pub struct FormatFieldConfig {
            $(
                pub [<$format:snake>]: bool,
            )*
        }

        impl Default for FormatFieldConfig {
            fn default() -> Self {
                let mut config = Self {
                    $(
                        [<$format:snake>]: false,
                    )*
                };
                $(
                    config.[<$format_default:snake>] = true;
                )*
                config
            }
        }

        /// Configuration for which fields to write to VCF
        #[derive(Debug, Clone, Default)]
        pub struct FieldConfig {
            pub info: InfoFieldConfig,
            pub format: FormatFieldConfig,
        }

        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct InfoFieldId(pub &'static cstr8::CStr8);

        #[allow(unused)]
        impl InfoFieldId {
            pub const ALL_IDS: &'static [&'static str] = &[
                $($info::ID.as_str()),*
            ];
        }

        impl std::str::FromStr for InfoFieldId {
            type Err = String;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                $(
                    if s == $info::ID.as_str() {
                        return Ok(InfoFieldId($info::ID));
                    }
                )*
                let available = vec![$($info::ID.as_str()),*].join(", ");
                Err(format!("Unknown INFO field: '{}'. Available: {}", s, available))
            }
        }

        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct FormatFieldId(pub &'static cstr8::CStr8);

        #[allow(unused)]
        impl FormatFieldId {
            pub const ALL_IDS: &'static [&'static str] = &[
                $($format::ID.as_str()),*
            ];
        }

        impl std::str::FromStr for FormatFieldId {
            type Err = String;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                $(
                    if s == $format::ID.as_str() {
                        return Ok(FormatFieldId($format::ID));
                    }
                )*
                let available = vec![$($format::ID.as_str()),*].join(", ");
                Err(format!("Unknown FORMAT field: '{}'. Available: {}", s, available))
            }
        }

        impl FieldConfig {
            pub fn with_field_ids(
                mut self,
                info_fields: &[InfoFieldId],
                format_fields: &[FormatFieldId],
            ) -> Self {
                for id in info_fields {
                    $(if id.0.as_str() == $info::ID.as_str() { self.info.[<$info:snake>] = true; })*
                }
                for id in format_fields {
                    $(if id.0.as_str() == $format::ID.as_str() { self.format.[<$format:snake>] = true; })*
                }
                self
            }

            #[allow(unused)]
            pub fn with_all_fields(mut self) -> Self {
                $( self.info.[<$info:snake>] = true; )*
                $( self.format.[<$format:snake>] = true; )*
                self
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
            pub samples: rastair_types::SmallVec<Format, $min_samples>,
        }

        impl $crate::WriteToVcf for Record {
            type Config = FieldConfig;

            fn write_header(header: &mut rust_htslib::bcf::Header) -> color_eyre::Result<()> {
                use color_eyre::eyre::WrapErr;

                Filters::write_header(header).wrap_err("Failed to write filters header")?;
                Info::write_header(header).wrap_err("Failed to write info header")?;
                Format::write_header(header).wrap_err("Failed to write format header")?;

                Ok(())
            }

            fn write(&self, record: &mut rust_htslib::bcf::Record, config: &Self::Config) -> color_eyre::Result<()> {
                use color_eyre::eyre::WrapErr;

                self.main.write(record)?;
                self.filters.write(record, config)
                    .wrap_err("Failed to write filters")?;
                self.info.write(record, config)
                    .wrap_err("Failed to write info")?;
                self.samples.iter().try_for_each(|sample| {
                    sample.write(record, config)
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

#[cfg(test)]
mod tests {
    // Test module for developing the new field config syntax
    mod field_config_test {
        use crate::*;

        filter!(PASS, "All filters pass");
        filter!(LowQual, "Low quality");

        info_field!(TestDepth(usize), "TD", "Test depth", InfoFieldNumber::Num(1));
        info_field!(TestQuality(f64), "TQ", "Test quality", InfoFieldNumber::Num(1));
        info_field!(TestFlag(bool), "TF", "Test flag", InfoFieldNumber::Num(0));

        format_field!(TestGenotype(usize), "TGT", "Test genotype", FormatFieldNumber::Num(1));
        format_field!(TestScore(f64), "TS", "Test score", FormatFieldNumber::Num(1));

        // Test with minimal fields - will expand macro to generate FieldConfig
        // TestDepth and TestGenotype are marked as default
        vcf_record!(
            filters: [PASS, LowQual],
            info: [TestDepth default, TestQuality, TestFlag],
            format: [TestGenotype default, TestScore],
            min_samples: 1
        );

        #[test]
        fn test_basic_types_exist() {
            // Smoke test - just verify the macro generates the expected types
            let filters = Filters::new();
            assert!(filters.is_empty());

            let info = Info::default();
            let _ = info;

            let format = Format::default();
            let _ = format;
        }

        #[test]
        fn test_field_config_generation() {
            // Test that FieldConfig structs are generated with correct defaults
            let config = FieldConfig::default();

            // TestDepth and TestGenotype should default to true (marked with "default")
            assert!(config.info.test_depth, "TestDepth should be in default set");
            assert!(!config.info.test_quality, "TestQuality should not be in default set");
            assert!(!config.info.test_flag, "TestFlag should not be in default set");

            assert!(config.format.test_genotype, "TestGenotype should be in default set");
            assert!(!config.format.test_score, "TestScore should not be in default set");

            // Test that we can modify the config
            let mut config = config;
            config.info.test_quality = true;
            assert!(config.info.test_quality);
        }

        #[test]
        fn test_with_field_ids() {
            let info_ids: Vec<InfoFieldId> =
                ["TQ", "TF"].iter().map(|s| s.parse().unwrap()).collect();
            let format_ids: Vec<FormatFieldId> =
                ["TS"].iter().map(|s| s.parse().unwrap()).collect();
            let config = FieldConfig::default().with_field_ids(&info_ids, &format_ids);

            // Default fields should still be enabled
            assert!(config.info.test_depth, "TestDepth should still be default");
            assert!(config.format.test_genotype, "TestGenotype should still be default");

            // Additional fields should now be enabled
            assert!(config.info.test_quality, "TestQuality should be enabled");
            assert!(config.info.test_flag, "TestFlag should be enabled");
            assert!(config.format.test_score, "TestScore should be enabled");
        }

        #[test]
        fn test_field_id_parsing_invalid() {
            let result = "INVALID".parse::<InfoFieldId>();
            assert!(result.is_err(), "Should error on invalid INFO field");

            let result = "INVALID".parse::<FormatFieldId>();
            assert!(result.is_err(), "Should error on invalid FORMAT field");
        }

        #[test]
        fn test_with_all_fields() {
            let config = FieldConfig::default().with_all_fields();

            // Default fields should still be enabled
            assert!(config.info.test_depth, "TestDepth should still be default");
            assert!(config.format.test_genotype, "TestGenotype should still be default");

            // Additional fields should now be enabled
            assert!(config.info.test_quality, "TestQuality should be enabled");
            assert!(config.info.test_flag, "TestFlag should be enabled");
            assert!(config.format.test_score, "TestScore should be enabled");
        }
    }
}
