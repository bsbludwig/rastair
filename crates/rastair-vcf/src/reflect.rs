//! Module for extraction description of VCF format and fields

use crate::{FormatFieldNumber, InfoFieldNumber};
use color_eyre::Result;
use smol_str::SmolStr;
use std::io::Write;

/// Informational structure for VCF records
#[derive(Debug)]
pub struct VcfDescription {
    /// Filters that can be applied to the VCF records
    pub filters: Vec<Filter>,
    /// Info fields containing metadata about the variants
    pub infos: Vec<Info>,
    /// Fields with calling-data for each sample
    pub formats: Vec<Format>,
}

/// Represents a filter in a VCF record
#[derive(Debug)]
pub struct Filter {
    /// The name of the filter
    pub name: SmolStr,
    /// The description of the filter
    pub description: SmolStr,
}

/// Represents an informational field in a VCF record
#[derive(Debug)]
pub struct Info {
    /// The name of the info field
    pub name: SmolStr,
    /// The description of the info field
    pub description: SmolStr,
    /// The occurrence of the info field
    pub number: InfoFieldNumber,
    /// The type of the info field
    pub field_type: SmolStr,
    /// The Rust type of the info field
    pub rust_type: SmolStr,
}

/// Represents a format field in a VCF record
#[derive(Debug)]
pub struct Format {
    /// The name of the info field
    pub name: SmolStr,
    /// The description of the info field
    pub description: SmolStr,
    /// The occurrence of the info field
    pub number: FormatFieldNumber,
    /// The type of the info field
    pub field_type: SmolStr,
    /// The Rust type of the info field
    pub rust_type: SmolStr,
}

impl VcfDescription {
    /// Write the VCF description to a markdown file
    pub fn to_markdown(&self, mut writer: impl Write) -> Result<()> {
        writeln!(writer, "# VCF Fields")?;

        writeln!(
            writer,
            "Rastair's output follows the [VCFv4.5 specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf)."
        )?;

        writeln!(writer, "## Filters")?;
        writeln!(writer, "| Name | Description |")?;
        writeln!(writer, "| -- | -- |")?;
        for filter in &self.filters {
            writeln!(writer, "| **`{}`** | {} |", filter.name, filter.description)?;
        }

        writeln!(writer, "## Info Fields")?;
        writeln!(writer, "| Name | Description | VCF Type | Rust Type | Occurance |")?;
        writeln!(writer, "| -- | -- | -- | -- | -- |")?;
        for field in &self.infos {
            let rust_type = field.rust_type.split("::").last().unwrap_or(&field.rust_type);
            writeln!(
                writer,
                "| **`{}`** | {} | `{}` | `{}` | {} |",
                field.name, field.description, field.field_type, rust_type, field.number
            )?;
        }

        writeln!(writer, "## Format Fields")?;
        writeln!(writer, "| Name | Description | VCF Type | Rust Type | Occurance |")?;
        writeln!(writer, "| -- | -- | -- | -- | -- |")?;
        for field in &self.formats {
            let rust_type = field.rust_type.split("::").last().unwrap_or(&field.rust_type);
            writeln!(
                writer,
                "| **`{}`** | {} | `{}` | `{}` | {} |",
                field.name, field.description, field.field_type, rust_type, field.number
            )?;
        }
        Ok(())
    }
}
