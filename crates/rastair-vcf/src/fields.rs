mod format_field;
mod impls_for_wrappers;
mod info_field;

pub use format_field::{FormatField, FormatFieldNumber, FormatFieldValue};
pub use info_field::{InfoField, InfoFieldNumber, InfoFieldValue, StrandSpecificInfoField};

/// A field that can be used in the header section.
pub trait HeaderField {
    /// Description of the field, used in the VCF header.
    const DESCRIPTION: &'static str;
}

/// A field that can be used in VCF.
pub trait VcfField: HeaderField {
    /// The ID of the field, used in the VCF header and record.
    const ID: &'static cstr8::CStr8;
}
