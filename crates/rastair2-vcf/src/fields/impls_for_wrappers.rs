use crate::{FormatFieldValue, InfoFieldValue};
use color_eyre::eyre::{Context as _, Result};
use rust_htslib::bcf::Record;
use smallvec::SmallVec;

impl InfoFieldValue for rastair2_types::rms::RootMeanSquare {
    const TYPE_NAME: &'static str = "Float";

    #[allow(clippy::cast_possible_truncation)] // Allow casting f64 to f32, which is common in VCF
    fn write(record: &mut Record, tag: &cstr8::CStr8, values: &[Self]) -> Result<()> {
        record
            .push_info_float(tag, &values.iter().map(|&n| *n as f32).collect::<SmallVec<f32, 5>>())
            .wrap_err("Failed to set field (RMS values)")
    }
}

impl FormatFieldValue for rastair2_types::rms::RootMeanSquare {
    const TYPE_NAME: &'static str = "Float";

    #[allow(clippy::cast_possible_truncation)] // Allow casting f64 to f32, which is common in VCF
    fn write(record: &mut Record, tag: &cstr8::CStr8, values: &[Self]) -> Result<()> {
        record
            .push_format_float(
                tag,
                &values.iter().map(|&n| *n as f32).collect::<SmallVec<f32, 5>>(),
            )
            .wrap_err("Failed to set field (RMS values)")
    }
}

impl InfoFieldValue for rastair2_types::phred::Phred {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &cstr8::CStr8, values: &[Self]) -> Result<()> {
        record
            .push_info_integer(
                tag,
                &values.iter().map(|&n| n.as_int()).collect::<SmallVec<i32, 5>>(),
            )
            .wrap_err("Failed to set field (Phred values)")
    }
}

impl FormatFieldValue for rastair2_types::phred::Phred {
    const TYPE_NAME: &'static str = "Integer";

    fn write(record: &mut Record, tag: &cstr8::CStr8, values: &[Self]) -> Result<()> {
        record
            .push_format_integer(
                tag,
                &values.iter().map(|&n| n.as_int()).collect::<SmallVec<i32, 5>>(),
            )
            .wrap_err("Failed to set field (Phred values)")
    }
}
