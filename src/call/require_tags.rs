use crate::utils::cli;
use seqair::bam::aux::{Aux, AuxValue, find_tag};
use seqair_types::SmolStr;
use std::{str::FromStr, sync::Arc};

/// A parsed `TAG=VALUE` filter argument.
#[derive(Debug, Clone)]
pub struct TagValue {
    tag: [u8; 2],
    value: SmolStr,
}

impl FromStr for TagValue {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let Some((tag_str, value)) = s.split_once('=') else {
            return Err(format!("expected TAG=VALUE, got `{s}`"));
        };
        if tag_str.len() != 2 || !tag_str.is_ascii() {
            return Err(format!("tag must be exactly 2 ASCII characters, got `{tag_str}`"));
        }
        let b = tag_str.as_bytes();
        Ok(TagValue { tag: [b[0], b[1]], value: SmolStr::new(value) })
    }
}

impl TagValue {
    /// `VALUE` is read as the type the tag has in the read: `Z` and `H` match
    /// byte for byte, `A` as its one character, the integer types (`cCsSiI`)
    /// by value, so `NM=0` and `NM=+0` match `NM:i:0`, and `f`/`d` by value
    /// at their own precision. A `VALUE` that does not parse as that type, and
    /// a `B` array, never match.
    fn matches(&self, value: &AuxValue<'_>) -> bool {
        let text = self.value.as_str();
        match *value {
            AuxValue::String(bytes) | AuxValue::Hex(bytes) => bytes == text.as_bytes(),
            AuxValue::Char(c) => text.as_bytes() == [c],
            AuxValue::I8(_)
            | AuxValue::U8(_)
            | AuxValue::I16(_)
            | AuxValue::U16(_)
            | AuxValue::I32(_)
            | AuxValue::U32(_) => text.parse::<i64>().ok() == value.as_i64(),
            AuxValue::Float(v) => text.parse::<f32>().ok() == Some(v),
            AuxValue::Double(v) => text.parse::<f64>().ok() == Some(v),
            AuxValue::ArrayI8(_)
            | AuxValue::ArrayU8(_)
            | AuxValue::ArrayI16(_)
            | AuxValue::ArrayU16(_)
            | AuxValue::ArrayI32(_)
            | AuxValue::ArrayU32(_)
            | AuxValue::ArrayFloat(_) => false,
        }
    }
}

impl serde::Serialize for TagValue {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let tag_str = std::str::from_utf8(&self.tag).unwrap_or("??");
        s.serialize_str(&format!("{tag_str}={}", self.value))
    }
}

#[derive(Debug, clap::Args, Clone, Default, serde::Serialize)]
pub struct RequireTagsParams {
    /// Require reads to have a specific SAM tag value
    ///
    /// Format: TAG=VALUE, e.g. `--require-tags RG=mygroup HP=1`. Accepts one or more values
    /// (space-separated). A read is kept only if it matches all of the specified
    /// tag=value pairs. VALUE is compared as the tag's SAM type: strings exactly,
    /// integers and floats by value (`NM=0` matches `NM:i:0`); array tags never match.
    #[arg(long, num_args(1..))]
    #[arg(help_heading = cli::sections::INPUT)]
    pub require_tags: Vec<TagValue>,
}

/// A pre-processed, cheaply cloneable form of the tag filter.
#[derive(Debug, Clone, Default)]
pub enum TagRequirement {
    /// No filter specified, so all reads pass.
    #[default]
    All,
    /// Only reads matching every one of these tag=value pairs pass.
    AllOf(Arc<[TagValue]>),
}

impl RequireTagsParams {
    pub fn filter(&self) -> TagRequirement {
        if self.require_tags.is_empty() {
            TagRequirement::All
        } else {
            TagRequirement::AllOf(self.require_tags.as_slice().into())
        }
    }
}

impl TagRequirement {
    /// Whether a read with these aux fields passes. `None` is an aux block that
    /// could not be read, which fails a requirement (closed, not open).
    pub fn allows(&self, aux: Option<Aux<'_>>) -> bool {
        match self {
            TagRequirement::All => true,
            TagRequirement::AllOf(filters) => aux.is_some_and(|aux| {
                filters
                    .iter()
                    .all(|f| find_tag(aux.as_bytes(), f.tag).is_some_and(|value| f.matches(&value)))
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use seqair::bam::aux_data::AuxData;

    const L001: &str = "mTet1-PyBr-16h-p1_S1_L001";
    const L002: &str = "mTet1-PyBr-16h-p1_S1_L002";

    fn filter(tag_values: &[&str]) -> TagRequirement {
        RequireTagsParams { require_tags: tag_values.iter().map(|s| s.parse().unwrap()).collect() }
            .filter()
    }

    fn aux(tags: &[(&[u8; 2], &str)]) -> AuxData {
        let mut aux = AuxData::new();
        aux.set_int(*b"NM", 3).unwrap();
        for (tag, value) in tags {
            aux.set_string(**tag, value.as_bytes());
        }
        aux
    }

    fn allows(f: &TagRequirement, aux: &AuxData) -> bool {
        f.allows(Some(Aux::new(aux.as_bytes())))
    }

    #[test]
    fn all_filter_allows_every_record() {
        let f = filter(&[]);
        assert!(matches!(f, TagRequirement::All));
        assert!(allows(&f, &aux(&[(b"RG", L001)])));
        assert!(allows(&f, &aux(&[])));
        assert!(f.allows(None), "no filter asks nothing of the aux block");
    }

    #[test]
    fn single_tag_allows_only_the_matching_value() {
        let f = filter(&[&format!("RG={L001}")]);
        assert!(allows(&f, &aux(&[(b"RG", L001)])));
        assert!(!allows(&f, &aux(&[(b"RG", L002)])));
    }

    #[test]
    fn multiple_tags_require_all_members() {
        let record = aux(&[(b"RG", L001), (b"XX", "foo"), (b"YY", "bar")]);
        assert!(allows(&filter(&[&format!("RG={L001}"), "XX=foo", "YY=bar"]), &record));
        // One of several different tags mismatches → rejected.
        assert!(!allows(&filter(&[&format!("RG={L001}"), "XX=foo", "YY=wrong"]), &record));
        // A required tag missing from the record → rejected.
        assert!(!allows(&filter(&[&format!("RG={L001}"), "XX=foo", "ZZ=anything"]), &record));
    }

    #[test]
    fn same_tag_with_multiple_values_is_unsatisfiable() {
        // A record's RG can equal only one value, so requiring two different RG
        // values at once never matches.
        let f = filter(&[&format!("RG={L001}"), &format!("RG={L002}")]);
        assert!(!allows(&f, &aux(&[(b"RG", L001)])));
        assert!(!allows(&f, &aux(&[(b"RG", L002)])));
    }

    #[test]
    fn a_record_without_the_tag_or_its_aux_block_is_rejected() {
        let f = filter(&[&format!("RG={L001}")]);
        assert!(!allows(&f, &aux(&[])));
        assert!(!f.allows(None));
    }

    fn allows_one(requirement: &str, set: impl FnOnce(&mut AuxData)) -> bool {
        let mut record = AuxData::new();
        set(&mut record);
        allows(&filter(&[requirement]), &record)
    }

    #[test]
    fn integer_tags_match_by_value_at_every_width() {
        for value in [0, 7, -1, 300, -300, 70_000, -70_000, i64::from(u32::MAX)] {
            let set = |r: &mut AuxData| r.set_int(*b"NM", value).unwrap();
            assert!(allows_one(&format!("NM={value}"), set), "NM={value}");
            assert!(!allows_one(&format!("NM={}", value + 1), set), "NM={value} vs +1");
        }
        assert!(allows_one("HP=+1", |r| r.set_int(*b"HP", 1).unwrap()));
        assert!(!allows_one("HP=one", |r| r.set_int(*b"HP", 1).unwrap()));
        assert!(!allows_one("HP=1.0", |r| r.set_int(*b"HP", 1).unwrap()));
    }

    #[test]
    fn float_tags_match_by_value_of_their_width() {
        assert!(allows_one("XF=0.1", |r| r.set_float(*b"XF", 0.1)));
        assert!(allows_one("XF=1", |r| r.set_float(*b"XF", 1.0)));
        assert!(!allows_one("XF=0.2", |r| r.set_float(*b"XF", 0.1)));
        assert!(allows_one("XD=0.1", |r| r.set_double(*b"XD", 0.1)));
        assert!(!allows_one("XD=0.1000001", |r| r.set_double(*b"XD", 0.1)));
    }

    #[test]
    fn char_tags_match_the_one_character() {
        assert!(allows_one("XA=+", |r| r.set_char(*b"XA", b'+').unwrap()));
        assert!(!allows_one("XA=-", |r| r.set_char(*b"XA", b'+').unwrap()));
        assert!(!allows_one("XA=++", |r| r.set_char(*b"XA", b'+').unwrap()));
    }

    #[test]
    fn string_and_hex_tags_match_byte_for_byte() {
        assert!(allows_one("XZ=7", |r| r.set_string(*b"XZ", b"7")));
        assert!(!allows_one("XZ=07", |r| r.set_string(*b"XZ", b"7")));
        assert!(allows_one("XH=1AE3", |r| r.set_hex(*b"XH", b"1AE3")));
        assert!(!allows_one("XH=1ae3", |r| r.set_hex(*b"XH", b"1AE3")));
    }

    #[test]
    fn array_tags_never_match() {
        assert!(!allows_one("XB=7", |r| r.set_array_u8(*b"XB", &[7]).unwrap()));
    }
}
