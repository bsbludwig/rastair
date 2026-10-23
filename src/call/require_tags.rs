use crate::utils::cli;
use seqair::bam::aux::Aux;
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
    /// Format: TAG=VALUE, e.g. `--require-tags RG=mygroup`. Accepts one or more values
    /// (space-separated). A read is kept only if it matches all of the specified
    /// tag=value pairs.
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
                    .all(|f| aux.get::<&[u8]>(f.tag).is_ok_and(|value| value == f.value.as_bytes()))
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

    #[test]
    fn a_tag_of_another_type_does_not_match() {
        let mut record = AuxData::new();
        record.set_int(*b"XX", 7).unwrap();
        assert!(!allows(&filter(&["XX=7"]), &record));
    }
}
