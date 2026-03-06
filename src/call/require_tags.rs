use crate::utils::cli;
use rastair_types::SmolStr;
use rust_htslib::bam;
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
    /// (space-separated). A read is kept if it matches any of the specified tag=value pairs.
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
    /// Only reads matching at least one of these tag=value pairs pass.
    AnyOf(Arc<[TagValue]>),
}

impl RequireTagsParams {
    pub fn filter(&self) -> TagRequirement {
        if self.require_tags.is_empty() {
            TagRequirement::All
        } else {
            TagRequirement::AnyOf(self.require_tags.as_slice().into())
        }
    }
}

impl TagRequirement {
    /// Returns `true` if a read passes the tag filter.
    pub fn allows(&self, record: &bam::RecordView<'_>) -> bool {
        match self {
            TagRequirement::All => true,
            TagRequirement::AnyOf(filters) => {
                let Some(aux) = record.raw_aux_data() else {
                    return false;
                };
                filters
                    .iter()
                    .any(|f| find_tag(aux, f.tag).is_some_and(|v| v == f.value.as_bytes()))
            }
        }
    }
}

/// Scan raw BAM auxiliary data for a tag and return its string (`Z`-typed) value as raw bytes.
fn find_tag(data: &[u8], target: [u8; 2]) -> Option<&[u8]> {
    let mut pos = 0;
    #[allow(clippy::indexing_slicing, reason = "loop condition checks bounds")]
    while pos + 3 <= data.len() {
        let tag0 = data[pos];
        let tag1 = data[pos + 1];
        let typ = data[pos + 2];
        pos += 3;

        if tag0 == target[0] && tag1 == target[1] {
            if typ != b'Z' {
                return None;
            }
            let start = pos;
            while pos < data.len() && data[pos] != 0 {
                pos += 1;
            }
            return Some(&data[start..pos]);
        }

        // Skip this tag's value.
        match typ {
            b'A' | b'c' | b'C' => pos += 1,
            b's' | b'S' => pos += 2,
            b'i' | b'I' | b'f' => pos += 4,
            b'd' => pos += 8,
            b'Z' | b'H' => {
                while pos < data.len() && data[pos] != 0 {
                    pos += 1;
                }
                pos += 1; // null terminator
            }
            b'B' => {
                if pos + 5 > data.len() {
                    return None;
                }
                let elem_size = match data[pos] {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return None,
                };
                let count = u32::from_le_bytes([
                    data[pos + 1],
                    data[pos + 2],
                    data[pos + 3],
                    data[pos + 4],
                ]) as usize;
                pos = (count.checked_mul(elem_size))
                    .and_then(|n| n.checked_add(5))
                    .and_then(|n| pos.checked_add(n))?;
            }
            _ => return None,
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use color_eyre::Result;
    use rust_htslib::bam::{self, FetchDefinition, Read as _, Record};
    use std::collections::HashMap;

    const L001: &str = "mTet1-PyBr-16h-p1_S1_L001";
    const L002: &str = "mTet1-PyBr-16h-p1_S1_L002";
    const L003: &str = "mTet1-PyBr-16h-p1_S1_L003";
    const L004: &str = "mTet1-PyBr-16h-p1_S1_L004";

    fn filter(tag_values: &[&str]) -> TagRequirement {
        RequireTagsParams { require_tags: tag_values.iter().map(|s| s.parse().unwrap()).collect() }
            .filter()
    }

    /// Read records from test.bam, returning one record per unique RG tag seen.
    fn one_record_per_group() -> Result<Vec<(String, Record)>> {
        let mut reader = bam::IndexedReader::from_path("tests/data/test.bam")?;
        reader.fetch(FetchDefinition::All)?;
        let mut seen = HashMap::new();
        let mut record = Record::new();
        while let Some(result) = reader.read(&mut record) {
            result?;
            if let Ok(bam::record::Aux::String(rg)) = record.aux(b"RG") {
                seen.entry(rg.to_owned()).or_insert_with(|| record.clone());
            }
        }
        Ok(seen.into_iter().collect())
    }

    fn view<'a>(record: &'a Record) -> bam::RecordView<'a> {
        // Safety: the RecordView borrows from the record and does not outlive it.
        unsafe { bam::RecordView::from_raw(record.inner() as *const _) }
    }

    #[test]
    fn all_filter_allows_every_record() -> Result<()> {
        let f = filter(&[]);
        assert!(matches!(f, TagRequirement::All));
        for (_, record) in one_record_per_group()? {
            assert!(f.allows(&view(&record)));
        }
        Ok(())
    }

    #[test]
    fn single_tag_allows_matching_record() -> Result<()> {
        let f = filter(&[&format!("RG={L001}")]);
        let records = one_record_per_group()?;
        let (_, l001) = records.iter().find(|(rg, _)| rg == L001).expect("L001 record");
        assert!(f.allows(&view(l001)));
        Ok(())
    }

    #[test]
    fn single_tag_rejects_non_matching_records() -> Result<()> {
        let f = filter(&[&format!("RG={L001}")]);
        let records = one_record_per_group()?;
        for (rg, record) in &records {
            if rg != L001 {
                assert!(!f.allows(&view(record)), "expected {rg} to be rejected");
            }
        }
        Ok(())
    }

    #[test]
    fn multiple_tags_allows_all_members() -> Result<()> {
        let f = filter(&[&format!("RG={L001}"), &format!("RG={L002}"), &format!("RG={L003}")]);
        let records = one_record_per_group()?;
        for (rg, record) in &records {
            let expected = rg != L004;
            assert_eq!(f.allows(&view(record)), expected, "unexpected result for {rg}");
        }
        Ok(())
    }

    #[test]
    fn record_without_tag_is_rejected_when_filter_active() -> Result<()> {
        let f = filter(&[&format!("RG={L001}")]);
        let records = one_record_per_group()?;
        let (_, base) = records.first().expect("at least one record");
        let mut stripped = base.clone();
        stripped.remove_aux(b"RG").ok();
        assert!(!f.allows(&view(&stripped)));
        Ok(())
    }

    #[test]
    fn record_without_tag_passes_when_no_filter() -> Result<()> {
        let f = filter(&[]);
        let records = one_record_per_group()?;
        let (_, base) = records.first().expect("at least one record");
        let mut stripped = base.clone();
        stripped.remove_aux(b"RG").ok();
        assert!(f.allows(&view(&stripped)));
        Ok(())
    }
}
