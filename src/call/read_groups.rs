use crate::utils::cli;
use rastair_types::SmolStr;
use rust_htslib::bam;
use std::sync::Arc;

#[derive(Debug, clap::Args, Clone, Default, serde::Serialize)]
pub struct ReadGroupsParams {
    /// The read group(s) to filter reads by
    ///
    /// Accepts one or more SAM RG tag values, space-separated. Can also be used with shell
    /// command substitution, e.g. `--read-groups $(cat groups.txt)`.
    #[arg(long, num_args(1..))]
    #[arg(help_heading = cli::sections::INPUT)]
    pub read_groups: Vec<SmolStr>,
}

/// A pre-processed, cheaply cloneable form of the read group filter.
#[derive(Debug, Clone, Default)]
pub enum ReadGroupFilter {
    /// No filter specified, so all reads pass.
    #[default]
    All,
    /// Only reads whose RG tag is in this set pass.
    Groups(Arc<[SmolStr]>),
}

impl ReadGroupsParams {
    pub fn filter(&self) -> ReadGroupFilter {
        if self.read_groups.is_empty() {
            ReadGroupFilter::All
        } else {
            ReadGroupFilter::Groups(self.read_groups.as_slice().into())
        }
    }
}

impl ReadGroupFilter {
    /// Returns `true` if a read with the given RG tag value should be included.
    pub fn allows(&self, record: &bam::RecordView<'_>) -> bool {
        match self {
            ReadGroupFilter::All => true,
            ReadGroupFilter::Groups(groups) => {
                if let Ok(bam::record::Aux::String(rg_str)) = record.aux(b"RG") {
                    groups.iter().any(|g| g.as_str() == rg_str)
                } else {
                    false
                }
            }
        }
    }
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

    fn filter(groups: &[&str]) -> ReadGroupFilter {
        ReadGroupsParams { read_groups: groups.iter().map(SmolStr::new).collect() }.filter()
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

    fn view(record: &Record) -> bam::RecordView<'_> {
        // Safety: the RecordView borrows from the record and does not outlive it.
        unsafe { bam::RecordView::from_raw(record.inner() as *const _) }
    }

    #[test]
    fn all_filter_allows_every_record() -> Result<()> {
        let f = filter(&[]);
        assert!(matches!(f, ReadGroupFilter::All));
        for (_, record) in one_record_per_group()? {
            assert!(f.allows(&view(&record)));
        }
        Ok(())
    }

    #[test]
    fn single_group_allows_matching_record() -> Result<()> {
        let f = filter(&[L001]);
        let records = one_record_per_group()?;
        let (_, l001) = records.iter().find(|(rg, _)| rg == L001).expect("L001 record");
        assert!(f.allows(&view(l001)));
        Ok(())
    }

    #[test]
    fn single_group_rejects_non_matching_records() -> Result<()> {
        let f = filter(&[L001]);
        let records = one_record_per_group()?;
        for (rg, record) in &records {
            if rg != L001 {
                assert!(!f.allows(&view(record)), "expected {rg} to be rejected");
            }
        }
        Ok(())
    }

    #[test]
    fn multiple_groups_allows_all_members() -> Result<()> {
        let f = filter(&[L001, L002, L003]);
        let records = one_record_per_group()?;
        for (rg, record) in &records {
            let expected = rg != L004;
            assert_eq!(f.allows(&view(record)), expected, "unexpected result for {rg}");
        }
        Ok(())
    }

    #[test]
    fn record_without_rg_tag_is_rejected_when_filter_active() -> Result<()> {
        let f = filter(&[L001]);
        // Construct a record with no RG tag by cloning one and stripping it.
        let records = one_record_per_group()?;
        let (_, base) = records.first().expect("at least one record");
        let mut stripped = base.clone();
        stripped.remove_aux(b"RG").ok();
        assert!(!f.allows(&view(&stripped)));
        Ok(())
    }

    #[test]
    fn record_without_rg_tag_passes_when_no_filter() -> Result<()> {
        let f = filter(&[]);
        let records = one_record_per_group()?;
        let (_, base) = records.first().expect("at least one record");
        let mut stripped = base.clone();
        stripped.remove_aux(b"RG").ok();
        assert!(f.allows(&view(&stripped)));
        Ok(())
    }
}
