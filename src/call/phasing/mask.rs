use color_eyre::eyre::{Result, WrapErr as _, ensure};
use std::{collections::BTreeMap, io::BufReader, ops::Range, path::PathBuf};

#[derive(Debug, Default)]
pub struct PhaseMask {
    intervals: BTreeMap<String, Vec<Range<u64>>>,
}

impl PhaseMask {
    /// BED coordinates are zero-based and half-open. Plain, gzip and BGZF files
    /// are accepted; the UCSC centromere table's `Chr Start End` header is allowed.
    pub fn load(paths: &[PathBuf]) -> Result<Self> {
        let mut mask = Self::default();
        for path in paths {
            let reader = rust_htslib::bgzf::Reader::from_path(path)
                .wrap_err_with(|| format!("Cannot open phasing mask {}", path.display()))?;
            mask.read(BufReader::new(reader))
                .wrap_err_with(|| format!("Invalid phasing mask {}", path.display()))?;
        }
        for intervals in mask.intervals.values_mut() {
            intervals.sort_by_key(|interval| (interval.start, interval.end));
            let mut merged: Vec<Range<u64>> = Vec::with_capacity(intervals.len());
            for interval in intervals.drain(..) {
                if let Some(previous) = merged.last_mut()
                    && interval.start <= previous.end
                {
                    previous.end = previous.end.max(interval.end);
                } else {
                    merged.push(interval);
                }
            }
            *intervals = merged;
        }
        Ok(mask)
    }

    fn read(&mut self, source: impl std::io::BufRead) -> Result<()> {
        for (index, line) in source.lines().enumerate() {
            let line = line.wrap_err_with(|| format!("Cannot read line {}", index + 1))?;
            let mut fields = line.split_whitespace();
            let Some(chrom) = fields.next() else { continue };
            if chrom.starts_with('#') || matches!(chrom, "track" | "browser") {
                continue;
            }
            let start = fields.next();
            let end = fields.next();
            if chrom == "Chr" && start == Some("Start") && end == Some("End") {
                continue;
            }
            let parse = || -> Result<Range<u64>> {
                use color_eyre::eyre::OptionExt as _;
                let start = start.ok_or_eyre("Missing interval start")?.parse::<u64>()?;
                let end = end.ok_or_eyre("Missing interval end")?.parse::<u64>()?;
                ensure!(start < end, "Expected 0 <= start < end");
                Ok(start..end)
            };
            self.intervals
                .entry(chrom.to_owned())
                .or_default()
                .push(parse().wrap_err_with(|| format!("Invalid interval on line {}", index + 1))?);
        }
        Ok(())
    }

    pub(crate) fn unmasked(&self, chrom: &str, region: Range<u64>) -> Vec<Range<u64>> {
        let mut result = Vec::new();
        let mut start = region.start;
        if let Some(intervals) = self.intervals.get(chrom) {
            let first = intervals.partition_point(|interval| interval.end <= start);
            for interval in intervals.iter().skip(first).take_while(|i| i.start < region.end) {
                if start < interval.start {
                    result.push(start..interval.start.min(region.end));
                }
                start = start.max(interval.end);
            }
        }
        if start < region.end {
            result.push(start..region.end);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    #[test]
    fn merges_files_and_respects_half_open_boundaries() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let centromere = dir.path().join("centromeres.txt");
        std::fs::write(&centromere, "Chr\tStart\tEnd\tGapType\nchr1\t20\t30\tcentromere\n")?;
        let bed = dir.path().join("blacklist.bed.gz");
        let mut writer = rust_htslib::bgzf::Writer::from_path(&bed)?;
        writeln!(writer, "# blacklist\nchr1\t40\t50\tlow map\nchr1\t10\t20\nchr1\t45\t60")?;
        drop(writer);
        let mask = PhaseMask::load(&[centromere, bed])?;
        assert_eq!(mask.unmasked("chr1", 0..70), vec![0..10, 30..40, 60..70]);
        assert_eq!(mask.unmasked("chr1", 10..30), Vec::<Range<u64>>::new());
        assert_eq!(mask.unmasked("chr1", 30..40), vec![30..40]);
        assert_eq!(mask.unmasked("chr2", 10..30), vec![10..30]);
        Ok(())
    }

    #[test]
    fn invalid_intervals_and_missing_files_are_errors() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("mask.bed");
        assert!(PhaseMask::load(std::slice::from_ref(&path)).is_err());
        for contents in ["chr1 2 2", "chr1 -1 10", "chr1 20 10", "chr1 5", "chr1 0 invalid"] {
            std::fs::write(&path, contents)?;
            assert!(PhaseMask::load(std::slice::from_ref(&path)).is_err(), "{contents}");
        }
        Ok(())
    }
}
