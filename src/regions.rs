//! Genomic intervals an analysis is restricted to.

use color_eyre::{
    Result,
    eyre::{Context as _, ensure},
};
use rustc_hash::FxHashMap;
use seqair_types::SmolStr;
use std::path::Path;
use tracing::{info, instrument, warn};

/// Half-open intervals that an analysis is restricted to, keyed by contig.
///
/// Used both to score only where a truth set makes a claim (`rastair verify`)
/// and to *train* only there (`rastair ml train`). The latter matters more than
/// it looks: outside a high-confidence BED a truth VCF asserts nothing, so a real
/// variant there is simply absent from it and would be labelled negative. Those
/// mislabelled negatives are not spread evenly — the excluded regions are the
/// repetitive ones, which is exactly where indels concentrate.
///
/// Intervals are sorted and merged on load so that membership is a binary search
/// and overlapping input lines cannot make the search miss a hit.
#[derive(Debug, Default)]
pub struct ConfidentRegions {
    by_contig: FxHashMap<SmolStr, Vec<(u64, u64)>>,
}

impl ConfidentRegions {
    #[instrument(level = "info", skip_all, fields(path = %path.display()))]
    pub fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .wrap_err_with(|| format!("Failed to read BED file: {}", path.display()))?;

        let mut by_contig: FxHashMap<SmolStr, Vec<(u64, u64)>> = FxHashMap::default();
        for (lineno, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }
            let mut fields = line.split('\t');
            let (Some(chrom), Some(start), Some(end)) =
                (fields.next(), fields.next(), fields.next())
            else {
                warn!(line = lineno + 1, "Skipping BED line with fewer than 3 columns");
                continue;
            };
            let (Ok(start), Ok(end)) = (start.parse::<u64>(), end.parse::<u64>()) else {
                warn!(line = lineno + 1, "Skipping BED line with unparsable coordinates");
                continue;
            };
            if end <= start {
                continue;
            }
            by_contig.entry(SmolStr::from(chrom)).or_default().push((start, end));
        }

        for intervals in by_contig.values_mut() {
            intervals.sort_unstable();
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(intervals.len());
            for &(start, end) in intervals.iter() {
                match merged.last_mut() {
                    Some(last) if start <= last.1 => last.1 = last.1.max(end),
                    _ => merged.push((start, end)),
                }
            }
            *intervals = merged;
        }

        let total: usize = by_contig.values().map(Vec::len).sum();
        ensure!(total > 0, "BED file `{}` contains no usable intervals", path.display());
        info!(contigs = by_contig.len(), intervals = total, "Loaded confident regions");

        Ok(Self { by_contig })
    }

    /// Whether `pos` (0-based, as `bcf::Record::pos` reports it) is covered.
    pub fn contains(&self, chrom: &str, pos: u64) -> bool {
        let Some(intervals) = self.by_contig.get(chrom) else {
            return false;
        };
        intervals
            .binary_search_by(|&(start, end)| {
                if end <= pos {
                    std::cmp::Ordering::Less
                } else if start > pos {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .is_ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn regions(intervals: &[(&str, u64, u64)]) -> ConfidentRegions {
        let mut by_contig: FxHashMap<SmolStr, Vec<(u64, u64)>> = FxHashMap::default();
        for &(chrom, start, end) in intervals {
            by_contig.entry(SmolStr::from(chrom)).or_default().push((start, end));
        }
        for v in by_contig.values_mut() {
            v.sort_unstable();
        }
        ConfidentRegions { by_contig }
    }

    #[test]
    fn confident_regions_are_half_open() {
        let r = regions(&[("chr1", 100, 200)]);
        assert!(!r.contains("chr1", 99));
        assert!(r.contains("chr1", 100));
        assert!(r.contains("chr1", 199));
        assert!(!r.contains("chr1", 200), "BED end is exclusive");
    }

    /// A contig absent from the BED is outside the confident set, not inside it —
    /// otherwise restricting to a chr12 BED would silently score all of chr1 too.
    #[test]
    fn unlisted_contig_is_not_confident() {
        let r = regions(&[("chr1", 100, 200)]);
        assert!(!r.contains("chr2", 150));
    }

    #[test]
    fn membership_finds_the_right_interval_among_many() {
        let r = regions(&[("chr1", 0, 10), ("chr1", 100, 200), ("chr1", 1000, 1010)]);
        for pos in [0, 9, 100, 199, 1000, 1009] {
            assert!(r.contains("chr1", pos), "{pos} should be covered");
        }
        for pos in [10, 99, 200, 999, 1010] {
            assert!(!r.contains("chr1", pos), "{pos} should not be covered");
        }
    }
}
