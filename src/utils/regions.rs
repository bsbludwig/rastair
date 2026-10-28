use std::{path::Path, str::FromStr};

use color_eyre::{
    Result,
    eyre::{Context as _, ContextCompat as _, bail, ensure},
};
use seqair_types::{Pos1, RegionString, SmolStr};

/// A parsed CLI region input, resolved at argument parse time.
///
/// Accepts either:
/// - A space-separated list of region strings: `chr`, `chr:start`, `chr:start-end`
/// - A single `@path` reference to a BED file (only when `@` is the first character): `@targets.bed`
///
/// When `@` appears anywhere other than the start of the argument, it is treated as
/// part of a chromosome name rather than a BED file reference.
#[derive(Debug, Clone)]
pub struct CliRegionInput(Vec<RegionString>);

impl CliRegionInput {
    /// Create a `CliRegionInput` from a single [`RegionString`].
    ///
    /// Useful in tests and code that already has a parsed region.
    pub fn from_region(region: RegionString) -> Self {
        CliRegionInput(vec![region])
    }

    /// Access the resolved regions.
    pub fn regions(&self) -> &[RegionString] {
        &self.0
    }
}

impl FromStr for CliRegionInput {
    type Err = color_eyre::eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();
        ensure!(!s.is_empty(), "region argument must not be empty");

        if let Some(bed_path) = s.strip_prefix('@') {
            // Special case @file: read regions from BED file
            ensure!(!bed_path.is_empty(), "empty path after '@' in region argument");
            read_bed_file(bed_path)
                .map(CliRegionInput)
                .wrap_err_with(|| format!("Failed to read regions from BED file: {}", bed_path))
        } else {
            let tokens: Vec<&str> = s.split_ascii_whitespace().collect();
            ensure!(!tokens.is_empty(), "region argument must not be empty");
            let regions = tokens
                .iter()
                .map(|t| RegionString::from_str(t))
                .collect::<Result<Vec<RegionString>, _>>()?;
            Ok(CliRegionInput(regions))
        }
    }
}

/// Read a standard BED file and convert each record to a [`RegionString`].
///
/// BED uses 0-based half-open coordinates: `[start, end)`.
/// We convert to 1-based inclusive: `[start + 1, end]`.
fn read_bed_file(path: &str) -> Result<Vec<RegionString>> {
    read_bed_intervals(Path::new(path))?
        .into_iter()
        .map(|(chrom, start, end)| {
            let pos = |n: u64| {
                u32::try_from(n).ok().and_then(Pos1::new).wrap_err_with(|| {
                    format!("BED position {n} on {chrom} is out of range for a 1-based position")
                })
            };
            Ok(RegionString {
                start: Some(pos(start + 1)?),
                end: Some(pos(end)?),
                chromosome: chrom,
            })
        })
        .collect()
}

/// The `(contig, start, end)` intervals of a BED file, 0-based half-open as
/// written. Lines starting with `#`, `track`, or `browser` are meta/comment
/// lines; extra columns beyond the first three are ignored; empty intervals
/// are skipped with a warning, and a line without three columns is an error.
pub fn read_bed_intervals(path: &Path) -> Result<Vec<(SmolStr, u64, u64)>> {
    use std::{
        fs::File,
        io::{BufRead as _, BufReader},
    };

    ensure!(path.exists(), "BED file does not exist");

    let file = File::open(path).wrap_err("Failed to open BED file")?;
    let mut intervals = Vec::new();
    let mut skipped_empty = 0u64;

    for (line_num, line) in BufReader::new(file).lines().enumerate() {
        let line_num = line_num + 1; // 1-based for error messages
        let line = line.wrap_err_with(|| format!("Failed to read line {line_num} of BED file"))?;
        let line = line.trim();

        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }

        let mut fields = line.splitn(4, '\t');
        let (Some(chrom), Some(start_str), Some(end_str)) =
            (fields.next().filter(|s| !s.is_empty()), fields.next(), fields.next())
        else {
            bail!("BED line {line_num} does not have three tab-separated fields: {line:?}");
        };

        let start: u64 = start_str.parse().wrap_err_with(|| {
            format!("Invalid start position on BED line {line_num}: {start_str:?}")
        })?;
        let end: u64 = end_str.parse().wrap_err_with(|| {
            format!("Invalid end position on BED line {line_num}: {end_str:?}")
        })?;

        if end <= start {
            tracing::warn!("BED line {line_num}: empty interval [{start}, {end}), skipping");
            skipped_empty += 1;
            continue;
        }

        intervals.push((SmolStr::from(chrom), start, end));
    }

    if skipped_empty > 0 {
        tracing::warn!("Skipped {skipped_empty} empty BED record(s)");
    }

    ensure!(!intervals.is_empty(), "BED file {} contains no valid records", path.display());

    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_single_region() {
        let input: CliRegionInput = "chr1".parse().unwrap();
        assert_eq!(input.regions().len(), 1);
        assert_eq!(input.regions()[0].to_string(), "chr1");
    }

    #[test]
    fn parse_multiple_regions() {
        let input: CliRegionInput = "chr1 chr2:100-200 chr3:500".parse().unwrap();
        assert_eq!(input.regions().len(), 3);
        assert_eq!(input.regions()[2].to_string(), "chr3:500");
    }

    #[test]
    fn parse_bed_file() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bed_path = dir.path().join("targets.bed");
        std::fs::write(&bed_path, "chr1\t99\t199\nchr2\t0\t500\n")?;

        let input: CliRegionInput = format!("@{}", bed_path.display()).parse()?;
        assert_eq!(input.regions().len(), 2);
        assert_eq!(input.regions()[0].to_string(), "chr1:100-199");
        assert_eq!(input.regions()[1].to_string(), "chr2:1-500");

        Ok(())
    }

    #[test]
    fn parse_bed_file_with_comment_lines() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bed_path = dir.path().join("targets.bed");
        std::fs::write(
            &bed_path,
            "# header comment\nchr3\t1000\t2000\n# mid comment\nchr4\t0\t100\n",
        )?;

        let input: CliRegionInput = format!("@{}", bed_path.display()).parse()?;
        assert_eq!(input.regions().len(), 2);
        assert_eq!(input.regions()[0].to_string(), "chr3:1001-2000");

        Ok(())
    }

    #[test]
    fn at_sign_only_special_as_prefix() {
        // @ only triggers BED file mode when it starts the whole argument.
        // Inside a space-separated list, `@targets.bed` is a chromosome name.
        let input: CliRegionInput = "chr1 @targets.bed".parse().unwrap();
        assert_eq!(input.regions().len(), 2);
        assert_eq!(input.regions()[0].to_string(), "chr1");
        assert_eq!(input.regions()[1].to_string(), "@targets.bed");
    }

    #[test]
    fn empty_input_is_error() {
        let err = "  ".parse::<CliRegionInput>().unwrap_err();
        assert!(err.to_string().contains("must not be empty"));
    }

    #[test]
    fn empty_at_path_is_error() {
        let err = "@".parse::<CliRegionInput>().unwrap_err();
        assert!(err.to_string().contains("empty path after '@'"));
    }

    #[test]
    fn invalid_region_is_error() {
        let err = ":100".parse::<CliRegionInput>().unwrap_err();
        assert!(err.to_string().contains("Malformed region string"));
    }

    #[test]
    fn test_read_bed_file() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bed_path = dir.path().join("test.bed");
        std::fs::write(
            &bed_path,
            "chr1\t99\t199\nchr2\t0\t500\n# comment line\nchr3\t1000\t2000\n",
        )?;

        let regions = read_bed_file(bed_path.to_str().unwrap())?;
        assert_eq!(regions.len(), 3);

        // chr1:99-199 (0-based half-open) → 1-based inclusive chr1:100-199
        assert_eq!(regions[0].chromosome.as_str(), "chr1");
        assert_eq!(regions[0].start, Pos1::new(100));
        assert_eq!(regions[0].end, Pos1::new(199));

        // chr2:0-500 (0-based half-open) → 1-based inclusive chr2:1-500
        assert_eq!(regions[1].chromosome.as_str(), "chr2");
        assert_eq!(regions[1].start, Pos1::new(1));
        assert_eq!(regions[1].end, Pos1::new(500));

        // chr3:1000-2000 (0-based half-open) → 1-based inclusive chr3:1001-2000
        assert_eq!(regions[2].chromosome.as_str(), "chr3");
        assert_eq!(regions[2].start, Pos1::new(1001));
        assert_eq!(regions[2].end, Pos1::new(2000));

        Ok(())
    }

    #[test]
    fn a_bed_line_without_three_fields_is_an_error() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bed_path = dir.path().join("short.bed");
        std::fs::write(&bed_path, "chr1\t99\t199\nchr2\t0\n")?;

        let err = read_bed_file(bed_path.to_str().unwrap()).unwrap_err();
        assert!(err.to_string().contains("BED line 2"), "{err}");

        Ok(())
    }

    #[test]
    fn test_empty_bed_file_is_error() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let bed_path = dir.path().join("empty.bed");
        std::fs::write(&bed_path, "# only a comment\n")?;

        let err = read_bed_file(bed_path.to_str().unwrap()).unwrap_err();
        assert!(err.to_string().contains("no valid records"));

        Ok(())
    }

    #[test]
    fn test_nonexistent_bed_file_is_error() {
        let err = "@nonexistent_file.bed".parse::<CliRegionInput>().unwrap_err();
        let str = format!("{err:#}");
        assert!(str.contains("does not exist"));
    }
}
