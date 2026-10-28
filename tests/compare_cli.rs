mod utils;
use std::collections::HashSet;
use std::fs;

use utils::*;

#[test]
#[ignore = "still some differences"]
fn compare_rastair_1_and_2() -> Result<()> {
    apply_common_filters!();
    let temp_dir = TempDir::new()?;
    let temp_file = temp_dir.path().join("test.bed");

    rastair()
        .args([
            "call",
            "-r",
            "tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--cpgs-only",
            "--bed",
        ])
        .arg(&temp_file)
        .succeeds()
        .wrap_err("running rastair")?;

    let (missing, extra) = compare_bed_files("./tests/data/rastair1.bed", temp_file)?;

    if !missing.is_empty() || !extra.is_empty() {
        if !missing.is_empty() {
            eprintln!("Lines in rastair1.bed but not in rastair2.bed:");
            for line in &missing {
                eprintln!("{}", line);
            }
        }
        if !extra.is_empty() {
            eprintln!("Lines in rastair2.bed but not in rastair1.bed:");
            for line in &extra {
                eprintln!("{}", line);
            }
        }
        panic!("BED files differ");
    }

    Ok(())
}

/// Compare BED files ignoring de-novos
fn compare_bed_files<P1, P2>(old_file: P1, new_file: P2) -> Result<(Vec<String>, Vec<String>)>
where
    P1: AsRef<std::path::Path>,
    P2: AsRef<std::path::Path>,
{
    let old_file = old_file.as_ref();
    let new_file = new_file.as_ref();

    let normalize_row = |line: &str, max_cols: usize| -> String {
        let cols: Vec<&str> = line.split('\t').collect();
        cols.iter().take(max_cols).copied().collect::<Vec<_>>().join("\t")
    };

    let is_denovo =
        |line: &str| -> bool { line.rsplit('\t').next().expect("content").contains("NEW") };

    let old_content = fs::read_to_string(old_file)
        .wrap_err_with(|| format!("Failed to read old file: {old_file:?}"))?;
    let new_content = fs::read_to_string(new_file)
        .wrap_err_with(|| format!("Failed to read new file: {new_file:?}"))?;

    let old_lines: Vec<_> = old_content.lines().filter(|l| !l.starts_with('#')).collect();
    let new_lines: Vec<_> = new_content.lines().filter(|l| !l.starts_with('#')).collect();

    let old_col_count = old_lines.first().map(|l| l.split('\t').count()).unwrap_or(0);

    let old_set: HashSet<String> = old_lines.into_iter().map(|l| l.to_string()).collect();
    let new_normalized: HashSet<String> = new_lines
        .into_iter()
        .filter(|l| !is_denovo(l))
        .map(|l| normalize_row(l, old_col_count))
        .collect();

    let in_old_not_new: Vec<String> = old_set.difference(&new_normalized).cloned().collect();
    let in_new_not_old: Vec<String> = new_normalized.difference(&old_set).cloned().collect();

    Ok((in_old_not_new, in_new_not_old))
}
