mod utils;
use utils::*;

#[test]
fn verify_help_shows_expected_options() -> Result<()> {
    let output = rastair().args(["verify", "--help"]).output()?;
    let stdout = output.stdout();
    assert!(stdout.contains("--truth"), "missing --truth in help");
    assert!(stdout.contains("--competitor"), "missing --competitor in help");
    assert!(stdout.contains("--region"), "missing --region in help");
    assert!(stdout.contains("--output-json"), "missing --output-json in help");
    assert!(stdout.contains("--output-html"), "missing --output-html in help");
    Ok(())
}

#[test]
fn verify_requires_truth_or_competitor() -> Result<()> {
    // Running verify with only predictions and no --truth or --competitor must fail
    let output = rastair().args(["verify", "nonexistent.vcf"]).output()?;
    assert!(!output.status.success(), "verify should fail without --truth or --competitor");
    let stderr = output.stderr();
    // Should mention the requirement
    assert!(
        stderr.contains("truth") || stderr.contains("competitor"),
        "error message should mention truth or competitor, got: {stderr}"
    );
    Ok(())
}

/// Creates a minimal Rastair-format BCF file in `dir` using rust-htslib's Writer,
/// returning the path. The file is NOT indexed (no regions filtering needed for this test).
fn create_minimal_bcf(dir: &TempDir, name: &str) -> Result<std::path::PathBuf> {
    use cstr8::cstr8;
    use rust_htslib::bcf::{Format, Writer, header::Header};

    let path = dir.path().join(name);

    let mut header = Header::new();
    header.push_record(b"##fileformat=VCFv4.2");
    header.push_record(b"##contig=<ID=chr1,length=248956422>");
    header.push_record(b"##FORMAT=<ID=M5mC,Number=.,Type=Float,Description=\"Methylation level\">");
    header.push_record(b"##FILTER=<ID=PASS,Description=\"All filters passed\">");
    header.push_sample(b"sample");

    let mut writer = Writer::from_path(&path, &header, false, Format::Bcf)
        .wrap_err("failed to create BCF writer")?;

    for (pos, beta) in [(100u64, 0.85f32), (200, 0.2), (300, 0.95)] {
        let mut record = writer.empty_record();
        let rid = writer.header().name2rid(b"chr1").wrap_err("chr1 not in header")?;
        record.set_rid(Some(rid));
        record.set_pos(pos as i64);
        record.set_alleles(&[b"C", b"T"]).wrap_err("set alleles")?;
        record.set_filters(&["PASS".as_bytes()]).wrap_err("set filter")?;
        record.push_format_float(cstr8!("M5mC"), &[beta]).wrap_err("push M5mC")?;
        writer.write(&record).wrap_err("write record")?;
    }

    Ok(path)
}

#[test]
fn verify_html_output_is_written() -> Result<()> {
    let dir = TempDir::new()?;
    let vcf_path = create_minimal_bcf(&dir, "test.bcf")?;
    let html_path = dir.path().join("report.html");

    let output = rastair()
        .args([
            "verify",
            vcf_path.to_str().unwrap(),
            "--competitor",
            vcf_path.to_str().unwrap(),
            "--output-html",
            html_path.to_str().unwrap(),
        ])
        .output()?;

    let stderr = output.stderr();
    assert!(output.status.success(), "verify failed: {stderr}");
    assert!(html_path.exists(), "HTML file was not created");

    let html = std::fs::read_to_string(&html_path)?;
    assert!(html.contains("Rastair Verify Report"), "title missing from HTML");
    assert!(html.contains("density"), "density chart section missing from HTML");
    assert!(html.contains("methylation"), "methylation section key missing from HTML");

    Ok(())
}
