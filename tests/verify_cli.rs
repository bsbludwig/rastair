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
    assert!(stdout.contains("--match-mode"), "missing --match-mode in help");
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

/// Pull the first value column out of one of `verify`'s markdown tables, e.g.
/// `metric(&stdout, "Precision")`. The tables are column-padded, so matching on the
/// rendered row text directly is too brittle.
#[expect(clippy::panic, reason = "test code")]
fn metric(stdout: &str, label: &str) -> String {
    stdout
        .lines()
        .filter_map(|line| {
            let mut cells = line.split('|').map(str::trim);
            cells.next()?;
            (cells.next()? == label).then(|| cells.next().map(str::to_owned))?
        })
        .next()
        .unwrap_or_else(|| panic!("no `{label}` row in:\n{stdout}"))
}

/// Writes a VCF with one record per `(pos, ref, alt, gt)`, where `gt` is the
/// genotype's allele indices. Used to exercise genotype-aware matching end to end.
fn create_genotyped_vcf(
    dir: &TempDir,
    name: &str,
    records: &[(u64, &str, &str, [i32; 2])],
) -> Result<std::path::PathBuf> {
    use rust_htslib::bcf::{
        Format, Writer,
        header::Header,
        record::GenotypeAllele::{Phased, Unphased},
    };

    let path = dir.path().join(name);

    let mut header = Header::new();
    header.push_record(b"##fileformat=VCFv4.2");
    header.push_record(b"##contig=<ID=chr1,length=248956422>");
    header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    header.push_record(b"##FILTER=<ID=PASS,Description=\"All filters passed\">");
    header.push_sample(b"sample");

    let mut writer = Writer::from_path(&path, &header, true, Format::Vcf)
        .wrap_err("failed to create VCF writer")?;

    for (pos, reference, alt, gt) in records {
        let mut record = writer.empty_record();
        let rid = writer.header().name2rid(b"chr1").wrap_err("chr1 not in header")?;
        record.set_rid(Some(rid));
        record.set_pos(*pos as i64);
        record.set_alleles(&[reference.as_bytes(), alt.as_bytes()]).wrap_err("set alleles")?;
        record.set_filters(&["PASS".as_bytes()]).wrap_err("set filter")?;
        record.push_genotypes(&[Unphased(gt[0]), Phased(gt[1])]).wrap_err("push genotypes")?;
        writer.write(&record).wrap_err("write record")?;
    }

    Ok(path)
}

/// A het call of a hom-alt truth allele is a genotype error, not a true positive —
/// and `--match-mode allele` reverts to scoring it as a match.
#[test]
fn verify_scores_genotype_mismatches_as_errors() -> Result<()> {
    let dir = TempDir::new()?;
    let truth = create_genotyped_vcf(
        &dir,
        "truth.vcf",
        &[
            (100, "C", "T", [1, 1]), // hom alt
            (200, "G", "A", [0, 1]), // het
        ],
    )?;
    // Same alleles, but chr1:100 is called het instead of hom alt.
    let predictions = create_genotyped_vcf(
        &dir,
        "pred.vcf",
        &[(100, "C", "T", [0, 1]), (200, "G", "A", [0, 1])],
    )?;

    let gt_run = rastair()
        .args(["verify", predictions.to_str().unwrap(), "--truth", truth.to_str().unwrap()])
        .output()?;
    assert!(gt_run.status.success(), "verify failed: {}", gt_run.stderr());
    let stdout = gt_run.stdout();
    assert!(stdout.contains("wrong zygosity"), "genotype breakdown missing from report:\n{stdout}");
    assert_eq!(metric(&stdout, "Precision"), "0.5000", "one of two calls should score:\n{stdout}");
    assert_eq!(metric(&stdout, "FP: wrong zygosity"), "1", "in:\n{stdout}");
    assert_eq!(metric(&stdout, "TP"), "1", "in:\n{stdout}");
    assert_eq!(metric(&stdout, "FN"), "1", "a zygosity error misses the truth call:\n{stdout}");

    let allele_run = rastair()
        .args([
            "verify",
            predictions.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "--match-mode",
            "allele",
        ])
        .output()?;
    assert!(allele_run.status.success(), "verify failed: {}", allele_run.stderr());
    let stdout = allele_run.stdout();
    assert_eq!(
        metric(&stdout, "Precision"),
        "1.0000",
        "allele mode should ignore zygosity:\n{stdout}"
    );
    assert!(
        !stdout.contains("wrong zygosity"),
        "allele mode should not report a genotype breakdown:\n{stdout}"
    );

    Ok(())
}

/// Indels padded with an extra anchor base describe the same event and must match.
#[test]
fn verify_matches_indels_across_representations() -> Result<()> {
    let dir = TempDir::new()?;
    // chr1:100 TC>TCAA is the same insertion as chr1:101 C>CAA.
    let truth = create_genotyped_vcf(&dir, "truth.vcf", &[(101, "C", "CAA", [0, 1])])?;
    let predictions = create_genotyped_vcf(&dir, "pred.vcf", &[(100, "TC", "TCAA", [0, 1])])?;

    let output = rastair()
        .args([
            "verify",
            predictions.to_str().unwrap(),
            "--truth",
            truth.to_str().unwrap(),
            "--experimental-indels",
        ])
        .output()?;
    assert!(output.status.success(), "verify failed: {}", output.stderr());
    let stdout = output.stdout();
    assert_eq!(
        (metric(&stdout, "Precision"), metric(&stdout, "Recall")),
        ("1.0000".to_owned(), "1.0000".to_owned()),
        "padded indel should match the trimmed truth record:\n{stdout}"
    );

    Ok(())
}

/// Alt alleles the genotype does not carry are not calls, so they must not count as
/// false positives — rastair emits those for methylation evidence.
#[test]
fn verify_ignores_alt_alleles_the_genotype_does_not_carry() -> Result<()> {
    let dir = TempDir::new()?;
    let truth = create_genotyped_vcf(&dir, "truth.vcf", &[(100, "C", "T", [0, 1])])?;
    let predictions = create_genotyped_vcf(
        &dir,
        "pred.vcf",
        &[
            (100, "C", "T", [0, 1]),
            (300, "C", "T", [0, 0]), // methylation evidence only, not a variant
        ],
    )?;

    let output = rastair()
        .args(["verify", predictions.to_str().unwrap(), "--truth", truth.to_str().unwrap()])
        .output()?;
    assert!(output.status.success(), "verify failed: {}", output.stderr());
    let stdout = output.stdout();
    assert_eq!(
        metric(&stdout, "Precision"),
        "1.0000",
        "0/0 alt should not count as a false positive:\n{stdout}"
    );
    assert!(
        stdout.contains("does not carry"),
        "report should say how many alleles it ignored:\n{stdout}"
    );

    Ok(())
}
