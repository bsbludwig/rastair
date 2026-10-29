#[path = "../src/pinned.rs"]
mod pinned;
mod utils;
use rust_htslib::bam::ext::BamRecordExtensions as _;
use rust_htslib::bam::{self, Read as BamRead};
use utils::*;

#[test]
fn bam_rewrite_legacy_adds_xr_xg_xm_tags() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    let output_bam = temp_dir.path().join("output.bam");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6103075-6103100",
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    let tabix_output = std::process::Command::new("tabix")
        .args(["-p", "bed", "-f"])
        .arg(&calls_bed)
        .output()
        .wrap_err("Failed to run tabix")?;
    ensure!(
        tabix_output.status.success(),
        "tabix failed: {}",
        String::from_utf8_lossy(&tabix_output.stderr)
    );

    rastair()
        .args([
            "bam",
            "legacy",
            "--fasta-file=tests/data/test.fasta.gz",
            "--region=chr19:6103075-6103100",
            "tests/data/test.bam",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&output_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite BAM")?;

    let mut bam_reader =
        rust_htslib::bam::Reader::from_path(&output_bam).wrap_err("Failed to open output BAM")?;

    let mut record = rust_htslib::bam::Record::new();
    let mut records_checked = 0;
    let mut flags_seen = std::collections::HashSet::new();

    while let Some(result) = bam_reader.read(&mut record) {
        result.wrap_err("Failed to read BAM record")?;
        records_checked += 1;
        let flag = record.flags();
        flags_seen.insert(flag);

        let xr = record.aux(b"XR").wrap_err(format!("Missing XR tag for flag {}", flag))?;
        let xr_value = match xr {
            rust_htslib::bam::record::Aux::String(s) => s,
            _ => panic!("XR tag is not a string for flag {}", flag),
        };

        let xg = record.aux(b"XG").wrap_err(format!("Missing XG tag for flag {}", flag))?;
        let xg_value = match xg {
            rust_htslib::bam::record::Aux::String(s) => s,
            _ => panic!("XG tag is not a string for flag {}", flag),
        };

        let xm = record.aux(b"XM").wrap_err(format!("Missing XM tag for flag {}", flag))?;
        let xm_value = match xm {
            rust_htslib::bam::record::Aux::String(s) => s,
            _ => panic!("XM tag is not a string for flag {}", flag),
        };

        assert_eq!(
            xm_value.len(),
            record.seq_len(),
            "XM tag length doesn't match sequence length for flag {}",
            flag
        );

        match flag {
            99 => {
                assert_eq!(xr_value, "CT", "XR tag mismatch for flag 99");
                assert_eq!(xg_value, "CT", "XG tag mismatch for flag 99");
            }
            147 => {
                assert_eq!(xr_value, "GA", "XR tag mismatch for flag 147");
                assert_eq!(xg_value, "CT", "XG tag mismatch for flag 147");
            }
            83 => {
                assert_eq!(xr_value, "CT", "XR tag mismatch for flag 83");
                assert_eq!(xg_value, "GA", "XG tag mismatch for flag 83");
            }
            163 => {
                assert_eq!(xr_value, "GA", "XR tag mismatch for flag 163");
                assert_eq!(xg_value, "GA", "XG tag mismatch for flag 163");
            }
            _ => {}
        }

        // Legacy mode should NOT produce MM/ML tags
        assert!(
            record.aux(b"MM").is_err(),
            "Legacy mode should not produce MM tag for flag {}",
            flag
        );
    }

    assert!(records_checked > 0, "No records were checked");
    assert!(
        flags_seen.contains(&83)
            || flags_seen.contains(&99)
            || flags_seen.contains(&147)
            || flags_seen.contains(&163),
        "Expected to see at least one of the standard paired-end flags"
    );

    Ok(())
}

#[test]
fn bam_rewrite_preserves_existing_tags() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    let output_bam = temp_dir.path().join("output.bam");

    // Generate calls
    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6103075-6103100",
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    // Index the calls file
    let tabix_output = std::process::Command::new("tabix")
        .args(["-p", "bed", "-f"])
        .arg(&calls_bed)
        .output()
        .wrap_err("Failed to run tabix")?;
    ensure!(
        tabix_output.status.success(),
        "tabix failed: {}",
        String::from_utf8_lossy(&tabix_output.stderr)
    );

    rastair()
        .args([
            "bam",
            "standard",
            "--fasta-file=tests/data/test.fasta.gz",
            "--region=chr19:6103075-6103100",
            "tests/data/test.bam",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&output_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite BAM")?;

    // Verify existing tags are preserved
    let mut bam_reader =
        rust_htslib::bam::Reader::from_path(&output_bam).wrap_err("Failed to open output BAM")?;

    let mut record = rust_htslib::bam::Record::new();
    let mut records_checked = 0;

    while let Some(result) = bam_reader.read(&mut record) {
        result.wrap_err("Failed to read BAM record")?;
        records_checked += 1;

        // Check that standard tags are preserved
        record.aux(b"NM").wrap_err("Missing NM tag (should be preserved)")?;
        record.aux(b"MD").wrap_err("Missing MD tag (should be preserved)")?;
        record.aux(b"AS").wrap_err("Missing AS tag (should be preserved)")?;
    }

    assert!(records_checked > 0, "No records were checked");

    Ok(())
}

/// A `--region` rewrite is sorted like its input and holds each read of the
/// region once; a duplicate would break the sort
/// order, and `samtools index` would fail.
#[test]
fn bam_rewrite_of_a_region_is_sorted_without_duplicates() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    let output_bam = temp_dir.path().join("output.bam");

    // Generate calls over a wider region so we get enough reads to span
    // multiple segments.
    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6103000-6106000",
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    let tabix_output = std::process::Command::new("tabix")
        .args(["-p", "bed", "-f"])
        .arg(&calls_bed)
        .output()
        .wrap_err("Failed to run tabix")?;
    ensure!(
        tabix_output.status.success(),
        "tabix failed: {}",
        String::from_utf8_lossy(&tabix_output.stderr)
    );

    rastair()
        .args([
            "bam",
            "standard",
            "--fasta-file=tests/data/test.fasta.gz",
            "--region=chr19:6103000-6106000",
            "tests/data/test.bam",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&output_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite BAM")?;

    // Count records in the input BAM for the same region (the ground truth).
    let input_count = {
        let mut bam_in = bam::IndexedReader::from_path("tests/data/test.bam")
            .wrap_err("Failed to open input BAM")?;
        bam_in.fetch(("chr19", 6103000, 6106000))?;
        let mut r = bam::Record::new();
        let mut n = 0u32;
        while let Some(res) = bam_in.read(&mut r) {
            res?;
            n += 1;
        }
        n
    };

    let mut reader = bam::Reader::from_path(&output_bam).wrap_err("Failed to open output BAM")?;

    let mut record = bam::Record::new();
    let mut prev_pos: i64 = -1;
    let mut output_count = 0u32;

    while let Some(result) = reader.read(&mut record) {
        result.wrap_err("Failed to read BAM record")?;
        output_count += 1;

        let pos = record.pos();
        ensure!(pos >= prev_pos, "Output BAM is not sorted: position {prev_pos} followed by {pos}");
        prev_pos = pos;
    }

    ensure!(output_count > 10, "Expected more than 10 records, got {output_count}");
    ensure!(
        output_count == input_count,
        "Record count mismatch: input has {input_count} but output has {output_count} \
         (duplicates introduced or records lost)"
    );

    Ok(())
}

/// Full CLI roundtrip: `rastair call` → `rastair bam legacy` + `standard` → verify XM content.
/// Ensures the XM tag only annotates real CpG positions (not every C/G) and
/// that both legacy and standard modes produce matching methylation calls.
#[test]
fn bam_rewrite_xm_content_matches_calls() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    let legacy_bam = temp_dir.path().join("legacy.bam");
    let standard_bam = temp_dir.path().join("standard.bam");

    let region = "chr19:6103075-6103200";

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            &format!("--region={region}"),
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    let tabix_output = std::process::Command::new("tabix")
        .args(["-p", "bed", "-f"])
        .arg(&calls_bed)
        .output()
        .wrap_err("Failed to run tabix")?;
    ensure!(
        tabix_output.status.success(),
        "tabix failed: {}",
        String::from_utf8_lossy(&tabix_output.stderr)
    );

    rastair()
        .args([
            "bam",
            "legacy",
            "--fasta-file=tests/data/test.fasta.gz",
            &format!("--region={region}"),
            "tests/data/test.bam",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&legacy_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite legacy BAM")?;

    rastair()
        .args([
            "bam",
            "standard",
            "--fasta-file=tests/data/test.fasta.gz",
            &format!("--region={region}"),
            "tests/data/test.bam",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&standard_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite standard BAM")?;

    let mut legacy_reader = bam::Reader::from_path(&legacy_bam).wrap_err("open legacy BAM")?;
    let mut standard_reader =
        bam::Reader::from_path(&standard_bam).wrap_err("open standard BAM")?;

    let mut legacy_rec = bam::Record::new();
    let mut standard_rec = bam::Record::new();
    let mut records_checked = 0u32;

    loop {
        let legacy_ok = legacy_reader.read(&mut legacy_rec);
        let standard_ok = standard_reader.read(&mut standard_rec);

        match (legacy_ok, standard_ok) {
            (Some(Ok(())), Some(Ok(()))) => {}
            (None, None) => break,
            _ => bail!("Legacy and standard BAMs have different record counts"),
        }

        records_checked += 1;
        let flag = legacy_rec.flags();

        ensure!(
            legacy_rec.pos() == standard_rec.pos() && flag == standard_rec.flags(),
            "Record mismatch at index {records_checked}"
        );

        // --- XM content checks ---
        let xm = match legacy_rec.aux(b"XM") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => s,
            _ => bail!("Missing/bad XM tag at record {records_checked}"),
        };

        ensure!(
            xm.len() == legacy_rec.seq_len(),
            "XM length {} != seq length {} for flag {flag}",
            xm.len(),
            legacy_rec.seq_len(),
        );

        // XM must only contain valid characters
        for c in xm.chars() {
            ensure!(
                matches!(c, '.' | 'z' | 'Z' | 'x' | 'X' | 'h' | 'H'),
                "Invalid XM character '{c}' in record flag={flag}: {xm}"
            );
        }

        // Annotations must be much fewer than total C/G bases (sanity check
        // that we're not marking every C/G as CpG)
        let annotation_count = xm.chars().filter(|c| *c != '.').count();
        let seq = legacy_rec.seq().as_bytes();
        let cg_count = seq.iter().filter(|&&b| b == b'C' || b == b'G').count();
        if cg_count > 4 {
            ensure!(
                annotation_count < cg_count,
                "Too many XM annotations ({annotation_count}) vs C/G bases ({cg_count}) \
                 for flag {flag}. XM: {xm}"
            );
        }

        // --- Cross-check with MM tag ---
        // MM/ML are absent when no reads show methylation evidence (absent = 0 positions).
        let xm_methylated_count = xm.chars().filter(|c| c.is_ascii_uppercase()).count();
        let mm_methylated_count = match standard_rec.aux(b"MM") {
            Ok(rust_htslib::bam::record::Aux::String(mm)) => {
                mm.trim_end_matches(';').split(',').skip(1).filter(|s| !s.is_empty()).count()
            }
            Ok(_) => bail!("MM tag is not a string at record {records_checked}"),
            Err(_) => 0,
        };

        ensure!(
            xm_methylated_count == mm_methylated_count,
            "Methylated count mismatch for flag {flag}: XM has {xm_methylated_count}, \
             MM has {mm_methylated_count}. XM: {xm}"
        );
    }

    ensure!(records_checked > 0, "No records were checked");

    Ok(())
}

/// Mixed-methylation sites in legacy XM output.
///
/// At chr19:6103215 (1-based), the CpG C base is at 0-based position 6103214.
/// In `tests/data/test.bam` there are 18 OT reads with evidence split:
/// 8 methylated (T => `Z`) and 10 unmethylated (C => `z`).
/// Legacy XM must preserve this per-read split even when the site-level beta is < 0.5.
#[test]
fn bam_rewrite_legacy_mixed_site_keeps_per_read_z_and_z_uppercase() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    let output_bam = temp_dir.path().join("legacy.bam");

    let region = "chr19:6103000-6103300";
    // User-facing coordinate is 1-based chr19:6103215; BAM/reference lookup is 0-based.
    let target_pos_0_based = 6_103_214_i64;

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            &format!("--region={region}"),
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    rastair()
        .args([
            "bam",
            "legacy",
            "--fasta-file=tests/data/test.fasta.gz",
            &format!("--region={region}"),
            "tests/data/test.bam",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&output_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite legacy BAM")?;

    let mut reader = bam::Reader::from_path(&output_bam).wrap_err("open legacy BAM")?;
    let mut record = bam::Record::new();

    let mut z_upper = 0usize;
    let mut z_lower = 0usize;

    while let Some(result) = reader.read(&mut record) {
        result.wrap_err("Failed to read BAM record")?;

        let xm = match record.aux(b"XM") {
            Ok(rust_htslib::bam::record::Aux::String(s)) => s,
            _ => bail!("Missing/bad XM tag for record with flag {}", record.flags()),
        };
        ensure!(
            xm.len() == record.seq_len(),
            "XM length {} != seq length {}",
            xm.len(),
            record.seq_len()
        );

        for [pos_in_read, pos_in_ref] in record.aligned_pairs_full() {
            let Some(pos_in_read) = pos_in_read else { continue };
            let Some(pos_in_ref) = pos_in_ref else { continue };

            if pos_in_ref != target_pos_0_based {
                continue;
            }

            let idx =
                usize::try_from(pos_in_read).wrap_err("read position does not fit in usize")?;
            let symbol = xm
                .as_bytes()
                .get(idx)
                .copied()
                .map(char::from)
                .ok_or_else(|| eyre!("XM index out of bounds"))?;

            match symbol {
                'Z' => z_upper += 1,
                'z' => z_lower += 1,
                '.' => {}
                other => bail!("Unexpected XM symbol '{other}' at target position"),
            }
        }
    }

    ensure!(z_upper == 8, "Expected 8 uppercase Z calls at chr19:6103215, got {z_upper}");
    ensure!(z_lower == 10, "Expected 10 lowercase z calls at chr19:6103215, got {z_lower}");
    ensure!(z_upper + z_lower == 18, "Expected 18 CpG annotations at chr19:6103215");

    Ok(())
}

/// `bam legacy` must accept a CRAM input file (using the required FASTA as the
/// decoding reference) and produce a valid BAM with XR/XG/XM tags.
#[test]
fn bam_legacy_accepts_cram_input() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    let output_bam = temp_dir.path().join("output.bam");

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6103075-6103100",
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    std::process::Command::new("tabix")
        .args(["-p", "bed", "-f"])
        .arg(&calls_bed)
        .status()
        .wrap_err("Failed to run tabix")?
        .succeeds()
        .wrap_err("tabix failed")?;

    rastair()
        .args([
            "bam",
            "legacy",
            "--fasta-file=tests/data/test.fasta.gz",
            "--region=chr19:6103075-6103100",
            "tests/data/test.cram",
        ])
        .arg(&calls_bed)
        .args(["-o"])
        .arg(&output_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite CRAM to legacy BAM")?;

    let mut reader = bam::Reader::from_path(&output_bam).wrap_err("Failed to open output BAM")?;
    let mut record = bam::Record::new();
    let mut records_checked = 0u32;

    while let Some(result) = reader.read(&mut record) {
        result.wrap_err("Failed to read record")?;
        records_checked += 1;

        record
            .aux(b"XR")
            .wrap_err_with(|| format!("Missing XR tag for flag {}", record.flags()))?;
        record
            .aux(b"XG")
            .wrap_err_with(|| format!("Missing XG tag for flag {}", record.flags()))?;
        let xm = match record.aux(b"XM") {
            Ok(bam::record::Aux::String(s)) => s,
            _ => bail!("Missing/bad XM tag for flag {}", record.flags()),
        };
        ensure!(
            xm.len() == record.seq_len(),
            "XM length {} != seq length {} for flag {}",
            xm.len(),
            record.seq_len(),
            record.flags()
        );
        // Output must be BAM, not CRAM
        assert!(
            record.aux(b"MM").is_err(),
            "Legacy mode should not produce MM tag (flag {})",
            record.flags()
        );
    }

    ensure!(records_checked > 0, "No records were written to output BAM");
    Ok(())
}

/// CRAM and BAM inputs must produce identical XR/XG/XM tag content for every
/// record when running `bam legacy` over the same region.
#[test]
fn bam_legacy_cram_matches_bam_output() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    let bam_out = temp_dir.path().join("from_bam.bam");
    let cram_out = temp_dir.path().join("from_cram.bam");

    let region = "chr19:6103075-6103200";

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            &format!("--region={region}"),
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    std::process::Command::new("tabix")
        .args(["-p", "bed", "-f"])
        .arg(&calls_bed)
        .status()
        .wrap_err("Failed to run tabix")?
        .succeeds()
        .wrap_err("tabix failed")?;

    for (input, output) in [("tests/data/test.bam", &bam_out), ("tests/data/test.cram", &cram_out)]
    {
        rastair()
            .args([
                "bam",
                "legacy",
                "--fasta-file=tests/data/test.fasta.gz",
                &format!("--region={region}"),
                input,
            ])
            .arg(&calls_bed)
            .args(["-o"])
            .arg(output)
            .silent()
            .succeeds()
            .wrap_err_with(|| format!("Failed to rewrite {input} to legacy BAM"))?;
    }

    let mut bam_reader = bam::Reader::from_path(&bam_out).wrap_err("open BAM-sourced output")?;
    let mut cram_reader = bam::Reader::from_path(&cram_out).wrap_err("open CRAM-sourced output")?;

    let mut bam_rec = bam::Record::new();
    let mut cram_rec = bam::Record::new();
    let mut records_compared = 0u32;

    loop {
        let bam_ok = bam_reader.read(&mut bam_rec);
        let cram_ok = cram_reader.read(&mut cram_rec);

        match (bam_ok, cram_ok) {
            (None, None) => break,
            (Some(Ok(())), Some(Ok(()))) => {}
            _ => bail!("BAM-sourced and CRAM-sourced outputs have different record counts"),
        }

        records_compared += 1;
        let flag = bam_rec.flags();

        ensure!(
            bam_rec.pos() == cram_rec.pos() && flag == cram_rec.flags(),
            "Record mismatch at index {records_compared}: pos {}/{}, flag {}/{}",
            bam_rec.pos(),
            cram_rec.pos(),
            flag,
            cram_rec.flags()
        );

        for tag in [b"XR", b"XG", b"XM"] {
            let bam_val = match bam_rec.aux(tag) {
                Ok(bam::record::Aux::String(s)) => s.to_owned(),
                _ => bail!(
                    "Missing/bad {} tag in BAM-sourced output (flag {flag})",
                    String::from_utf8_lossy(tag)
                ),
            };
            let cram_val = match cram_rec.aux(tag) {
                Ok(bam::record::Aux::String(s)) => s.to_owned(),
                _ => bail!(
                    "Missing/bad {} tag in CRAM-sourced output (flag {flag})",
                    String::from_utf8_lossy(tag)
                ),
            };
            ensure!(
                bam_val == cram_val,
                "{} mismatch at record {records_compared} (flag {flag}): BAM={bam_val:?} CRAM={cram_val:?}",
                String::from_utf8_lossy(tag)
            );
        }
    }

    ensure!(records_compared > 0, "No records were compared");
    Ok(())
}

/// A CpG no read shows methylation for has an empty `beta_est` in the calls
/// BED. `bam` and `per-read` still read it as a CpG: per-read tags come from
/// each read's base, not from the position's β.
///
/// The reference CpG at 300 is covered only by top-strand reads that end on its
/// C, so no read pairs the C with its G, and β is empty. The de-novo CpG made
/// by the `T>C` at 600 is covered only by bottom-strand reads, which cannot
/// show the C side's methylation.
#[test]
fn a_cpg_without_beta_stays_a_cpg_downstream() -> Result<()> {
    use color_eyre::eyre::OptionExt as _;
    use utils::scenario::{Hap, READ_LEN, Scenario, Strand};
    const CPG: usize = 300;
    const DENOVO: usize = 600;

    let scenario = Scenario::new(1_000)
        .reference_at(CPG, b"CG")
        .reference_at(DENOVO, b"TG")
        .hom(DENOVO, b'C')
        .methylated(Hap::One)
        .pairs(Hap::One, Strand::Top, CPG + 1 - READ_LEN, 3)
        .pairs(Hap::Two, Strand::Top, CPG + 1 - READ_LEN, 3)
        .pairs(Hap::One, Strand::Bottom, DENOVO - 50, 3)
        .pairs(Hap::Two, Strand::Bottom, DENOVO - 50, 3);
    let dir = TempDir::new()?;
    let fasta = scenario.write_fasta(dir.path())?;
    let reads = dir.path().join("reads.bam");
    scenario.write_bam(&reads)?;
    let calls = dir.path().join("calls.bed.gz");
    rastair()
        .arg("call")
        .arg("--fasta-file")
        .arg(&fasta)
        .arg(&reads)
        .arg(scenario.region())
        .arg("--bed")
        .arg(&calls)
        .silent()
        .succeeds()?;

    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut rust_htslib::bgzf::Reader::from_path(&calls).wrap_err("open calls BED")?,
        &mut text,
    )?;
    let beta = |pos: usize| {
        let pos = pos.to_string();
        text.lines().find(|line| line.split('\t').nth(1) == Some(pos.as_str()))?.split('\t').nth(4)
    };
    ensure!(beta(CPG) == Some(""), "the CpG's β is not empty:\n{text}");
    ensure!(beta(DENOVO) == Some(""), "the de-novo CpG's β is not empty:\n{text}");

    let legacy = dir.path().join("legacy.bam");
    let output = rastair()
        .arg("bam")
        .arg("legacy")
        .arg("--fasta-file")
        .arg(&fasta)
        .arg(&reads)
        .arg(&calls)
        .arg("--output")
        .arg(&legacy)
        .output()?;
    ensure!(output.status.success(), "bam legacy failed: {}", output.stderr());
    ensure!(!output.stderr().contains("WARN"), "bam legacy warned: {}", output.stderr());

    let mut reader = bam::Reader::from_path(&legacy).wrap_err("open legacy BAM")?;
    let mut record = bam::Record::new();
    let (mut methylated, mut unmethylated) = (0, 0);
    while let Some(result) = reader.read(&mut record) {
        result?;
        let Ok(bam::record::Aux::String(xm)) = record.aux(b"XM") else {
            bail!("no XM on {}", record.pos())
        };
        let xm = xm.as_bytes().to_vec();
        for [in_read, in_ref] in record.aligned_pairs_full() {
            if in_ref != Some(i64::try_from(CPG)?) {
                continue;
            }
            let in_read = usize::try_from(in_read.ok_or_eyre("no read position")?)?;
            match xm.get(in_read) {
                Some(b'Z') => methylated += 1,
                Some(b'z') => unmethylated += 1,
                other => bail!("XM at the CpG is {other:?}"),
            }
        }
    }
    ensure!(
        (methylated, unmethylated) == (3, 3),
        "XM at the CpG: {methylated} Z, {unmethylated} z"
    );

    let output = rastair()
        .arg("per-read")
        .arg("--fasta-file")
        .arg(&fasta)
        .arg(&reads)
        .arg("--calls")
        .arg(&calls)
        .output()?;
    ensure!(output.status.success(), "per-read failed: {}", output.stderr());
    ensure!(!output.stderr().contains("WARN"), "per-read warned: {}", output.stderr());
    Ok(())
}

/// β is 0 both at an unmethylated CpG and at one both chromosomes lose, so
/// `bam` must not take a β of 0 for an SNP: the CpG at 300 has no methylated
/// read at all and every read there gets `z`. At 900 one methylated read among
/// 210 rounds β to `0.00`; that read keeps its MM/ML entry.
#[test]
fn an_unmethylated_cpg_keeps_its_tags() -> Result<()> {
    use color_eyre::eyre::OptionExt as _;
    use utils::scenario::{Hap, Scenario, Strand};
    const UNMETHYLATED: usize = 300;
    const NEARLY: usize = 900;

    let scenario = Scenario::new(1_500)
        .reference_at(UNMETHYLATED, b"CG")
        .reference_at(NEARLY, b"CG")
        .methylated(Hap::One)
        .pairs(Hap::Two, Strand::Top, UNMETHYLATED - 50, 3)
        .pairs(Hap::Two, Strand::Bottom, UNMETHYLATED - 50, 3)
        .pairs(Hap::One, Strand::Top, NEARLY - 50, 1)
        .pairs(Hap::Two, Strand::Top, NEARLY - 50, 210)
        .pairs(Hap::Two, Strand::Bottom, NEARLY - 50, 3);
    let dir = TempDir::new()?;
    let fasta = scenario.write_fasta(dir.path())?;
    let reads = dir.path().join("reads.bam");
    scenario.write_bam(&reads)?;
    let calls = dir.path().join("calls.bed.gz");
    rastair()
        .arg("call")
        .arg("--fasta-file")
        .arg(&fasta)
        .arg(&reads)
        .arg(scenario.region())
        .arg("--bed")
        .arg(&calls)
        .silent()
        .succeeds()?;

    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut rust_htslib::bgzf::Reader::from_path(&calls).wrap_err("open calls BED")?,
        &mut text,
    )?;
    let beta = |pos: usize| {
        let pos = pos.to_string();
        text.lines().find(|line| line.split('\t').nth(1) == Some(pos.as_str()))?.split('\t').nth(4)
    };
    for pos in [UNMETHYLATED, UNMETHYLATED + 1, NEARLY] {
        ensure!(beta(pos) == Some("0.00"), "β at {pos} is not 0.00:\n{text}");
    }

    let rewrite = |mode: &str| -> Result<std::path::PathBuf> {
        let out = dir.path().join(format!("{mode}.bam"));
        let output = rastair()
            .arg("bam")
            .arg(mode)
            .arg("--fasta-file")
            .arg(&fasta)
            .arg(&reads)
            .arg(&calls)
            .arg("--output")
            .arg(&out)
            .output()?;
        ensure!(output.status.success(), "bam {mode} failed: {}", output.stderr());
        Ok(out)
    };

    let mut reader = bam::Reader::from_path(rewrite("legacy")?).wrap_err("open legacy BAM")?;
    let mut record = bam::Record::new();
    let mut at_unmethylated = Vec::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        let Ok(bam::record::Aux::String(xm)) = record.aux(b"XM") else {
            bail!("no XM on {}", record.pos())
        };
        let xm = xm.as_bytes().to_vec();
        for [in_read, in_ref] in record.aligned_pairs_full() {
            let Some(in_ref) = in_ref else { continue };
            if ![UNMETHYLATED, UNMETHYLATED + 1].contains(&usize::try_from(in_ref)?) {
                continue;
            }
            let in_read = usize::try_from(in_read.ok_or_eyre("no read position")?)?;
            at_unmethylated.push(xm.get(in_read).copied().map(char::from));
        }
    }
    let unmethylated = at_unmethylated.iter().filter(|&&xm| xm == Some('z')).count();
    ensure!(
        unmethylated == 6 && !at_unmethylated.contains(&Some('Z')),
        "XM at the unmethylated CpG: {at_unmethylated:?}"
    );

    let mut reader = bam::Reader::from_path(rewrite("standard")?).wrap_err("open standard BAM")?;
    let mut with_mm = 0;
    while let Some(result) = reader.read(&mut record) {
        result?;
        if let Ok(bam::record::Aux::String(mm)) = record.aux(b"MM")
            && mm.trim_end_matches(';').split(',').skip(1).any(|s| !s.is_empty())
        {
            with_mm += 1;
        }
    }
    ensure!(with_mm == 1, "{with_mm} reads have MM entries, expected the one methylated read");
    Ok(())
}

/// `bam` streams the whole input, so every record comes out in input order,
/// the unplaced unmapped reads at its end (`*` in `samtools idxstats`)
/// included.
#[test]
fn bam_keeps_every_record() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    rastair()
        .args(CALL_TEST_BAM)
        .args(["--cpgs-only", "--bed"])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    let names = |path: &Path| -> Result<Vec<(Vec<u8>, u16)>> {
        let mut reader = bam::Reader::from_path(path).wrap_err("Failed to open BAM")?;
        reader
            .records()
            .map(|record| {
                let record = record.wrap_err("Failed to read BAM record")?;
                Ok((record.qname().to_owned(), record.flags()))
            })
            .collect()
    };
    let input = names(Path::new("tests/data/test.bam"))?;
    ensure!(input.len() == 29_358, "the test BAM has {} records", input.len());

    for mode in ["standard", "legacy"] {
        let output_bam = temp_dir.path().join(format!("{mode}.bam"));
        rastair()
            .args(["bam", mode, "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam"])
            .arg(&calls_bed)
            .arg("-o")
            .arg(&output_bam)
            .silent()
            .succeeds()
            .wrap_err_with(|| format!("Failed to rewrite BAM in {mode} mode"))?;
        ensure!(names(&output_bam)? == input, "{mode} should keep every record, in order");
    }
    Ok(())
}

mod generated_inputs {
    //! `bam` on generated inputs: mapped, placed-unmapped and unplaced
    //! records, sorted or not, on either side of a batch boundary.

    use super::*;
    use hegel::{TestCase, generators as gs};
    use seqair::bam::{
        AuxData, BamHeader, CigarOp, RawBamRecord, aligned_pairs::AlignedPairs,
        writer::BamWriterBuilder,
    };
    use seqair::reader::AlignmentStream;
    use seqair_types::Strand;
    use std::fmt::Write as _;
    use std::io::Write as _;

    const CONTIGS: [&str; 2] = ["c1", "c2"];
    const CONTIG_LEN: usize = 200;
    /// The tags `bam` writes; every other byte of a record stays.
    const EDITED: [[u8; 2]; 5] = [*b"MM", *b"ML", *b"XR", *b"XG", *b"XM"];

    #[derive(Debug, Clone)]
    enum Placement {
        Mapped { tid: usize, pos: usize, cigar: Vec<(char, usize)> },
        PlacedUnmapped { tid: usize, pos: usize },
        Unplaced,
    }

    #[derive(Debug, Clone)]
    struct Read {
        placement: Placement,
        /// FLAG bits besides 0x4
        flags: u16,
        seq: Vec<u8>,
        /// SAM aux fields, a stale methylation tag among them at times
        tags: Vec<&'static str>,
    }

    #[derive(Debug, Clone)]
    struct Call {
        tid: usize,
        pos: usize,
        /// `+`, the C side, or `-`, the G side
        plus: bool,
        denovo: bool,
    }

    #[derive(Debug, Clone)]
    struct Input {
        reference: Vec<Vec<u8>>,
        reads: Vec<Read>,
        /// Each read is written this many times, to cross a batch boundary
        copies: usize,
        sorted: bool,
        calls: Vec<Call>,
        /// `RASTAIR_BAM_BATCH_RECORDS` and `RASTAIR_BAM_LOOKUP_SPAN`: small,
        /// so batches and lookup runs end inside the input, among the calls
        batch_records: usize,
        lookup_span: u64,
        /// The read whose aux block gets a tag without its NUL, which `bam`
        /// must fail on, naming it
        malformed: Option<usize>,
    }

    #[hegel::composite]
    fn reads(tc: &TestCase) -> Read {
        let tid = tc.draw(gs::integers::<usize>().max_value(CONTIGS.len() - 1));
        let kind = tc.draw(gs::integers::<u8>().max_value(2));
        let (placement, query_len) = match kind {
            0 => {
                let n_ops = tc.draw(gs::integers::<usize>().min_value(1).max_value(5));
                let mut cigar =
                    vec![('M', tc.draw(gs::integers::<usize>().min_value(1).max_value(20)))];
                for _ in 1..n_ops {
                    let op = tc.draw(gs::sampled_from(['M', 'I', 'D'].as_slice()));
                    cigar.push((op, tc.draw(gs::integers::<usize>().min_value(1).max_value(6))));
                }
                if tc.draw(gs::booleans()) {
                    cigar.push(('S', tc.draw(gs::integers::<usize>().min_value(1).max_value(6))));
                }
                let ref_len: usize =
                    cigar.iter().filter(|(op, _)| "MD".contains(*op)).map(|(_, n)| n).sum();
                let query_len =
                    cigar.iter().filter(|(op, _)| "MIS".contains(*op)).map(|(_, n)| n).sum();
                let pos = tc.draw(gs::integers::<usize>().max_value(CONTIG_LEN - ref_len));
                (Placement::Mapped { tid, pos, cigar }, query_len)
            }
            1 => {
                let pos = tc.draw(gs::integers::<usize>().max_value(CONTIG_LEN - 1));
                let len = tc.draw(gs::integers::<usize>().min_value(1).max_value(20));
                (Placement::PlacedUnmapped { tid, pos }, len)
            }
            _ => (Placement::Unplaced, tc.draw(gs::integers::<usize>().min_value(1).max_value(20))),
        };
        let bits = [0x1, 0x10, 0x20, 0x40, 0x80];
        let flags = bits.iter().filter(|_| tc.draw(gs::booleans())).sum();
        let seq = tc.draw(
            gs::vecs(gs::sampled_from(b"ACGTN".as_slice())).min_size(query_len).max_size(query_len),
        );
        let optional = ["NM:i:1", "RG:Z:rg1", "XM:Z:stale", "MM:Z:C+m,0;", "XS:i:-3"];
        let tags = optional.into_iter().filter(|_| tc.draw(gs::booleans())).collect();
        Read { placement, flags, seq, tags }
    }

    #[hegel::composite]
    fn inputs(tc: &TestCase) -> Input {
        let reference = CONTIGS
            .iter()
            .map(|_| {
                tc.draw(
                    gs::vecs(gs::sampled_from(b"ACGT".as_slice()))
                        .min_size(CONTIG_LEN)
                        .max_size(CONTIG_LEN),
                )
            })
            .collect();
        let n_reads = tc.draw(gs::integers::<usize>().max_value(12));
        let mut reads: Vec<Read> = (0..n_reads).map(|_| tc.draw(reads())).collect();
        let malformed = (n_reads > 0 && tc.draw(gs::integers::<u8>().max_value(3)) == 0)
            .then(|| tc.draw(gs::integers::<usize>().max_value(n_reads - 1)));
        // Mostly an unplaced one: a batch of only those has no lookup run
        if let Some(read) = malformed.and_then(|index| reads.get_mut(index))
            && tc.draw(gs::booleans())
        {
            read.placement = Placement::Unplaced;
        }
        let copies = tc.draw(gs::sampled_from([1, 2, 300].as_slice()));
        let n_calls = tc.draw(gs::integers::<usize>().max_value(40));
        let calls = (0..n_calls)
            .map(|_| Call {
                tid: tc.draw(gs::integers::<usize>().max_value(CONTIGS.len() - 1)),
                pos: tc.draw(gs::integers::<usize>().max_value(CONTIG_LEN - 1)),
                plus: tc.draw(gs::booleans()),
                denovo: tc.draw(gs::booleans()),
            })
            .collect();
        Input {
            reference,
            reads,
            copies,
            sorted: tc.draw(gs::booleans()),
            calls,
            batch_records: tc.draw(gs::integers::<usize>().min_value(1).max_value(40)),
            lookup_span: tc.draw(gs::integers::<u64>().max_value(60)),
            malformed,
        }
    }

    impl Input {
        fn sam(&self) -> Result<String> {
            let order = if self.sorted { "coordinate" } else { "unsorted" };
            let mut sam = format!("@HD\tVN:1.6\tSO:{order}\n");
            for name in CONTIGS {
                writeln!(sam, "@SQ\tSN:{name}\tLN:{CONTIG_LEN}")?;
            }
            sam.push_str("@RG\tID:rg1\n");
            let mut reads: Vec<(usize, &Read)> = self.reads.iter().enumerate().collect();
            if self.sorted {
                reads.sort_by_key(|(_, read)| match read.placement {
                    Placement::Mapped { tid, pos, .. } | Placement::PlacedUnmapped { tid, pos } => {
                        (tid, pos)
                    }
                    Placement::Unplaced => (usize::MAX, 0),
                });
            }
            for (index, read) in reads {
                for copy in 0..self.copies {
                    let seq = String::from_utf8_lossy(&read.seq);
                    let (flags, rname, pos, cigar) = match &read.placement {
                        Placement::Mapped { tid, pos, cigar } => {
                            let mut text = String::new();
                            for (op, n) in cigar {
                                write!(text, "{n}{op}")?;
                            }
                            (read.flags, CONTIGS.get(*tid).copied().unwrap_or("*"), pos + 1, text)
                        }
                        Placement::PlacedUnmapped { tid, pos } => (
                            read.flags | 0x4,
                            CONTIGS.get(*tid).copied().unwrap_or("*"),
                            pos + 1,
                            "*".into(),
                        ),
                        Placement::Unplaced => (read.flags | 0x4, "*", 0, "*".into()),
                    };
                    let mut line = format!(
                        "r{index}_{copy}\t{flags}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{}",
                        "I".repeat(read.seq.len())
                    );
                    for tag in &read.tags {
                        line.push('\t');
                        line.push_str(tag);
                    }
                    sam.push_str(&line);
                    sam.push('\n');
                }
            }
            Ok(sam)
        }

        fn write(&self, dir: &Path) -> Result<()> {
            let mut stream =
                AlignmentStream::from_reader(std::io::Cursor::new(self.sam()?.into_bytes()))?;
            let header = BamHeader::from_template(stream.header());
            let mut writer = BamWriterBuilder::to_path(&dir.join("in.bam"), &header).build()?;
            let mut record = RawBamRecord::default();
            while stream.read_record(&mut record)? {
                if self.is_malformed(&record) {
                    let mut payload = record.as_payload().to_vec();
                    payload.extend_from_slice(b"XXZab");
                    writer.write_raw(&RawBamRecord::from_payload(payload, &header)?)?;
                } else {
                    writer.write_raw(&record)?;
                }
            }
            writer.finish()?;

            let mut fasta = String::new();
            let mut fai = String::new();
            for (name, seq) in CONTIGS.iter().zip(&self.reference) {
                let offset = fasta.len() + name.len() + 2;
                writeln!(fasta, ">{name}\n{}", String::from_utf8_lossy(seq))?;
                writeln!(fai, "{name}\t{CONTIG_LEN}\t{offset}\t{CONTIG_LEN}\t{}", CONTIG_LEN + 1)?;
            }
            std::fs::write(dir.join("ref.fa"), fasta)?;
            std::fs::write(dir.join("ref.fa.fai"), fai)?;

            let mut calls = self.calls.clone();
            calls.sort_by_key(|call| (call.tid, call.pos));
            calls.dedup_by_key(|call| (call.tid, call.pos));
            let path = dir.join("calls.bed.gz");
            let mut bed = rust_htslib::bgzf::Writer::from_path(&path)?;
            for call in &calls {
                let (strand, base) = if call.plus { ("+", "C") } else { ("-", "G") };
                let origin = if call.denovo { "NEW" } else { "REF" };
                let name = CONTIGS.get(call.tid).copied().unwrap_or("*");
                writeln!(
                    bed,
                    "{name}\t{}\t{}\t.\t0.50\t{strand}\t1\t1\t2\t0\t2\t{base}/{base}\t7\t7\t{origin}",
                    call.pos,
                    call.pos + 1
                )?;
            }
            drop(bed);
            rust_htslib::tbx::build_index(&path, rust_htslib::tbx::TabixFormat::Bed, 0)?;
            Ok(())
        }

        /// Whether `record` is a copy of the [`Self::malformed`] read
        fn is_malformed(&self, record: &RawBamRecord) -> bool {
            let name = String::from_utf8_lossy(record.qname()).into_owned();
            self.malformed.is_some_and(|index| name.starts_with(&format!("r{index}_")))
        }

        fn called(&self, tid: usize, pos: usize) -> bool {
            self.calls.iter().any(|call| call.tid == tid && call.pos == pos)
        }
    }

    fn read_bam(path: &Path) -> Result<Vec<RawBamRecord>> {
        let mut stream = AlignmentStream::from_reader(std::fs::File::open(path)?)?;
        let mut records = Vec::new();
        let mut record = RawBamRecord::default();
        while stream.read_record(&mut record)? {
            records.push(record.clone());
        }
        Ok(records)
    }

    fn unedited_tags(record: &RawBamRecord) -> Vec<u8> {
        let mut aux = AuxData::from_bytes(record.aux().as_bytes().to_vec());
        for tag in EDITED {
            aux.remove(tag);
        }
        aux.as_bytes().to_vec()
    }

    fn string_tag(record: &RawBamRecord, tag: [u8; 2]) -> Option<Vec<u8>> {
        seqair::bam::aux::find_tag(record.aux().as_bytes(), tag)
            .and_then(|value| value.as_str())
            .map(<[u8]>::to_vec)
    }

    /// The stored-SEQ indices MM marks, decoded against the read's SEQ in its
    /// original orientation as the SAM spec counts them; checks ML has one
    /// probability for each.
    fn mm_positions(record: &RawBamRecord) -> Result<Vec<usize>> {
        let Some(mm) = string_tag(record, *b"MM") else { return Ok(Vec::new()) };
        let mm = String::from_utf8(mm)?;
        let mut fields =
            mm.strip_suffix(';').ok_or_else(|| eyre!("MM {mm} lacks its ;"))?.split(',');
        let base =
            fields.next().and_then(|head| head.bytes().next()).ok_or_else(|| eyre!("MM {mm}"))?;
        let stored = seq(record);
        let reverse = record.flags().is_reverse();
        let oriented: Vec<u8> = if reverse {
            stored
                .iter()
                .rev()
                .map(|b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => b'N',
                })
                .collect()
        } else {
            stored.clone()
        };
        let of_base: Vec<usize> =
            oriented.iter().enumerate().filter(|(_, b)| **b == base).map(|(i, _)| i).collect();
        let mut next = 0;
        let mut positions = Vec::new();
        for skip in fields {
            next += skip.parse::<usize>()?;
            let at = *of_base
                .get(next)
                .ok_or_else(|| eyre!("MM {mm} skips past the {}s", base as char))?;
            positions.push(if reverse { stored.len() - 1 - at } else { at });
            next += 1;
        }
        let ml = match seqair::bam::aux::find_tag(record.aux().as_bytes(), *b"ML") {
            Some(seqair::bam::aux::AuxValue::ArrayU8(ml)) => ml.len(),
            other => bail!("ML is {other:?}"),
        };
        ensure!(ml == positions.len(), "ML has {ml} values for {} MM positions", positions.len());
        positions.sort_unstable();
        Ok(positions)
    }

    fn seq(record: &RawBamRecord) -> Vec<u8> {
        let mut seq = Vec::new();
        record.decode_seq_into(&mut seq);
        seq
    }

    /// Every input record comes out once, in order, with every byte but SEQ
    /// and the methylation tags as it was. Legacy XM marks exactly the read
    /// bases at called positions that show the strand's C (or G) or its
    /// conversion, in upper case where converted, whatever batch or lookup
    /// run the read falls in, and standard mode turns exactly those converted
    /// bases back. (`bam` once cut the genome into 10 kb segments and looked
    /// up only a segment's own calls, so a read reaching past a segment's end
    /// lost its tags there.)
    #[hegel::test(test_cases = pinned::cases(32))]
    fn bam_rewrites_every_record_and_only_its_methylation(tc: TestCase) {
        let input = tc.draw(inputs());
        let check = || -> Result<()> {
            let dir = TempDir::new()?;
            input.write(dir.path())?;
            let records = read_bam(&dir.path().join("in.bam"))?;
            let first_malformed = records.iter().find(|record| input.is_malformed(record));
            let mut outputs = Vec::new();
            for mode in ["standard", "legacy"] {
                let out = dir.path().join(format!("{mode}.bam"));
                let output = rastair()
                    .env("RASTAIR_BAM_BATCH_RECORDS", input.batch_records.to_string())
                    .env("RASTAIR_BAM_LOOKUP_SPAN", input.lookup_span.to_string())
                    .args(["bam", mode, "-@", "3", "--fasta-file"])
                    .arg(dir.path().join("ref.fa"))
                    .arg(dir.path().join("in.bam"))
                    .arg(dir.path().join("calls.bed.gz"))
                    .arg("-o")
                    .arg(&out)
                    .output()?;
                let stderr = String::from_utf8_lossy(&output.stderr);
                if let Some(record) = first_malformed {
                    let name = String::from_utf8_lossy(record.qname());
                    ensure!(!output.status.success(), "{mode} should fail on {name}");
                    let named = format!("failed to rewrite record {name}");
                    ensure!(stderr.contains(&named), "{mode} should name {name}: {stderr}");
                    ensure!(!out.exists(), "{mode} should commit no output");
                    continue;
                }
                ensure!(output.status.success(), "{mode}: {stderr}");
                let written = read_bam(&out)?;
                ensure!(
                    written.len() == records.len(),
                    "{mode}: {} of {} records",
                    written.len(),
                    records.len()
                );
                outputs.push(written);
            }
            if first_malformed.is_some() {
                return Ok(());
            }
            let (standard, legacy) = (&outputs[0], &outputs[1]);
            for ((input_record, standard), legacy) in records.iter().zip(standard).zip(legacy) {
                let name = String::from_utf8_lossy(input_record.qname()).into_owned();
                for output in [standard, legacy] {
                    ensure!(output.qname() == input_record.qname(), "order changed at {name}");
                    ensure!(output.flags() == input_record.flags(), "{name}: FLAG");
                    ensure!(
                        (output.tid(), output.pos()) == (input_record.tid(), input_record.pos()),
                        "{name}: position"
                    );
                    ensure!(output.mapq() == input_record.mapq(), "{name}: MAPQ");
                    ensure!(
                        output.stored_cigar_bytes() == input_record.stored_cigar_bytes(),
                        "{name}: CIGAR"
                    );
                    ensure!(output.qual() == input_record.qual(), "{name}: QUAL");
                    ensure!(output.mate_fields() == input_record.mate_fields(), "{name}: mate");
                    ensure!(unedited_tags(output) == unedited_tags(input_record), "{name}: tags");
                }
                let original = seq(input_record);
                ensure!(seq(legacy) == original, "{name}: legacy changed SEQ");

                let strand = Strand::from_flags(input_record.flags());
                let mut expected = vec![b'.'; original.len()];
                if let (Some(strand), Some(tid), Some(pos)) =
                    (strand, input_record.tid(), input_record.pos())
                {
                    let (unmod, conv) = match strand {
                        Strand::OT => (b'C', b'T'),
                        Strand::OB => (b'G', b'A'),
                    };
                    let cigar: Vec<CigarOp> = input_record.cigar().collect();
                    for aligned in AlignedPairs::new(pos, &cigar).matches_only() {
                        let qpos = aligned.qpos.as_usize();
                        if !input.called(tid.as_usize(), aligned.rpos.as_usize()) {
                            continue;
                        }
                        let (Some(&base), Some(mark)) =
                            (original.get(qpos), expected.get_mut(qpos))
                        else {
                            continue;
                        };
                        if base == unmod {
                            *mark = b'z';
                        } else if base == conv {
                            *mark = b'Z';
                        }
                    }
                }
                let xm = string_tag(legacy, *b"XM").ok_or_else(|| eyre!("{name}: no XM"))?;
                // De-novo CpGs take their letter (z, x or h) from the context
                let class = |c: &u8| match c {
                    b'.' => b'.',
                    c if c.is_ascii_uppercase() => b'Z',
                    _ => b'z',
                };
                ensure!(
                    xm.iter().map(class).collect::<Vec<_>>() == expected,
                    "{name}: XM {} where {} was expected",
                    String::from_utf8_lossy(&xm),
                    String::from_utf8_lossy(&expected)
                );

                let rewritten = seq(standard);
                ensure!(rewritten.len() == original.len(), "{name}: SEQ length");
                for (i, ((&before, &after), &mark)) in
                    original.iter().zip(&rewritten).zip(&expected).enumerate()
                {
                    let converted_back = matches!((before, after), (b'T', b'C') | (b'A', b'G'));
                    ensure!(
                        before == after || (mark == b'Z' && converted_back),
                        "{name}: SEQ {i} {} became {}",
                        before as char,
                        after as char
                    );
                    ensure!(mark != b'Z' || converted_back, "{name}: SEQ {i} not converted back");
                }
                ensure!(
                    string_tag(standard, *b"MM").is_some() == expected.contains(&b'Z'),
                    "{name}: MM present without methylation, or missing"
                );
                let converted: Vec<usize> = expected
                    .iter()
                    .enumerate()
                    .filter(|(_, m)| **m == b'Z')
                    .map(|(i, _)| i)
                    .collect();
                ensure!(mm_positions(standard)? == converted, "{name}: MM/ML marks other bases");

                let first = input_record.flags().is_first_in_template();
                let xr: &[u8] = if first { b"CT" } else { b"GA" };
                let xg: &[u8] = if strand == Some(Strand::OB) { b"GA" } else { b"CT" };
                ensure!(string_tag(legacy, *b"XR").as_deref() == Some(xr), "{name}: XR");
                ensure!(string_tag(legacy, *b"XG").as_deref() == Some(xg), "{name}: XG");
            }
            Ok(())
        };
        let result = check();
        assert!(result.is_ok(), "{result:?}");
    }
}

/// A coordinate-sorted BAM of three reads in `c1`; with `invalid` the second
/// has a CIGAR longer than its sequence. The records are as long either way,
/// so an index built on the valid file also fits the invalid one.
fn write_three_reads_bam(dir: &std::path::Path, invalid: bool) -> Result<std::path::PathBuf> {
    use seqair::bam::RawBamRecord;
    use seqair::reader::AlignmentStream;
    use std::io::Write as _;

    let sam = "@SQ\tSN:c1\tLN:100\n\
               a\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
               b\t0\tc1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
               c\t0\tc1\t3\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    let mut stream = AlignmentStream::from_reader(std::io::Cursor::new(sam.as_bytes().to_vec()))?;
    let header = stream.header().clone();
    let mut payloads = Vec::new();
    let mut record = RawBamRecord::default();
    while stream.read_record(&mut record)? {
        payloads.push(record.as_payload().to_vec());
    }
    // The second record's CIGAR becomes 5M: the first op follows the
    // 32-byte fixed fields and the NUL-terminated name.
    let cigar_at = 32 + 2;
    let op = payloads
        .get_mut(1)
        .and_then(|payload| payload.get_mut(cigar_at..cigar_at + 4))
        .ok_or_else(|| eyre!("no second record"))?;
    op.copy_from_slice(&(if invalid { 5_u32 } else { 4 } << 4).to_le_bytes());

    let input = dir.join("in.bam");
    let mut bam_bytes = b"BAM\x01".to_vec();
    let text = header.header_text().as_bytes();
    bam_bytes.extend_from_slice(&u32::try_from(text.len())?.to_le_bytes());
    bam_bytes.extend_from_slice(text);
    bam_bytes.extend_from_slice(&1_u32.to_le_bytes());
    bam_bytes.extend_from_slice(&3_u32.to_le_bytes());
    bam_bytes.extend_from_slice(b"c1\0");
    bam_bytes.extend_from_slice(&100_u32.to_le_bytes());
    for payload in &payloads {
        bam_bytes.extend_from_slice(&u32::try_from(payload.len())?.to_le_bytes());
        bam_bytes.extend_from_slice(payload);
    }
    let mut writer = rust_htslib::bgzf::Writer::from_path(&input)?;
    writer.write_all(&bam_bytes)?;
    drop(writer);

    Ok(input)
}

/// A record that is not valid BAM (here its CIGAR is longer than its SEQ)
/// stops the run with an error naming it, and nothing is committed: `bam`
/// rewrites losslessly or not at all.
#[test]
fn an_invalid_record_fails_the_run() -> Result<()> {
    use std::io::Write as _;

    let temp_dir = TempDir::new()?;
    let input = write_three_reads_bam(temp_dir.path(), true)?;

    let fasta = write_fasta_with_contigs(temp_dir.path(), &["c1"])?;
    let calls = temp_dir.path().join("calls.bed.gz");
    let mut bed = rust_htslib::bgzf::Writer::from_path(&calls)?;
    writeln!(bed, "c1\t1\t2\t.\t0.50\t+\t1\t1\t2\t0\t2\tC/C\t7\t7\tREF")?;
    drop(bed);
    rust_htslib::tbx::build_index(&calls, rust_htslib::tbx::TabixFormat::Bed, 0)?;

    for mode in ["standard", "legacy"] {
        let out = temp_dir.path().join(format!("{mode}.bam"));
        let output = rastair()
            .args(["bam", mode, "--fasta-file"])
            .arg(&fasta)
            .arg(&input)
            .arg(&calls)
            .arg("-o")
            .arg(&out)
            .output()?;
        let stderr = String::from_utf8_lossy(&output.stderr);
        ensure!(!output.status.success(), "{mode} should fail: {stderr}");
        ensure!(stderr.contains("BAM record 2 is invalid"), "{mode}: {stderr}");
        ensure!(!out.exists(), "{mode} should commit no output");
        let partial = temp_dir.path().join(format!("{mode}.bam.partial"));
        let written = bam::Reader::from_path(&partial)?.records().count();
        ensure!(written == 1, "{mode}: the records before it are left, got {written}");
    }
    Ok(())
}

/// Through the index, an invalid record fails the run with its BGZF virtual
/// offset, which `samtools view -X` and a hex dump can find.
#[test]
fn an_invalid_record_under_region_reports_its_offset() -> Result<()> {
    use std::io::Write as _;

    let temp_dir = TempDir::new()?;
    let input = write_three_reads_bam(temp_dir.path(), false)?;
    bam::index::build(&input, None, bam::index::Type::Bai, 1)?;
    write_three_reads_bam(temp_dir.path(), true)?;
    let fasta = write_fasta_with_contigs(temp_dir.path(), &["c1"])?;
    let calls = temp_dir.path().join("calls.bed.gz");
    let mut bed = rust_htslib::bgzf::Writer::from_path(&calls)?;
    writeln!(bed, "c1\t1\t2\t.\t0.50\t+\t1\t1\t2\t0\t2\tC/C\t7\t7\tREF")?;
    drop(bed);
    rust_htslib::tbx::build_index(&calls, rust_htslib::tbx::TabixFormat::Bed, 0)?;

    let output = rastair()
        .args(["bam", "legacy", "--fasta-file"])
        .arg(&fasta)
        .arg(&input)
        .arg(&calls)
        .args(["--region", "c1", "-o"])
        .arg(temp_dir.path().join("out.bam"))
        .output()?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    ensure!(!output.status.success(), "should fail: {stderr}");
    ensure!(stderr.contains("virtual offset 0x"), "{stderr}");
    Ok(())
}

/// `--region` with several regions keeps each overlapping read once, also one
/// reaching into two of them, as `samtools view -M` does.
#[test]
fn a_read_in_two_regions_is_written_once() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    rastair()
        .args(CALL_TEST_BAM)
        .args(["--region=chr19:6103000-6103300", "--cpgs-only", "--bed"])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    // 1-based closed regions; reads are about 80 bases, so many span the gaps
    let regions = [(6103000, 6103050), (6103060, 6103100), (6103080, 6103300)];
    let mut expected = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut bam = bam::IndexedReader::from_path("tests/data/test.bam")?;
    for (start, end) in regions {
        bam.fetch(("chr19", start - 1, end))?;
        for record in bam.records() {
            let record = record?;
            let key = (record.qname().to_owned(), record.flags());
            if seen.insert(key.clone()) {
                expected.push(key);
            }
        }
    }
    expected.sort_by_key(|(name, flags)| (name.clone(), *flags));

    let output_bam = temp_dir.path().join("out.bam");
    let region = regions.map(|(start, end)| format!("chr19:{start}-{end}")).join(" ");
    rastair()
        .args(["bam", "legacy", "--fasta-file=tests/data/test.fasta.gz", "--region", &region])
        .arg("tests/data/test.bam")
        .arg(&calls_bed)
        .arg("-o")
        .arg(&output_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite BAM")?;
    let mut written: Vec<(Vec<u8>, u16)> = bam::Reader::from_path(&output_bam)?
        .records()
        .map(|record| {
            Ok::<_, color_eyre::Report>(record.map(|r| (r.qname().to_owned(), r.flags()))?)
        })
        .collect::<Result<_>>()?;
    let spanning = written.len();
    written.sort();
    ensure!(written.len() > 10, "too few reads to test");
    ensure!(written == expected, "{spanning} reads written, {} expected", expected.len());
    Ok(())
}

/// `bam` on its own output adds a second `@PG` with an ID of its own, whose
/// `PP` is the first.
#[test]
fn bam_twice_gives_unique_pg_ids() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    rastair()
        .args(CALL_TEST_BAM)
        .args(["--region=chr19:6103000-6103300", "--cpgs-only", "--bed"])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;
    let mut input = Path::new("tests/data/test.bam").to_path_buf();
    for run in ["first", "second"] {
        let output = temp_dir.path().join(format!("{run}.bam"));
        rastair()
            .args(["bam", "legacy", "--fasta-file=tests/data/test.fasta.gz"])
            .arg(&input)
            .arg(&calls_bed)
            .arg("-o")
            .arg(&output)
            .silent()
            .succeeds()
            .wrap_err_with(|| format!("{run} run failed"))?;
        input = output;
    }
    let header = seqair::bam::BamHeader::from_bam_path(&input)?;
    let pgs: Vec<Vec<&str>> = header
        .header_text()
        .lines()
        .filter(|line| line.starts_with("@PG"))
        .map(|line| line.split('\t').collect())
        .collect();
    let field =
        |pg: &[&str], key: &str| pg.iter().find_map(|f| f.strip_prefix(key)).map(str::to_owned);
    let ids: Vec<String> = pgs.iter().filter_map(|pg| field(pg, "ID:")).collect();
    let unique: std::collections::HashSet<&String> = ids.iter().collect();
    ensure!(unique.len() == ids.len(), "duplicate @PG IDs: {ids:?}");
    let last = pgs.last().ok_or_else(|| eyre!("no @PG"))?;
    ensure!(field(last, "ID:").as_deref() == Some("rastair.rewrite_bam.1"), "{last:?}");
    ensure!(field(last, "PP:").as_deref() == Some("rastair.rewrite_bam"), "{last:?}");
    Ok(())
}

/// A name-sorted CRAM (`SO:queryname`), as aligners and `samtools sort -n`
/// write it, is rewritten whole, in its own order.
#[test]
fn a_name_sorted_cram_keeps_its_order() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let cram = temp_dir.path().join("by_name.cram");
    let sorted = Command::new("samtools")
        .args(["sort", "-n", "-O", "cram", "--reference", "tests/data/test.fasta.gz", "-o"])
        .arg(&cram)
        .arg("tests/data/test.bam")
        .output()
        .wrap_err("Failed to run samtools")?;
    ensure!(sorted.status.success(), "samtools: {}", String::from_utf8_lossy(&sorted.stderr));

    let calls_bed = temp_dir.path().join("calls.bed.gz");
    rastair()
        .args(CALL_TEST_BAM)
        .args(["--region=chr19:6103000-6103300", "--cpgs-only", "--bed"])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;
    let output_bam = temp_dir.path().join("out.bam");
    rastair()
        .args(["bam", "legacy", "--fasta-file=tests/data/test.fasta.gz"])
        .arg(&cram)
        .arg(&calls_bed)
        .arg("-o")
        .arg(&output_bam)
        .silent()
        .succeeds()
        .wrap_err("Failed to rewrite the name-sorted CRAM")?;

    let names = |path: &Path| -> Result<Vec<(Vec<u8>, u16)>> {
        let mut reader = bam::Reader::from_path(path)?;
        reader.set_reference("tests/data/test.fasta.gz")?;
        reader
            .records()
            .map(|record| Ok(record.map(|r| (r.qname().to_owned(), r.flags()))?))
            .collect()
    };
    let input = names(&cram)?;
    ensure!(input.len() == 29_358, "the CRAM has {} records", input.len());
    ensure!(names(&output_bam)? == input, "every record, in the CRAM's order");
    Ok(())
}
