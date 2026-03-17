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
            "--no-ml",
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
    ensure!(tabix_output.status.success(), "tabix failed: {}", String::from_utf8_lossy(&tabix_output.stderr));

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
            "--no-ml",
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
    ensure!(tabix_output.status.success(), "tabix failed: {}", String::from_utf8_lossy(&tabix_output.stderr));

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

/// Rewritten BAM must be sorted by position and contain no duplicate reads,
/// even when `--segment-max-length` forces processing across multiple segments.
/// Reads overlapping a segment boundary used to be emitted by both the segment
/// they start in and the next one, causing duplicates and sort-order violations
/// that made `samtools index` fail.
#[test]
fn bam_rewrite_is_sorted_with_small_segments() -> Result<()> {
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
            "--no-ml",
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
    ensure!(tabix_output.status.success(), "tabix failed: {}", String::from_utf8_lossy(&tabix_output.stderr));

    // Use a very small segment length to force many segment boundaries,
    // maximising the chance that reads straddle them.
    rastair()
        .args([
            "bam",
            "standard",
            "--fasta-file=tests/data/test.fasta.gz",
            "--region=chr19:6103000-6106000",
            "--segment-max-length=200",
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
            "--no-ml",
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
    ensure!(tabix_output.status.success(), "tabix failed: {}", String::from_utf8_lossy(&tabix_output.stderr));

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

/// Regression test for mixed-methylation sites in legacy XM output.
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
            "--no-ml",
            &format!("--region={region}"),
            "--cpgs-only",
            "-o",
        ])
        .arg(&calls_bed)
        .silent()
        .succeeds()
        .wrap_err("Failed to generate calls")?;

    // `rastair call` already writes the tabix index for BED output.
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
