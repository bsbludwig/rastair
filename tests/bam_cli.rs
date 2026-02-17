mod utils;
use rust_htslib::bam::Read as BamRead;
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

    std::process::Command::new("tabix")
        .args(["-p", "bed"])
        .arg(&calls_bed)
        .output()
        .wrap_err("Failed to index calls file")?;

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
    std::process::Command::new("tabix")
        .args(["-p", "bed"])
        .arg(&calls_bed)
        .output()
        .wrap_err("Failed to index calls file")?;

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
