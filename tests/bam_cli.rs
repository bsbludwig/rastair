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

/// Rewritten BAM must be sorted by position and contain no duplicate reads,
/// even when `--segment-max-length` forces processing across multiple segments.
/// A read overlapping a segment boundary is emitted by the segment it starts
/// in only; emitting it from the next one too would duplicate it and break the
/// sort order, and `samtools index` would fail.
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

    // Use the smallest segment length to force several segment boundaries,
    // maximising the chance that reads straddle them.
    rastair()
        .args([
            "bam",
            "standard",
            "--fasta-file=tests/data/test.fasta.gz",
            "--region=chr19:6103000-6106000",
            "--segment-max-length=1000",
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
