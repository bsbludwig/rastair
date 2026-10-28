mod utils;
use std::{fs, path::PathBuf};

use utils::*;

#[test]
fn simple_per_read_call() -> Result<()> {
    apply_common_filters!();

    assert_cmd_snapshot!(rastair().args([
        "per-read",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105900-6105950",
    ]), @r#"
    success: true
    exit_code: 0
    ----- stdout -----
    #chr	start	end	read_id	mapq	orientation	insert_size	read_length	flag	num_cpg	num_mod	mod_cpgs	unmod_cpgs	snp_cpgs	mod_denovos	unmod_denovos
    chr19	6105906	6105986	NB502094:69:HN2H2BGX5:2:12112:7718:7791	60	-	233	80	147	3	3	33,46,77				
    chr19	6105906	6105984	NB502094:69:HN2H2BGX5:2:12307:1537:12140	60	-	164	78	147	3	2	33,46	77			
    chr19	6105906	6105984	NB502094:69:HN2H2BGX5:4:22407:12172:14935	60	-	132	78	83	2	2	34,47				
    chr19	6105908	6105985	NB502094:69:HN2H2BGX5:2:12310:16034:14314	60	+	196	77	99	3	3	31,44,75				
    chr19	6105914	6105993	NB502094:69:HN2H2BGX5:4:21611:24189:11108	60	-	255	79	83	3	2	39,70	26			
    chr19	6105915	6105995	NB502094:69:HN2H2BGX5:1:12305:15843:12087	60	+	225	80	163	3	1	38	25,69			
    chr19	6105921	6106000	NB502094:69:HN2H2BGX5:4:13611:7898:19397	60	-	160	79	83	4	3	32,63,75	19			
    chr19	6105925	6106005	NB502094:69:HN2H2BGX5:3:11507:14836:15120	60	-	236	80	147	4	4	14,27,58,70				
    chr19	6105926	6106006	NB502094:69:HN2H2BGX5:4:13411:10524:1389	60	+	282	80	99	4	4	13,26,57,69				
    chr19	6105926	6106006	NB502094:69:HN2H2BGX5:3:23508:14923:3979	60	-	160	80	83	4	4	14,27,58,70				
    chr19	6105935	6106015	NB502094:69:HN2H2BGX5:2:13204:9897:11002	60	+	316	80	163	4	2	18,61	5,49			
    chr19	6105940	6106020	NB502094:69:HN2H2BGX5:3:21401:25012:7155	60	-	249	80	147	3	1	12	43,55			
    chr19	6105941	6106020	NB502094:69:HN2H2BGX5:4:21411:10562:9167	60	+	190	79	163	3	3	12,43,55				
    chr19	6105947	6106026	NB502094:69:HN2H2BGX5:1:12302:4502:11914	60	-	263	79	83	3	2	6,37	49			
    chr19	6105947	6106027	NB502094:69:HN2H2BGX5:2:12207:14946:2414	60	-	241	80	83	3	3	6,37,49				

    ----- stderr -----
    [TIME] INFO rastair::runtime::partial_output: Wrote BED output file="-"
    [TIME] INFO rastair: Calling reads finished [DURATION]
    "#);

    Ok(())
}

#[test]
fn per_read_with_guess_read_orientation() -> Result<()> {
    apply_common_filters!();

    assert_cmd_snapshot!(rastair()
        .args([
            "per-read",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6105900-6105950",
            "--guess-read-orientation",
        ]), @r#"
    success: true
    exit_code: 0
    ----- stdout -----
    #chr	start	end	read_id	mapq	orientation	insert_size	read_length	flag	num_cpg	num_mod	mod_cpgs	unmod_cpgs	snp_cpgs	mod_denovos	unmod_denovos
    chr19	6105906	6105986	NB502094:69:HN2H2BGX5:2:12112:7718:7791	60	-	233	80	147	3	3	33,46,77				
    chr19	6105906	6105984	NB502094:69:HN2H2BGX5:2:12307:1537:12140	60	-	164	78	147	3	2	33,46	77			
    chr19	6105906	6105984	NB502094:69:HN2H2BGX5:4:22407:12172:14935	60	-	132	78	83	2	2	34,47				
    chr19	6105908	6105985	NB502094:69:HN2H2BGX5:2:12310:16034:14314	60	+	196	77	99	3	3	31,44,75				
    chr19	6105914	6105993	NB502094:69:HN2H2BGX5:4:21611:24189:11108	60	-	255	79	83	3	2	39,70	26			
    chr19	6105915	6105995	NB502094:69:HN2H2BGX5:1:12305:15843:12087	60	+	225	80	163	3	1	38	25,69			
    chr19	6105921	6106000	NB502094:69:HN2H2BGX5:4:13611:7898:19397	60	-	160	79	83	4	3	32,63,75	19			
    chr19	6105925	6106005	NB502094:69:HN2H2BGX5:3:11507:14836:15120	60	-	236	80	147	4	4	14,27,58,70				
    chr19	6105926	6106006	NB502094:69:HN2H2BGX5:4:13411:10524:1389	60	+	282	80	99	4	4	13,26,57,69				
    chr19	6105926	6106006	NB502094:69:HN2H2BGX5:3:23508:14923:3979	60	-	160	80	83	4	4	14,27,58,70				
    chr19	6105935	6106015	NB502094:69:HN2H2BGX5:2:13204:9897:11002	60	+	316	80	163	4	2	18,61	5,49			
    chr19	6105940	6106020	NB502094:69:HN2H2BGX5:3:21401:25012:7155	60	-	249	80	147	3	1	12	43,55			
    chr19	6105941	6106020	NB502094:69:HN2H2BGX5:4:21411:10562:9167	60	+	190	79	163	3	3	12,43,55				
    chr19	6105947	6106026	NB502094:69:HN2H2BGX5:1:12302:4502:11914	60	-	263	79	83	3	2	6,37	49			
    chr19	6105947	6106027	NB502094:69:HN2H2BGX5:2:12207:14946:2414	60	-	241	80	83	3	3	6,37,49				

    ----- stderr -----
    [TIME] INFO rastair::runtime::partial_output: Wrote BED output file="-"
    [TIME] INFO rastair: Calling reads finished [DURATION]
    "#);

    Ok(())
}

#[test]
fn enhance_with_calls_bed_file_to_add_denovo_counts() -> Result<()> {
    apply_common_filters!();

    // Create explicit temp dir because its name is gonna be part of the snapshot
    let tmp_dir = PathBuf::from("tmp/tests");
    fs::create_dir_all(&tmp_dir)?;
    let calls_bed = tmp_dir.join("calls.bed.gz");
    match fs::remove_file(&calls_bed) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e.into()),
    }

    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "-c",
            "--region=chr19:6105900-6105950",
            "--bed",
        ])
        .arg(&calls_bed)
        .succeeds()
        .wrap_err("Failed to run call")?;

    assert_cmd_snapshot!(
        rastair()
            .args([
                "per-read",
                "--fasta-file=tests/data/test.fasta.gz",
                "tests/data/test.bam",
                "--region=chr19:6105900-6105950",
                "--calls",
            ])
            .arg(&calls_bed)
    );

    fs::remove_file(calls_bed).wrap_err("Failed to remove calls.bed.gz")?;

    Ok(())
}

#[test]
fn fails_when_given_a_calls_file_as_stdin() -> Result<()> {
    apply_common_filters!();

    assert_cmd_snapshot!(rastair()
        .args([
            "per-read",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6105900-6105950",
            "--calls=-",
        ])
        .pass_stdin("#chr"),
        @r#"
    success: false
    exit_code: 1
    ----- stdout -----

    ----- stderr -----
    Error: 
       0: Failed to read calls from BED file "-"
       1: Failed to open tabix file
       2: No such file or directory (os error 2)

    Note: Rastair can only read calls from local .bed.gz files as it also needs to load the tabix index
    "#
    );

    Ok(())
}

#[test]
fn fails_when_given_a_calls_file_with_no_tbi() -> Result<()> {
    apply_common_filters!();

    assert_cmd_snapshot!(rastair()
        .args([
            "per-read",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "--region=chr19:6105900-6105950",
            "--calls=tests/data/rastair1.bed", // this file has no .tbi
        ])
        .pass_stdin("#chr"),
        @r#"
    success: false
    exit_code: 1
    ----- stdout -----

    ----- stderr -----
    Error: 
       0: Failed to read calls from BED file "tests/data/rastair1.bed"
       1: Failed to open tabix file
       2: No such file or directory (os error 2)

    Suggestion: You can compress and index the calls file with:

    bgzip tests/data/rastair1.bed
    tabix tests/data/rastair1.bed.gz
    "#
    );

    Ok(())
}

#[test]
fn can_tabix_files() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let bed_file = temp_dir.path().join("test.bed.gz");

    rastair()
        .args(["per-read", "--fasta-file=tests/data/test.fasta.gz", "tests/data/test.bam", "--bed"])
        .arg(&bed_file)
        .succeeds()
        .wrap_err("Failed to run per-read")?;

    assert_cmd_snapshot!(Command::new("tabix")
        .args(["-f", "-p", "bed"])
        .arg(&bed_file)
        .current_dir(temp_dir.path()),
        @r"
        success: true
        exit_code: 0
        ----- stdout -----

        ----- stderr -----
        "
    );

    Ok(())
}

/// Reads starting exactly at a segment boundary are reported once: the BAM
/// fetch is one base wider than the segment, so both bounds are checked.
/// Small segments put many reads on a boundary, so the output must still match
/// a run that processes everything in one segment.
#[test]
fn reads_at_segment_boundaries_are_reported_once() -> Result<()> {
    let run = |segment_max_length: &str| -> Result<Vec<String>> {
        let output = rastair()
            .args([
                "per-read",
                "--fasta-file=tests/data/test.fasta.gz",
                "tests/data/test.bam",
                "--all-reads",
                "--region=chr19:6105000-6205000",
                "--segment-max-length",
                segment_max_length,
            ])
            .output()
            .wrap_err("Failed to run per-read")?;
        output.clone().succeeds()?;
        Ok(output.stdout().lines().map(str::to_owned).collect())
    };

    let chunked = run("1000")?;
    let single = run("1000000")?;

    let duplicates: Vec<_> = {
        let mut seen = BTreeSet::new();
        chunked.iter().filter(|line| !seen.insert((*line).clone())).collect()
    };
    ensure!(duplicates.is_empty(), "duplicated read rows: {duplicates:#?}");
    ensure!(chunked.len() > 100, "expected a non-trivial number of reads, got {}", chunked.len());
    ensure!(chunked == single, "segmentation changed the reported reads");

    Ok(())
}

/// The de-novo columns go by each read's base, like `bam`'s tags: a read on
/// the strand that converts the CpG side lists it as methylated where it
/// shows the converted base and unmethylated where it shows the CpG base; a
/// read on the other strand lists it in neither. Here β is 5/8, but three
/// reads on each strand are unmethylated.
#[test]
fn denovo_columns_follow_each_reads_base() -> Result<()> {
    let stdout = per_read_over_a_denovo_cpg(&["--all-reads"])?;
    let (methylated, unmethylated) = (listed(&stdout, 14), listed(&stdout, 15));
    ensure!(
        (methylated, unmethylated) == (10, 6),
        "mod_denovos {methylated}, unmod_denovos {unmethylated}:\n{stdout}"
    );
    Ok(())
}

/// A read whose only CpG is de-novo is reported without `--all-reads`, and
/// counts it in `num_cpg`.
#[test]
fn a_read_with_only_a_denovo_cpg_is_reported() -> Result<()> {
    let stdout = per_read_over_a_denovo_cpg(&[])?;
    let (methylated, unmethylated) = (listed(&stdout, 14), listed(&stdout, 15));
    ensure!(
        (methylated, unmethylated) == (10, 6),
        "mod_denovos {methylated}, unmod_denovos {unmethylated}:\n{stdout}"
    );
    let rows: Vec<_> = stdout.lines().filter(|line| !line.starts_with('#')).collect();
    ensure!(
        rows.len() == 16 && rows.iter().all(|row| row.split('\t').nth(9) == Some("1")),
        "expected 16 reads with num_cpg 1:\n{stdout}"
    );
    let num_mod = |value| rows.iter().filter(|row| row.split('\t').nth(10) == Some(value)).count();
    ensure!(
        (num_mod("1"), num_mod("0")) == (10, 6),
        "expected num_mod 1 on the 10 methylated reads:\n{stdout}"
    );
    Ok(())
}

/// The number of positions all rows list in `column`.
fn listed(stdout: &str, column: usize) -> usize {
    stdout
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.split('\t').nth(column))
        .map(|positions| positions.split(',').filter(|p| !p.is_empty()).count())
        .sum()
}

/// `per-read` with `args` over a homozygous de-novo CpG with β 5/8, given
/// the calls `call --bed` makes there.
fn per_read_over_a_denovo_cpg(args: &[&str]) -> Result<String> {
    use utils::scenario::{Hap, Scenario, Strand};
    const DENOVO: usize = 600;

    let scenario = Scenario::new(1_000)
        .reference_at(DENOVO, b"TG")
        .hom(DENOVO, b'C')
        .methylated(Hap::One)
        .pairs(Hap::One, Strand::Top, DENOVO - 50, 5)
        .pairs(Hap::Two, Strand::Top, DENOVO - 50, 3)
        .pairs(Hap::One, Strand::Bottom, DENOVO - 50, 5)
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
    let row = |pos: usize| {
        let pos = pos.to_string();
        text.lines().find(|line| line.split('\t').nth(1) == Some(pos.as_str()))
    };
    for pos in [DENOVO, DENOVO + 1] {
        let fields: Vec<_> = row(pos).unwrap_or_default().split('\t').collect();
        ensure!(
            fields.get(4) == Some(&"0.62") && fields.get(14) == Some(&"NEW"),
            "no de-novo CpG with β 0.62 at {pos}:\n{text}"
        );
    }

    let output = rastair()
        .arg("per-read")
        .arg("--fasta-file")
        .arg(&fasta)
        .arg(&reads)
        .args(args)
        .arg("--calls")
        .arg(&calls)
        .output()?;
    ensure!(output.status.success(), "per-read failed: {}", output.stderr());
    Ok(String::from_utf8(output.stdout)?)
}
