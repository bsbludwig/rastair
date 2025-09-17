mod utils;
use std::fs;

use utils::*;

#[test]
fn simple_per_read_call() -> Result<()> {
    apply_common_filters!();

    assert_cmd_snapshot!(rastair().args([
        "per-read",
        "--fasta-file=tests/data/test.fasta.gz",
        "tests/data/test.bam",
        "--region=chr19:6105900-6105950",
    ]), @r"
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
    [TIME] INFO rastair: Calling reads finished [DURATION]
    ");

    Ok(())
}

#[test]
fn enhance_with_calls_bed_file_to_add_denovo_counts() -> Result<()> {
    apply_common_filters!();

    let temp_dir = TempDir::new()?;
    let calls_bed = temp_dir.path().join("calls.bed.gz");
    rastair()
        .args([
            "call",
            "--fasta-file=tests/data/test.fasta.gz",
            "tests/data/test.bam",
            "-c",
            "--region=chr19:6105900-6105950",
            "--thresholds",
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

    fs::remove_file(temp_dir.path().join("test.bed.gz.tbi"))
        .wrap_err("Failed to remove existing tabix index")?;

    assert_cmd_snapshot!(Command::new("tabix")
            .args(["-p", "bed"])
            .arg(&bed_file)
            .current_dir(temp_dir.path()), @r"
    success: true
    exit_code: 0
    ----- stdout -----

    ----- stderr -----
    ");

    Ok(())
}
