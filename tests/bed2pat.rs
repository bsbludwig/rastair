use assert_cmd::prelude::*;
// Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs
use tempfile::NamedTempFile;

fn stage_read_bed(region: &str) -> Result<NamedTempFile, Box<dyn std::error::Error>>
{
    // Run rastair per-read on some region
    let mut cmd = Command::cargo_bin("rastair")?;
    let file = NamedTempFile::new()?;
    let write_handle = file.reopen()?;

    cmd.arg("per-read");
    cmd.args(["--fasta-file", "test_data/test.fasta.gz"]);
    if region.len() > 0 {
        cmd.args(["-l", region]);
    }
    cmd.arg("test_data/test.bam")
       .stdout(write_handle)
       .status()
       .expect("Failed to create per-read file needed for tests");

    Ok(file)
}

#[test]
fn missing_bed() -> Result<(), Box<dyn std::error::Error>>
{
    // stage input file
    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("bed2pat");
    cmd.args(["-r", "test_data/test.fasta.gz"]);
    cmd.arg("/path/to/nonexistent/file.bed");
    cmd.assert()
       .failure()
       .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn missing_fasta() -> Result<(), Box<dyn std::error::Error>>
{
    let file = stage_read_bed("bacteriophage_lambda_CpG:1-1000")?;

    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("bed2pat");
    cmd.args(["-r", "test_data/test_which_doesnt_exist.fasta"]);
    cmd.arg(file.path());
    cmd.assert()
       .failure()
       .stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn can_create_pat() -> Result<(), Box<dyn std::error::Error>>
{
    let file = stage_read_bed("bacteriophage_lambda_CpG:1-200")?;

    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("bed2pat");
    cmd.args(["-r", "test_data/test.fasta.gz"]);
    cmd.arg(file.path());
    cmd.assert().success();

    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let line_count = output_str.lines().count();
    assert!(line_count > 0);

    // there should be 24 reads reported in total
    let mut read_count: i32 = output_str.lines()
                                        .map(|line| {
                                            let elems: Vec<&str> = line.split("\t").collect();
                                            // this is safe because I filtered for len>9 above
                                            if elems.len() >= 4 {
                                                elems[3].parse::<i32>().unwrap_or_default()
                                            } else {
                                                0
                                            }
                                        })
                                        .sum();
    // There are two reads in this region that contain an indel in the
    // CpG, but their mate is outside the 1-200 window, so only 22 out of
    // 24 read pairs are reported
    assert_eq!(read_count, 22);

    // There should be 6 reads that had TTCCCCCCC - that's because 3 had
    // a pair with a deletion on the CpG, and 3 had a pair that's outside
    // the 1-200 window so were not in the bam file
    let ttccc_line = output_str.lines()
                               .find(|line| {
                                   let elems: Vec<&str> = line.split("\t").collect();
                                   // this is safe because I filtered for len>9 above
                                   if elems.len() >= 4 && elems[2] == "TTCCCCCCC" {
                                       true
                                   } else {
                                       false
                                   }
                               })
                               .unwrap_or_default();
    read_count = ttccc_line.split("\t")
                           .nth(3)
                           .unwrap_or("0")
                           .parse::<i32>()
                           .unwrap_or_default();
    assert_eq!(read_count, 6);
    assert_eq!(ttccc_line.split("\t").nth(1).unwrap_or_default(), "1");

    // there should be a read starting at position 12
    let p12_line = output_str.lines().find(|line| {
                                         let elems: Vec<&str> = line.split("\t").collect();
                                         // this is safe because I filtered for len>9 above
                                         if elems.len() >= 4 && elems[1] == "12" {
                                             true
                                         } else {
                                             false
                                         }
                                     });
    assert!(p12_line.is_some(),
            "Could not find read starting at position 12");
    Ok(())
}

#[test]
fn test_in_genomic_region() -> Result<(), Box<dyn std::error::Error>>
{
    let file = stage_read_bed("chr19:6109125-6109567")?;

    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("bed2pat");
    cmd.args(["-r", "test_data/test.fasta.gz"]);
    cmd.arg(file.path());
    cmd.assert().success();

    let output = cmd.output().unwrap();
    let output_str = String::from_utf8_lossy(&output.stdout);
    let line_count = output_str.lines().count();
    assert!(line_count > 0);

    // there should be 30 reads reported in total
    let read_count: i32 = output_str.lines()
                                    .map(|line| {
                                        let elems: Vec<&str> = line.split("\t").collect();
                                        // this is safe because I filtered for len>9 above
                                        if elems.len() >= 4 {
                                            elems[3].parse::<i32>().unwrap_or_default()
                                        } else {
                                            0
                                        }
                                    })
                                    .sum();
    assert_eq!(read_count, 30);

    // just compare the whole thing:
    assert_eq!(output_str.to_string(), "chr19\t42369\tC\t8\nchr19\t42369\tT\t4\nchr19\t42370\tC\t8\nchr19\t42370\tT\t7\nchr19\t42371\tC\t2\nchr19\t42371\tT\t1\n");
    Ok(())
}

#[test]
fn test_with_clipping() -> Result<(), Box<dyn std::error::Error>>
{
    let file = stage_read_bed("chr19")?;

    let mut cmd = Command::cargo_bin("rastair")?;

    cmd.arg("bed2pat");
    cmd.args(["-r",
              "test_data/test.fasta.gz",
              "--nOT",
              "0,0,150,150",
              "--nOB",
              "0,0,150,150"]);
    cmd.arg(file.path());
    cmd.assert().success();

    let mut output = cmd.output().unwrap();
    let mut output_str = String::from_utf8_lossy(&output.stdout);
    let mut line_count = output_str.lines().count();
    assert!(line_count > 0);

    // there should be 30 reads reported in total
    let mut read_count: i32 = output_str.lines()
                                        .map(|line| {
                                            let elems: Vec<&str> = line.split("\t").collect();
                                            // this is safe because I filtered for len>9 above
                                            if elems.len() >= 4 {
                                                elems[3].parse::<i32>().unwrap_or_default()
                                            } else {
                                                0
                                            }
                                        })
                                        .sum();
    // I see that there are 3483 r1 reads on chr19:
    // rastair per-read -f 67 -F 3980 -l chr19 -r test_data/test.fasta test_data/test.bam | tail -n +2 | wc -l
    // I checked manually that 1 read is skipped because of indels
    // 5 read contain only 1 CpG which is a SNP (neither C nor T observed)
    // I therefore expect 3483 - 6 = 3477 reads reported
    assert_eq!(read_count, 3477);

    cmd = Command::cargo_bin("rastair")?;
    cmd.arg("bed2pat");
    cmd.args(["-r",
              "test_data/test.fasta.gz",
              "--nOT",
              "0,0,150,0",
              "--nOB",
              "0,0,150,0"]);
    cmd.arg(file.path());
    cmd.assert().success();

    output = cmd.output().unwrap();
    output_str = String::from_utf8_lossy(&output.stdout);
    line_count = output_str.lines().count();
    assert!(line_count > 0);

    // there should be 30 reads reported in total
    read_count = output_str.lines()
                           .map(|line| {
                               let elems: Vec<&str> = line.split("\t").collect();
                               // this is safe because I filtered for len>9 above
                               if elems.len() >= 4 {
                                   elems[3].parse::<i32>().unwrap_or_default()
                               } else {
                                   0
                               }
                           })
                           .sum();
    assert_eq!(read_count, 3477);

    cmd = Command::cargo_bin("rastair")?;
    cmd.arg("bed2pat");
    cmd.args(["-r",
              "test_data/test.fasta.gz",
              "--nOT",
              "0,0,0,150",
              "--nOB",
              "0,0,0,150"]);
    cmd.arg(file.path());
    cmd.assert().success();

    output = cmd.output().unwrap();
    output_str = String::from_utf8_lossy(&output.stdout);
    line_count = output_str.lines().count();
    assert!(line_count > 0);

    // there should be 30 reads reported in total
    read_count = output_str.lines()
                           .map(|line| {
                               let elems: Vec<&str> = line.split("\t").collect();
                               // this is safe because I filtered for len>9 above
                               if elems.len() >= 4 {
                                   elems[3].parse::<i32>().unwrap_or_default()
                               } else {
                                   0
                               }
                           })
                           .sum();
    assert_eq!(read_count, 3477);

    Ok(())
}
