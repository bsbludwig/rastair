//! Test big files with various flag combinations, very slow and depends on local test data
//!
//! Run this with
//!
//!     cargo test -p rastair -- --test-threads=1 --include-ignored big_file`
//!
//! and fetch a cup of coffee.

#![allow(non_snake_case)]
mod utils;
use utils::*;

macro_rules! test_files {
    ([ $(($name:ident, $bam:expr, $fasta:expr)),* $(,)? ]) => {
        $(
            pastey::paste! {
                mod [<big_file_ $name>] {
                    use super::*;

                    fn call(args: &[&str]) -> Result<()> {
                        std::process::Command::new("cargo")
                            .stdout(std::process::Stdio::null())
                            .args(["run", "--release", "--"])
                            .args(["call", $bam, &format!("--fasta-file={}", $fasta)])
                            .args(args)
                            .succeeds()
                    }

                    #[test]
                    #[ignore = "big file test"]
                    fn plain() -> Result<()> {
                        call(&[])
                    }

                    #[test]
                    #[ignore = "big file test"]
                    fn c() -> Result<()> {
                        call(&["-c"])
                    }

                    #[test]
                    #[ignore = "big file test"]
                    fn all() -> Result<()> {
                        call(&["--all"])
                    }

                    #[test]
                    #[ignore = "big file test"]
                    fn thresholds() -> Result<()> {
                        call(&["--thresholds"])
                    }

                    #[test]
                    #[ignore = "big file test"]
                    fn c_all() -> Result<()> {
                        call(&["-c", "--all"])
                    }

                    #[test]
                    #[ignore = "big file test"]
                    fn c_thresholds() -> Result<()> {
                        call(&["-c", "--thresholds"])
                    }

                    #[test]
                    #[ignore = "big file test"]
                    fn all_thresholds() -> Result<()> {
                        call(&["--all", "--thresholds"])
                    }

                    #[test]
                    #[ignore = "big file test"]
                    fn c_all_thresholds() -> Result<()> {
                        call(&["-c", "--all", "--thresholds"])
                    }
                }
            }
        )*
    };
}

test_files!([(test_bam, "tests/data/test.bam", "tests/data/test.fasta.gz"),]);
test_files!([(
    NA12878_aa_chr12,
    "tmp/taps/NA12878_aa_chr12.bam",
    "tmp/na12878/GRCh38_full_analysis_set_plus_decoy_hla.fa"
),]);
test_files!([(
    HG00096_chrom20,
    "tmp/1000genomes/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam",
    "tmp/1000genomes/hs37d5.fa.gz"
),]);
