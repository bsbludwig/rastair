// Externals
use log::debug;
use clap::{arg, command, value_parser};

use std::io::{stdout, Write};
use std::path::PathBuf;

use rastair::operations::{VariantCounterConfig, VariantCounter};

fn main()
{
    reset_sigpipe();

    let matches = command!() // requires `cargo` feature
        .arg(
            arg!(<BAM_FILE> "Path to bam file")
            .required(true)
            .value_parser(value_parser!(PathBuf))
        )
        .arg(
            arg!(
                -f --fasta <FASTA_FILE> "Reference fasta file. Must be indexed (expecting .fai file to exist)"
            )
            // We don't have syntax yet for optional options, so manually calling `required`
            .required(true)
            .value_parser(value_parser!(PathBuf)),
        )
        .arg(arg!(
            -v --verbosity ... "Increase output verbosity. Can be repeated multiple times for more detail"
        ))
        .get_matches();

    /* Initialise logging */
    let verbose = matches
    .get_one::<u8>("verbosity")
    .expect("Count's are defaulted");

    stderrlog::new()
        .module(module_path!())
        .verbosity(verbose.clone() as usize)
        .timestamp(stderrlog::Timestamp::Second)
        .init()
        .unwrap();

    /* Read fasta index, and open fasta file for tokenising */
    let fasta_path = matches.get_one::<PathBuf>("fasta").expect("Argument missing");
    debug!("Reading fasta and index from {}", fasta_path.display());
    /* Open the bam file */
    let bam_path= matches.get_one::<PathBuf>("BAM_FILE").expect("Error getting BAMFILE param");
    let config = VariantCounterConfig::with_paths(fasta_path, bam_path).unwrap();
    let counter = VariantCounter::with_config(config).unwrap();

    let mut lock = stdout().lock();
    for cpgs in counter
    {
        for cpg in cpgs
        {
            if cpg.ref_base == b'C'
            { // C
                writeln!(lock, "{}\t{}\t{}\t{}\t{}", cpg.contig, cpg.pos, cpg.pos+1, cpg.top.c, cpg.top.t).unwrap();
            }
            else
            { // G
                writeln!(lock, "{}\t{}\t{}\t{}\t{}", cpg.contig, cpg.pos, cpg.pos+1, cpg.bottom.g, cpg.bottom.a).unwrap();
            }
        }
    }
}

/*
 * some super-hacky sh*t to make this behave like a normal unix program and quit when the pipe ends
 */
#[cfg(unix)]
fn reset_sigpipe()
{
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}

#[cfg(not(unix))]
fn reset_sigpipe()
{
    // no-op
}