// Externals
use log::{info, debug, trace, warn};
use clap::{arg, command, value_parser};
use rust_htslib::bam::IndexedReader;
use std::path::PathBuf;
use std::process;

// Library imports
use rastair::sequence_segment::SequenceSegmentIterator;

fn main() {
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
    let bam = match matches.get_one::<PathBuf>("BAM_FILE") {
        Some(bam_path) => IndexedReader::from_path(bam_path).expect("Error reading bam file"),
        None    => {
            warn!("Piping bam input not yet implemented");
            process::exit(1);
        }
    };

    /* Step through all chromosomes, at a fixed step size, and 
     * identify the locations of all CpG positions
     */ 
    let seq_iter = SequenceSegmentIterator::from_file_with_stepsize(fasta_path, 1000).expect("Error creating FASTA iterator");
    for segment in seq_iter {
        println!("{}\t{}\t{}\t{}", segment.contig, segment.start, segment.stop, std::str::from_utf8(&segment.sequence).unwrap() );
    }
    /* Fetch the pileup for the region from the bam file, and go
     * through all CpG positions, performing whatever calculation
     * needs to be performed. Stream the output to a writer that
     * writes the results to STDOUT or some file.
     */
}
