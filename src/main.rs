// Externals
use log::error;
use clap::{arg, command, value_parser, Parser, Subcommand};
use std::path::PathBuf;

use rastair::operations::count_variants::run_caller;

#[derive(Parser)]
#[command(author="Benjamin Schuster-Boeckler", version, about, long_about=None, arg_required_else_help = true)]
struct Cli 
{
    #[arg(short, long, action=clap::ArgAction::Count)]
    verbosity: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}


#[derive(Subcommand)]
enum Commands 
{
    /// Call methylation on a bam file
    Call 
    {
        /// A sorted and indexed bam file
        #[arg(value_name="BAM_FILE", value_parser=value_parser!(PathBuf))]
        bam_file: PathBuf,

        /// A sorted and indexed (via samtools faidx) fasta file. Note that bgzip compressed files are NOT currently supported
        #[arg(short='r', long, value_name="FASTA_FILE", required=true, value_parser=value_parser!(PathBuf))]
        fasta_file: PathBuf,

        /// Minimum mapping quality per aligned read [default: 1]
        #[arg(short='q', long)]
        min_mapq: Option<u8>,
        
        /// Minimum base quality per base in a read [default: 10]
        #[arg(short='Q', long, value_parser = clap::value_parser!(u8).range(1..))]
        min_baseq: Option<u8>,
        
        /// Limit depth at highly covered positions to improve performance [default: 500]
        #[arg(short='x', long, value_parser = clap::value_parser!(u8).range(1..))]
        max_depth: Option<u32>,
        
        /// number of reference positions processed in-memory at once [default: 100000]
        #[arg(short='s', long, value_parser = clap::value_parser!(u8).range(1..))]
        chunk_size: Option<usize>,
        
        /// Include reads that match all of these bit-flags (as decimal) [default: 3]
        #[arg(short='f', long)]
        required_flags: Option<u16>,
        
        /// Exclude reads matching any of these bit-flags (as decimal) [default: 3852]
        #[arg(short='F', long)]
        excluded_flags: Option<u16>,

        /// Soft-trim the OT by that many bases from either end [default: 0,0,0,0]
        #[arg(long="nOT")]
        n_ot: Option<String>,
        
        /// Soft-trim the OT by that many bases from either end [default: 0,0,0,0]
        #[arg(long="nOB")]
        n_ob: Option<String>,

        /// Number of threads to use inside htslib for decompression [default: 1]
        #[arg(long)]
        read_threads: Option<usize>,
        
    },
}

fn main()
{
    reset_sigpipe();

    let cli = Cli::parse();

    /* Initialise logging */
    stderrlog::new()
        .module(module_path!())
        .verbosity(cli.verbosity as usize)
        .timestamp(stderrlog::Timestamp::Second)
        .init()
        .unwrap();    
    
    match &cli.command
    {
        Some(Commands::Call { 
            bam_file, 
            fasta_file, 
            min_mapq, 
            min_baseq, 
            max_depth, 
            chunk_size, 
            required_flags, 
            excluded_flags,
            n_ot,
            n_ob,
            read_threads }) => 
            {
                match run_caller(bam_file, 
                    fasta_file, 
                    min_mapq, 
                    min_baseq, 
                    max_depth, 
                    chunk_size, 
                    required_flags, 
                    excluded_flags,
                    n_ot,
                    n_ob,
                    read_threads ) 
                {
                    Ok(()) => (),
                    Err(e)  => 
                    {
                        error!("Error running caller: {}", e)
                    }
                }
            },
        None => ()
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