// Externals
use log::{debug, error};
use clap::{arg, command, value_parser, Parser, Subcommand};
use std::error::Error;
use std::io::{stdout, Write};
use std::path::PathBuf;

use rastair::operations::{VariantCounterConfig, VariantCounter};

#[derive(Parser)]
#[command(author, version, about, long_about=None, arg_required_else_help = true)]
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
    Call 
    {
        #[arg(value_name="BAM_FILE", value_parser=value_parser!(PathBuf))]
        bam_file: PathBuf,

        #[arg(short='r', long, value_name="FASTA_FILE", required=true, value_parser=value_parser!(PathBuf))]
        fasta_file: PathBuf,

        #[arg(short='q', long)]
        min_mapq: Option<u8>,
        #[arg(short='Q', long, value_parser = clap::value_parser!(u8).range(1..))]
        min_baseq: Option<u8>,
        #[arg(short='x', long, value_parser = clap::value_parser!(u8).range(1..))]
        max_depth: Option<u32>,
        #[arg(short='s', long, value_parser = clap::value_parser!(u8).range(1..))]
        chunk_size: Option<usize>,
        #[arg(short='f', long)]
        required_flags: Option<u16>,
        #[arg(short='F', long)]
        excluded_flags: Option<u16>
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
            excluded_flags }) => 
            {
                match run_caller(bam_file, 
                    fasta_file, 
                    min_mapq, 
                    min_baseq, 
                    max_depth, 
                    chunk_size, 
                    required_flags, 
                    excluded_flags ) 
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

fn run_caller(
    bam_path: &PathBuf,
    fasta_path: &PathBuf,
    mapq_option: &Option<u8>,
    baseq_option: &Option<u8>,
    max_depth_option: &Option<u32>,
    chunk_size_option: &Option<usize>,
    req_flags_option: &Option<u16>,
    excl_flags_option: &Option<u16>) -> Result<(), Box<dyn Error>> 
{
    /* Read fasta index, and open fasta file for tokenising */
    debug!("Reading fasta and index from {}", fasta_path.display());
    
    let mut config = VariantCounterConfig::with_paths(fasta_path, bam_path).unwrap();
    if let Some(min_mapq) = mapq_option {
        config.min_mapq = *min_mapq;
    }
    if let Some(min_baseq) = baseq_option {
        config.min_baseq = *min_baseq;
    }
    if let Some(max_depth) = max_depth_option {
        config.max_depth = *max_depth;
    }
    if let Some(cs) = chunk_size_option {
        config.chunk_size = *cs;
    }
    if let Some(flags) = req_flags_option {
        config.required_flags = *flags;
    }
    if let Some(flags) = excl_flags_option {
        config.excluded_flags = *flags;
    }
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
    Ok(())
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