// Externals
use log::error;
use clap::{arg, command, value_parser, Parser, Subcommand};
use clio::*;

use rastair::operations::{count_variants, count_reads};
use rastair::sequence_segment::run_finder;

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
        #[arg(value_name="BAM_FILE", value_parser=value_parser!(ClioPath).exists().is_file())]
        bam_file: ClioPath,

        /// A sorted and indexed (via samtools faidx) fasta file. Note that bgzip compressed files are NOT currently supported
        #[arg(short='r', long, value_name="FASTA_FILE", required=true, value_parser=value_parser!(ClioPath).exists().is_file())]
        fasta_file: ClioPath,

        /// Minimum mapping quality per aligned read [default: 1]
        #[arg(short='q', long)]
        min_mapq: Option<u8>,

        /// Minimum base quality per base in a read [default: 10]
        #[arg(short='Q', long, value_parser = clap::value_parser!(u8).range(1..))]
        min_baseq: Option<u8>,

        /// Limit depth at highly covered positions to improve performance [default: 500]
        #[arg(short='x', long, value_parser = clap::value_parser!(u32).range(1..))]
        max_depth: Option<u32>,

        /// number of reference positions processed in-memory at once [default: 100000]
        #[arg(short='s', long, value_parser = clap::value_parser!(u32).range(1..))]
        chunk_size: Option<u32>,

        /// Include reads that match all of these bit-flags (as decimal) [default: 3]
        #[arg(short='f', long)]
        required_flags: Option<u16>,

        /// Exclude reads matching any of these bit-flags (as decimal) [default: 3852]
        #[arg(short='F', long)]
        excluded_flags: Option<u16>,

        /// Exclude reads where the orientation cannot be unambiguously determined [default: false]
        #[arg(long, action=clap::ArgAction::SetTrue)]
        exclude_ambiguous: Option<bool>,

        /// Soft-trim the OT by that many bases from either end [default: 0,0,0,0]
        #[arg(long="nOT")]
        n_ot: Option<String>,

        /// Soft-trim the OT by that many bases from either end [default: 0,0,0,0]
        #[arg(long="nOB")]
        n_ob: Option<String>,

        /// Number of threads to use inside htslib for decompression [default: 1]
        #[arg(long, value_parser = clap::value_parser!(u8).range(1..))]
        read_threads: Option<u8>,

        /// Number of threads to use for variant calling
        #[arg(short='@', long, value_parser = clap::value_parser!(u8).range(1..))]
        threads: Option<u8>,
    },
    /// Call methylation per read. This will produce a bed file that list the methylation status of all CpGs
    /// in every read that overlaps a CpG, plus some other metadata
    PerRead
    {
        /// A sorted and indexed bam file
        #[arg(value_name="BAM_FILE", value_parser=value_parser!(ClioPath).exists().is_file())]
        bam_file: ClioPath,

        /// A sorted and indexed (via samtools faidx) fasta file. Note that bgzip compressed files are NOT currently supported
        #[arg(short='r', long, value_name="FASTA_FILE", required=true, value_parser=value_parser!(ClioPath).exists().is_file())]
        fasta_file: ClioPath,

        /// Minimum mapping quality per aligned read [default: 1]
        #[arg(short='q', long)]
        min_mapq: Option<u8>,

        /// number of reference positions processed in-memory at once [default: 100000]
        #[arg(short='s', long, value_parser = clap::value_parser!(u32).range(1..))]
        chunk_size: Option<u32>,

        /// Include reads that match all of these bit-flags (as decimal) [default: 3]
        #[arg(short='f', long)]
        required_flags: Option<u16>,

        /// Exclude reads matching any of these bit-flags (as decimal) [default: 3852]
        #[arg(short='F', long)]
        excluded_flags: Option<u16>,

        /// Report reads with no CpGs in them
        #[arg(short='A', long, action=clap::ArgAction::SetTrue)]
        all_reads: Option<bool>,

        /// Number of threads to use inside htslib for decompression [default: 1]
        #[arg(long, value_parser = clap::value_parser!(u8).range(1..))]
        read_threads: Option<u8>,
    },
    /// print a map of all CpGs in a fasta file
    MapCpgs
    {
        /// An indexed fasta file
        #[arg(value_name="FASTA_FILE", value_parser=value_parser!(ClioPath).exists().is_file())]
        fasta_file: ClioPath,

        /// number of reference positions processed in-memory at once [default: 100000]
        #[arg(short='s', long, value_parser = clap::value_parser!(u32).range(1..))]
        chunk_size: Option<u32>,
    },
}

fn main()
{
    reset_sigpipe();

    let exit_code = real_main();
    std::process::exit(exit_code);
}

fn real_main() -> i32
{
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
            exclude_ambiguous,
            n_ot,
            n_ob,
            read_threads,
            threads }) =>
            {
                // TODO move the unboxing of Options to here
                // instead of doing that inside run_caller
                match count_variants::run_caller(&bam_file.to_path_buf(),
                    &fasta_file.to_path_buf(),
                    min_mapq,
                    min_baseq,
                    max_depth,
                    chunk_size,
                    required_flags,
                    excluded_flags,
                    exclude_ambiguous,
                    n_ot,
                    n_ob,
                    read_threads,
                    threads )
                {
                    Ok(()) => 0,
                    Err(e)  =>
                    {
                        error!("Error running caller: {}", e);
                        1
                    }
                }
            },
        Some(Commands::MapCpgs {
                fasta_file,
                chunk_size }) =>
                {
                    let step_size = chunk_size.unwrap_or_default();

                    match run_finder(&fasta_file.to_path_buf(), step_size as usize)
                    {
                        Ok(_) => 0,
                        Err(e) => {
                            error!("Failed to run cpg_finder: {}", e);
                            1
                        }
                    }
                }
        Some(Commands::PerRead {
                bam_file,
                fasta_file,
                min_mapq,
                chunk_size,
                required_flags,
                excluded_flags,
                all_reads,
                read_threads }) =>
                {
                    match count_reads::run_caller(
                        &bam_file.to_path_buf(),
                        &fasta_file.to_path_buf(),
                        min_mapq,
                        chunk_size,
                        required_flags,
                        excluded_flags,
                        all_reads,
                        read_threads,)
                    {
                        Ok(()) => 0,
                        Err(e)  =>
                        {
                            error!("Error running caller: {}", e);
                            1
                        }
                    }
                }
        None => 0
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