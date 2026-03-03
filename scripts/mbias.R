#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  library(rmarkdown)
  library(argparser)
})

# Define RMarkdown file path (assumes it's in the same directory as this script)
rmd_file <- "QC_report.Rmd"

get_script_dir <- function(command_line_args) {
  command_line = paste(command_line_args, collapse=" ")
  dir = gsub(".*--file=(.*) --args.*", "\\1", as.character(command_line), perl = TRUE)
  return(dir)
}

# Create argument parser
parser <- arg_parser("Render RMarkdown document with genomic analysis parameters")

# Add optional arguments
parser <- add_argument(parser, "--bed", 
                       help = "Input per-read bed.gz file (must be tabix indexed)")

parser <- add_argument(parser, "--vcf", 
                       help = "Input vcf file (for cpg/gc bias plots)")

parser <- add_argument(parser, "--bam", 
                       help = "Input bam file (slower: will run rastair on the bam file to derive other data)")

parser <- add_argument(parser, "--no-vbias", 
                       help = "Do not generate vbias (faster)",
                       flag=TRUE)

parser <- add_argument(parser, "--no-gc", 
                       help = "Do not generate gc/cpg bias plots",
                       flag=TRUE)

parser <- add_argument(parser, "--reference", 
                       help = "Reference fasta sequence, fai indexed (required for vbias and gc plots)",
                       default = NULL)

parser <- add_argument(parser, "--region", 
                       help = "Genomic region (optional)",
                       default = NULL)

parser <- add_argument(parser, "--include-flag", 
                       help = "Include bitflag as integer (optional)",
                       type = "integer",
                       default = 3)

parser <- add_argument(parser, "--exclude-flag", 
                       help = "Exclude bitflag as integer (optional)", 
                       type = "integer",
                       default = 3852)

parser <- add_argument(parser, "--read-length", 
                       help = "Read length as integer (optional)",
                       type = "integer", 
                       default = NULL)

parser <- add_argument(parser, "--tabix-path", 
                       help = "Path to tabix executable (optional)",
                       default = "tabix")

parser <- add_argument(parser, "--bcftools-path", 
                       help = "Path to bedtools executable (optional)",
                       default = "bcftools")

parser <- add_argument(parser, "--rastair-path", 
                       help = "Path to rastair executable (optional)",
                       default = "rastair")

parser <- add_argument(parser, "--output-prefix", 
                       help = "Output path prefix (optional)",
                       default = ".")

parser <- add_argument(parser, "--threads", 
                       help = "Number of threads to use (only relevant when using bam input)",
                       default = 1)

parser <- add_argument(parser, "--wgbs", 
                       help = "Treat the input as inverted, ie mod=unmod and unmod=mod",
                       flag = TRUE)

# Parse arguments
args <- parse_args(parser)

# Validate bed file exists
if ((is.na(args$bed) || is.null(args$bed) || !file.exists(args$bed)) && (is.na(args$bam) || is.null(args$bam) || !file.exists(args$bam))) {
  cat("Error: either bam or bed input required\n")
  quit(status = 1)
}

# Create parameters list for RMarkdown
params_list <- list(
  is_WGBS = args$wgbs,
  region = args$region,
  output_dir = normalizePath(args$output_prefix),
  include_flags = args$include_flag,
  exclude_flags = args$exclude_flag,
  read_len = is.na(args$read_length),
  threads = args$threads,
  tabix = args$tabix_path,
  rastair = args$rastair_path,
  bcftools = args$bcftools_path,
  reference = ifelse(is.na(args$reference)||is.null(args$reference), "", normalizePath(args$reference)),
  plot_vbias = is.null(args$no_vbias)||args$no_vbias==FALSE,
  plot_gc = is.null(args$no_gc)||args$no_gc==FALSE,
  input_bgz = ifelse(is.na(args$bed)||is.null(args$bed), "", normalizePath(args$bed)),
  input_vcf = ifelse(is.na(args$vcf)||is.null(args$vcf), "", normalizePath(args$vcf)),
  input_bam = ifelse(is.na(args$bam)||is.null(args$bam), "", normalizePath(args$bam))
)

script_dir = get_script_dir(commandArgs())
script_dir = dirname(normalizePath(script_dir))
rmd_file = paste0(script_dir, "/", rmd_file)

cat("Rmd_file: ", rmd_file, "\n")
# Check if RMarkdown file exists
if (!file.exists(rmd_file)) {
  cat("Error: RMarkdown file not found:", rmd_file, "\n")
  cat("Make sure 'QC_report.Rmd' is in the same directory as this script.\n")
  quit(status = 1)
}

# Define output file path
output_file <- paste0(normalizePath(args$output_prefix), "/qc_report.html")
str(params_list)
# Render the RMarkdown document
cat("Rendering RMarkdown document...\n")
tryCatch({
  rmarkdown::render(
    input = rmd_file,
    output_file = output_file,
    params = params_list,
    intermediates_dir = getwd(),
    clean = TRUE,
    quiet = FALSE
  )
  
  cat("Report successfully generated:", output_file, "\n")
  
}, error = function(e) {
  cat("Error rendering RMarkdown document:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})

cat("Done!\n")