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

# Add positional argument for bed file
parser <- add_argument(parser, "bed_file", 
                       help = "Input bed.gz file (required)")

# Add optional arguments
parser <- add_argument(parser, "--no-vbias", 
                       help = "Do not generate vbias (faster)",
                       flag=TRUE)

parser <- add_argument(parser, "--reference", 
                       help = "Reference fasta sequence, fai indexed (required for vbias)",
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

parser <- add_argument(parser, "--output-prefix", 
                       help = "Output path prefix (optional)",
                       default = ".")

# Parse arguments
args <- parse_args(parser)

# Validate bed file exists
if (!file.exists(args$bed_file)) {
  cat("Error: Input bed.gz file does not exist:", args$bed_file, "\n")
  quit(status = 1)
}

# Print summary of parameters
cat("=== RMarkdown Report Generation ===\n")
cat("Parameters:\n")
cat("  Reference:", ifelse(is.null(args$reference), "Not specified", args$reference), "\n")
cat("  Plot vbias:", ifelse(is.null(args$no_vbias)||args$no_vbias==FALSE, "Yes", "No"), "\n")
cat("  Genomic region:", ifelse(is.null(args$region), "Not specified", args$region), "\n")
cat("  Include bitflag:", ifelse(is.null(args$include_flag), "Not specified", args$include_flag), "\n")
cat("  Exclude bitflag:", ifelse(is.null(args$exclude_flag), "Not specified", args$exclude_flag), "\n")
cat("  Read length:", ifelse(is.null(args$read_length), "Not specified", args$read_length), "\n")
cat("  Tabix path:", args$tabix_path, "\n")
cat("  Output prefix:", args$output_prefix, "\n")
cat("  Input bed.gz file:", normalizePath(args$bed_file), "\n")
cat("====================================\n\n")

# Create parameters list for RMarkdown
params_list <- list(
  region = args$region,
  output_dir = normalizePath(args$output_prefix),
  include_flags = args$include_flag,
  exclude_flags = args$exclude_flag,
  read_len = is.na(args$read_length),
  tabix = args$tabix_path,
  reference = params$reference,
  plot_vbias = is.null(params$no_vbias)||params$no_vbias==FALSE,
  input_bgz = normalizePath(args$bed_file)
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
output_file <- paste0(normalizePath(args$output_prefix), "/mbias.html")
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