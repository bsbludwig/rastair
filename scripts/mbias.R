#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  library(rmarkdown)
  library(optparse)
})

# Define RMarkdown file path (assumes it's in the same directory as this script)
rmd_file <- "QC_report.Rmd"

# Define command line options
option_list <- list(
  make_option(c("-r", "--region"), 
              type = "character", 
              default = NULL,
              help = "Genomic region (optional) [default: %default]"),
  
  make_option(c("-i", "--include_flag"), 
              type = "integer", 
              default = NULL,
              help = "Include bitflag as integer (optional) [default: %default]"),
  
  make_option(c("-e", "--exclude_flag"), 
              type = "integer", 
              default = NULL,
              help = "Exclude bitflag as integer (optional) [default: %default]"),
  
  make_option(c("-l", "--read_length"), 
              type = "integer", 
              default = NULL,
              help = "Read length as integer (optional) [default: %default]"),
  
  make_option(c("-t", "--tabix_path"), 
              type = "character", 
              default = "tabix",
              help = "Path to tabix executable (optional) [default: %default]"),
  
  make_option(c("-o", "--output_prefix"), 
              type = "character", 
              default = ".",
              help = "Output path prefix (optional) [default: %default]"),
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list,
                           description = "Generate a methylation bias report as HTML")
opt <- parse_args2(opt_parser)

bed_file = args[1]
# Check for required arguments
if (is.null(bed_file)) {
  cat("Error: Input bed.gz file is required!\n")
  print_help(opt_parser)
  quit(status = 1)
}

# Validate bed file exists
if (!file.exists(bed_file)) {
  cat("Error: Input bed.gz file does not exist:", bed_file, "\n")
  quit(status = 1)
}

# Print summary of parameters
cat("=== RMarkdown Report Generation ===\n")
cat("Parameters:\n")
cat("  Genomic region:", ifelse(is.null(opt$region), "Not specified", opt$region), "\n")
cat("  Include bitflag:", ifelse(is.null(opt$include_flag), "Not specified", opt$include_flag), "\n")
cat("  Exclude bitflag:", ifelse(is.null(opt$exclude_flag), "Not specified", opt$exclude_flag), "\n")
cat("  Read length:", ifelse(is.null(opt$read_length), "Not specified", opt$read_length), "\n")
cat("  Tabix path:", opt$tabix_path, "\n")
cat("  Output prefix:", opt$output_prefix, "\n")
cat("  Input bed.gz file:", opt$bed_file, "\n")
cat("====================================\n\n")

# Create parameters list for RMarkdown
params_list <- list(
  genomic_region = opt$region,
  include_flags = opt$include_flag,
  exclude_flags = opt$exclude_flag,
  read_len = opt$read_length,
  tabix = opt$tabix_path,
  output_dir = opt$output_prefix,
  input_bgz = opt$bed_file
)

# Check if RMarkdown file exists
if (!file.exists(rmd_file)) {
  cat("Error: RMarkdown file not found:", rmd_file, "\n")
  cat("Make sure 'genomic_analysis_report.Rmd' is in the same directory as this script.\n")
  quit(status = 1)
}

# Define output file path
output_file <- paste0(opt$output_prefix, "_report.html")

# Render the RMarkdown document
cat("Rendering RMarkdown document...\n")
tryCatch({
  rmarkdown::render(
    input = rmd_file,
    output_file = output_file,
    params = params_list,
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