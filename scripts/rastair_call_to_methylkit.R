#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

#' Convert rastair call output to MethylKit format and save to file
#'
#' @param df A data frame from rastair call
#' @param output_filename Path to output TSV file
#'
#' @return Nothing, writes to disk

save_as_methylkit <- function(df, output_filename) {
  if (!all(c("#chr", "start", "mod", "unmod", "coverage") %in% colnames(df))) {
    stop("Input dataframe must contain columns: #chr, start, mod, unmod, coverage")
  }

# If file has no rows, exit gracefully
  if (nrow(df) == 0) {
    message(paste("Input df (input file) is empty — skipping conversion."))
    # Create an empty output file to satisfy downstream expectations, if needed:
    write_tsv(tibble(), output_filename)
    quit(status = 0)
  }

  # Create 'chrBase' as 'chr:start'
  df <- df %>%
    mutate(
      chrBase = paste0(`#chr`, ":", start),
      chr = `#chr`,
      base = start,
      strand = case_when(
        strand == "-" ~ "R",
        strand == "+" ~ "F",
        TRUE ~ "."
      ),
      coverage = coverage,
      freqT = 100 * unmod / (mod + unmod),
      freqC = 100 * mod / (mod + unmod)
    ) %>%
    select(chrBase, chr, base, strand, coverage, freqC, freqT)

  # Write to file
  write_tsv(df, output_filename)
}

#' Convert rastair call file to MethylKit format file
#'
#' @param input_file Path to input rastair call file
#' @param output_file Path to output MethylKit TSV file
#'
#' @return Nothing, writes to disk
rastair_to_methylkit <- function(input_file, output_file) {
  if (!file.exists(input_file)) {
    stop(paste("Input file", input_file, "does not exist."))
  }

  df <- read_tsv(input_file, col_types = cols(), show_col_types = FALSE)

  save_as_methylkit(df, output_file)
}

#' Command-line interface wrapper
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1 || length(args) > 2) {
    stop("Usage: Rscript rastair_to_methylkit.R input_file [output_file]")
  }

  input_file <- args[1]
  output_file <- if (length(args) == 2) {
    args[2]
  } else {
    paste0(basename(input_file), ".methylkit.tsv")
  }

  rastair_to_methylkit(input_file, output_file)
}

# Run main() if script is executed directly
if (sys.nframe() == 0) {
  main()
}
