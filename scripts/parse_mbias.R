#!/usr/bin/env Rscript

#' Parse mbias output table and write a transposed version
#'
#' @param input_file Path to input tab-delimited file with 'left' and 'right' columns
#' @param output_file Path to output CSV file (no headers, no rownames)
#'
#' @return Nothing, writes to disk

empty_file_message = "0,0,0,0\n0,0,0,0"

print_and_exit <- function(message, output_file) {
  writeLines(message, output_file)
  quit(status = 0)
}

parse_mbias <- function(input_file, output_file) {
  if (!file.exists(input_file)) {
    print(paste("Input file", input_file, "does not exist."))
    print_and_exit(empty_file_message, output_file)
  }

  if (file.info(input_file)$size == 0) {
    print(paste("Input file", input_file, "is empty (0 bytes)."))
    print_and_exit(empty_file_message, output_file)
  }

  # Read the input file
  df <- read.table(input_file, sep = "\t", header = TRUE)

  if (!all(c("left", "right") %in% colnames(df))) {
    print("Input file must contain 'left' and 'right' columns.")
    print_and_exit(empty_file_message, output_file)
  }

  # Select and transpose the specified columns
  df_transposed <- t(df[, c("left", "right")])

  # Create final matrix with 2 rows: left and right values
  df_final <- matrix(df_transposed, nrow = 2, byrow = TRUE)

  # Write to the specified output file
  write.table(df_final, output_file, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#' Run script as command-line tool
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 2) {
    stop("Usage: parse_mbias.R <input_file> <output_file>")
  }

  input_file <- args[1]
  output_file <- args[2]

  parse_mbias(input_file, output_file)
}

# Run main() if called as script
if (sys.nframe() == 0) {
  main()
}
