#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Function to open a file given command-line args
file_from_args <- function(args) {
  if (length(args) == 0 || args[length(args)] == "-") {
    input <- file("stdin")
  } else {
    fileName <- args[length(args)]
    if (grepl("\\.gz$", fileName)) {
      input <- gzfile(fileName)
    } else {
      input <- file(fileName)
    }
  }
  return(input)
}

# Function to parse command line args - I'm trying to avoid unnecessary dependencies, so not
# using optargs or other packages here
value_from_args <- function(args, option_strings, default="", as.number=FALSE, is.boolean=FALSE) {
  if (!is.vector(args) || length(args) == 0) {
    stop("Missing options")
  }
  if (length(option_strings) == 0) {
    stop("Missing options")
  }

  value = ifelse(is.boolean, FALSE, default)
  while(length(args) > 0) {
    next_arg = args[1]
    args <- args[-1]
    if (next_arg %in% option_strings)
    {
      if (is.boolean) {
        return(TRUE)
      }

      if (length(args) > 0)
      {
        value = args[1]
        if (grepl("^\\-", value))
        {
          warn(paste("Missing value for", option_strings[1], "option\n"));
        }
      }
      break
    }
  }
  return(ifelse(as.number, as.numeric(value), value))
}

calculate_summary <- function(methylation) {
  summary_beta <- summary(methylation$beta_est)
  names(summary_beta) <- c("min", "q1", "median","mean", "q3", "max", "NAs")
  summary_beta["sd"] <- sd(methylation$beta_est, na.rm = TRUE)
  summary_beta["sum_mean"] <- sum(methylation$mod)/(sum(methylation$mod)+sum(methylation$unmod))
  return(summary_beta)
}

main <- function() {
  methylation <- read.delim(file_from_args(c("tests/data/lambda_calls.bed.gz")), stringsAsFactors=TRUE, na.strings = c("NA","NaN",".",""))
  summary_beta <- calculate_summary(methylation)
  cat(paste(names(summary_beta), sep="\t"))
  cat(summary_beta, sep="\t")
}

# Run as script unless explicitly asked not to, for unit testing
if (getOption('run.main', default=TRUE)) {
  main()
}
