#!/usr/bin/env Rscript

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

annotate_with_context <- function(methylation, fasta_file) {
  suppressMessages(suppressWarnings(library("Rsamtools")))
  # use only the G in the CpG - coordinates are in bed format originally,
  # whereas GRanges uses 1-offset, right closed intervals. This should now
  # refer to each the CpG pair
  subset_regions <- makeGRangesFromDataFrame(methylation[methylation$strand=="-",],
                                             keep.extra.columns = FALSE,
                                             ignore.strand = TRUE,
                                             seqnames.field = "X.chr",
                                             seqinfo = seqinfo(scanFaIndex(fasta_file)))
  # extend the regions and cut out from fasta file
  start(subset_regions) <- start(subset_regions)-1
  end(subset_regions) <- end(subset_regions)+1
  # ensure all regions are inside the chgromosome
  subset_regions <- trim(subset_regions)
  values(subset_regions) <- DataFrame( sequence=as.vector(scanFa(fasta_file, subset_regions, as="DNAStringSet")))

  # summarise beta per CpG pair
  c_meth <- subset(methylation, strand=="+")
  g_meth <- subset(methylation, strand=="-")
  # make coordinates the same
  c_meth$start <- c_meth$start + 1
  c_meth$end <- c_meth$end + 1
  merged_meth <- merge(c_meth, g_meth, by=c("X.chr", "start", "end"))
  merged_meth <- merged_meth[, c("X.chr", "start", "end", "mod.x", "unmod.x", "mod.y", "unmod.y")]
  merged_meth$beta <- (merged_meth$mod.x + merged_meth$mod.y)/(merged_meth$mod.x + merged_meth$mod.y+merged_meth$unmod.x + merged_meth$unmod.y)
  subset_regions <- as.data.frame(subset_regions)
  subset_regions$start <- subset_regions$start+1
  subset_regions$end <- subset_regions$end-1
  subset_regions <- merge(subset_regions, merged_meth, by.x=c("seqnames", "start", "end"), by.y=c("X.chr", "start", "end"))
  names(subset_regions) <- gsub("seqnames", "X.chr", names(subset_regions))

  return(subset_regions[, c("X.chr", "start", "end", "sequence", "beta")])
}

summarise_by_context <- function(context_table) {

}

plot_context_stats <- function(context_table, args=c()) {
  # Only load if necessary
  suppressMessages(suppressWarnings(library("ggplot2", "gtable")))

  out_filename <- value_from_args(args, c("-o", "--output-file"), default = "conversion_context.png")
  fig_width=value_from_args(args, c("-x", "--width"), default = 1200, as.number = TRUE)
  fig_height=value_from_args(args, c("-y", "--height"), default = 800, as.number = TRUE)

  if (value_from_args(args, c("-j", "--jpeg"), is.boolean = TRUE)) {
    jpg(out_filename, width=fig_width, height=fig_height)
  } else if (value_from_args(args, c("--pdf"), is.boolean = TRUE)) {
    resolution=150
    pdf(out_filename, width=fig_width/resolution, height=fig_height/resolution, onefile = TRUE, bg="white")
  } else if (value_from_args(args, c("--svg"), is.boolean = TRUE)) {
    resolution=150
    svg(out_filename, width=fig_width/resolution, height=fig_height/resolution, onefile = TRUE)
  }
  else {
    png(out_filename, width=fig_width, height=fig_height)
  }

  g_box <- ggplot(subset(context_table, nchar(sequence)==4), aes(x=sequence, y=beta, color=substr(sequence,1,1))) + geom_boxplot(outlier.size = 1, outlier.alpha = 0.5)
  g_box <- g_box + labs(color = "5' Context") +
           xlab(NULL) +
           theme_light()

  g_bar <- ggplot(subset(context_table, nchar(sequence)==4), aes(x=sequence, fill=substr(sequence,1,1))) + geom_bar()
  g_bar <- g_bar + labs(fill = "5' Context") +
           xlab(NULL) +
           theme_light()

  gtbl <- gtable_col(name="conversion", grobs=list(Coverage=ggplotGrob(g_bar), Conversion=ggplotGrob(g_box)), heights=unit(c(1,2)/3, "npc"), width=unit(1, "npc"))
  grid.newpage()
  grid.draw(gtbl)
  invisible(dev.off())
}

calculate_summary <- function(methylation) {
  summary_beta <- summary(methylation$beta_est)
  if(length(summary_beta)==6) {
    names(summary_beta) <- c("min", "q1", "median","mean", "q3", "max")
    summary_beta["missing"] <- 0
  } else {
    names(summary_beta) <- c("min", "q1", "median","mean", "q3", "max", "missing")
  }
  summary_beta["sd"] <- sd(methylation$beta_est, na.rm = TRUE)
  summary_beta["sum_mean"] <- sum(methylation$mod)/(sum(methylation$mod)+sum(methylation$unmod))
  return(summary_beta)
}

main <- function(args) {
  if (value_from_args(args, c("-h", "--help"), is.boolean = TRUE))
  {
    cat(
      "Summarise a methylation call table into averages and ranges to estimate conversion statistics.",
      "",
      "Usage: calculate_conversion.R [OPTIONS] <INFILE>",
      "",
      "Arguments:",
      "\t<INFILE> Path to a file generated by rastair call, or `-` to read from STDIN",
      "",
      "Options:",
      "\t-c/--plot-context Generate an boxplot of beta values by sequence context surrounding the CpG",
      "\t-r/--fasta-file REFERENCE Path to the indexed reference fasta file (can be bgzip-compressed). Only needed with --plot-context",
      "\t-o/--output-file OUTFILENAME Name of the output image file if the context plot is requested [conversion_context.png]",
      "\t-x/--width WIDTH Width of output image in pixels [1200]",
      "\t-y/--height HEIGHT Height of output image in pixels [800]",
      "\t-j/--jpeg Generate JPEG output",
      "\t--pdf Generate PDF output",
      "\t--svg Generate SVG output",
      "\t-h/--help Show this help message",
      sep = "\n"
    )
    quit(save="no", status=0)
  }

  methylation <- read.delim(file_from_args(args), stringsAsFactors=TRUE, na.strings = c("NA","NaN",".",""))

  # Calculate the summary over the input data and print
  summary_beta <- calculate_summary(methylation)
  cat(names(summary_beta), sep="\t", "\n")
  cat(summary_beta, sep="\t", "\n")

  # if requested, also generate a figure of per-context conversion
  if(value_from_args(args, c("-c", "--plot-context"), is.boolean = TRUE)) {
    fasta_file = value_from_args(args, c("-r", "--fasta-file"))
    if (nchar(fasta_file) == 0) {
      stop("Missing fasta file, use -r")
    }
    annotated_table <- annotate_with_context(methylation, fasta_file)
    plot_context_stats(annotated_table, args)
  }
}

# Run as script unless explicitly asked not to, for unit testing
if (getOption('run.main', default=TRUE)) {
  args <- commandArgs(trailingOnly = TRUE)
  main(args)
}
