#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

pop <- function(args) {
  if (is.vector(args) && length(args) > 0)
  {
    return(args[-length(args)])
  }
  return(c())
}

`%between%` <- function(val, range) {
  if (is.null(val) || is.null(range) || !is.vector(range) || length(range) != 2 || range[2]<=range[1])
  {
    return(NULL)
  }
  return(val >= range[1] && val <= range[2])
}

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

generate_figure <- function(mbias_table, args, cutoffs=NULL) {
  # Only load if necessary
  library(ggplot2)

  out_filename <- value_from_args(args, c("-o", "--output-file"), default = "mbias.png")
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
  max_obs = max(mbias_table$X.pos)
  g <- ggplot(mbias_table, aes(x=X.pos, y=beta, color=read_pair)) +
    geom_line() +
    geom_ribbon(aes(ymin=lower_ci, ymax=upper_ci, fill=read_pair), alpha=0.2) +
    labs(color = "Read") +
    xlab("Position in read (5'->3')") +
    scale_color_discrete(labels=c('First in pair', 'Second in pair')) +
    scale_fill_discrete(guide="none") +
    theme_light()
  if(!is.null(cutoffs)) {
    g = g + geom_vline(aes(xintercept=left, color=read_pair), data=cutoffs)
    g = g + geom_vline(aes(xintercept=max_obs-right, color=read_pair), data=cutoffs)
  }
  g = g + facet_wrap(vars(strand), nrow = 2, as.table = FALSE)

  print(g)
  invisible(dev.off())
}

# For a subset of observations of read orientation on one strand,
# calculate the suggested cutoff points
calculate_cutoff <- function(mbias_table) {
  # Ensure we process things in orders
  mbias_table <- mbias_table[order(mbias_table$X.pos),]
  # get the mean beta of the middle 60% of the data
  n_obs <- dim(mbias_table)[1]
  if (n_obs < 3) {
    stop("Not enough data points to calculate cutoffs")
  }
  upper_lim <- min(ceiling(n_obs * 0.8), n_obs-1)
  lower_lim <- max(floor(n_obs * 0.2), 2)
  mean_middle <- mean(mbias_table[lower_lim:upper_lim,]$beta)
  upper_ci_max <- max(mbias_table[lower_lim:upper_lim,]$upper_ci)
  lower_ci_min <- min(mbias_table[lower_lim:upper_lim,]$lower_ci)

  # Search from the inside out for the first point left/right that is
  # either outside the 95% CI bounds of the middle, or the mean of the
  # middle is outside the CI of the point, or it's more than 0.05 away
  # from the mean of the middle
  check_match <- function(row, mean_middle, lower_ci_min, upper_ci_max) {
    if (!row$beta %between% c(lower_ci_min, upper_ci_max)
        && abs(row$beta - mean_middle) > 0.05
        && !mean_middle %between% c(row$lower_ci, row$upper_ci)) {
      return(TRUE)
    }
    return(FALSE)
  }

  left_cutoff = 0
  right_cutoff = 0

  middle = ceiling(n_obs/2)
  for (pos in middle:1) {
    if (check_match(mbias_table[pos,], mean_middle, lower_ci_min, upper_ci_max)) {
      left_cutoff = pos
      break
    }
  }
  for (pos in 1:(n_obs-middle)) {
    if (check_match(mbias_table[middle+pos,], mean_middle, lower_ci_min, upper_ci_max)) {
      right_cutoff = n_obs - (middle+pos) + 1
      break
    }
  }
  return(c(left=left_cutoff, right=right_cutoff))
}

post_process_mbias_table <- function(mbias_table) {
  mbias_table$strand <- as.factor(t(simplify2array(strsplit(as.character(mbias_table$type), "/", fixed=TRUE)))[,1])
  mbias_table$read_pair <- as.factor(t(simplify2array(strsplit(as.character(mbias_table$type), "/", fixed=TRUE)))[,2])

  # Flip reverse reads so that the order of beta is always 5' to 3' relative to the sequencing direction
  mbias_table$X.pos[which(mbias_table$strand=="OT" & mbias_table$read_pair=="2")] <- rev(subset(mbias_table, strand=="OT" & read_pair=="2")$X.pos)
  mbias_table$X.pos[which(mbias_table$strand=="OB" & mbias_table$read_pair=="1")] <- rev(subset(mbias_table, strand=="OB" & read_pair=="1")$X.pos)
  return(mbias_table)
}

add_CI <- function(mbias_table) {
  conf_ints <- as.data.frame(t(apply(subset(mbias_table, select=c("mod", "unmod")), 1, function(row){ binom.test(row[1], row[1]+row[2], conf.level = 0.999)$conf.int })))
  colnames(conf_ints) <- c("lower_ci", "upper_ci")
  mbias_table <- cbind(mbias_table, conf_ints)

  return(mbias_table)
}

main <- function() {
  if (value_from_args(args, c("-h", "--help"), is.boolean = TRUE))
  {
    cat(
      "Generate a figure from the output of rastair mbias and calculate suggested -nOT/-nOB parameters.",
      "",
      "Usage: plot_mbias.R [OPTIONS] -o <OUTFILENAME> <INFILE>",
      "",
      "Arguments:",
      "\t<INFILE> Path to an output file generated by rastair mbias, or `-` to read from STDIN",
      "",
      "Options:",
      "\t-o/--output-file OUTFILENAME Name of the output image file [mbias.png]",
      "\t-n/--no-plot Do not generate image output, just calculate soft-mask suggestions",
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

  input <- file_from_args(args)
  args <- pop(args)

  mbias_table <- read.delim(input, stringsAsFactors = TRUE)

  if (dim(mbias_table)[1] == 0) {
    message("Empty input table")
    quit(save="no", status=2)
  }
  # Add strand and read_pair columns
  mbias_table <- post_process_mbias_table(mbias_table)

  # add confidence intervals
  mbias_table <- add_CI(mbias_table)

  # Calculate optimal cutoffs
  cutoffs <- rbind(data.frame(strand="OT", read_pair="1", as.list(calculate_cutoff(subset(mbias_table, strand=="OT" & read_pair == "1")))),
                   data.frame(strand="OT", read_pair="2", as.list(calculate_cutoff(subset(mbias_table, strand=="OT" & read_pair == "2")))),
                   data.frame(strand="OB", read_pair="1", as.list(calculate_cutoff(subset(mbias_table, strand=="OB" & read_pair == "1")))),
                   data.frame(strand="OB", read_pair="2", as.list(calculate_cutoff(subset(mbias_table, strand=="OB" & read_pair == "2")))))
  cutoffs$strand <- as.factor(cutoffs$strand)
  cutoffs$read_pair <- as.factor(cutoffs$read_pair)

  if (!value_from_args(args, c("-n", "--no-plot"), is.boolean = TRUE)) {
    generate_figure(mbias_table, args, cutoffs)
  }

  write.table(cutoffs, stdout(), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
}

# Run as script unless explicitly asked not to, for unit testing
if (getOption('run.main', default=TRUE)) {
  main()
}
