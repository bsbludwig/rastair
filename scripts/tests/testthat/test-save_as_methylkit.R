options(run.main=FALSE)
source("../../rastair_call_to_methylkit.R")

test_that("save_as_methylkit correctly converts example input", {

  # Create example dataframe as given
  df <- data.frame(
    "#chr" = rep("J02459.1", 7),
    start = c(3, 4, 6, 7, 12, 13, 14),
    end = c(4, 5, 7, 8, 13, 14, 15),
    name = ".",
    beta_est = c(1.00, 0.94, 0.50, 0.95, 0.75, 0.93, 1.00),
    strand = c("+", "-", "+", "-", "+", "-", "+"),
    unmod = c(0, 1, 1, 1, 1, 2, 0),
    mod = c(1, 15, 1, 21, 3, 27, 5),
    no_snp = c(0, 2, 22, 2, 29, 4, 29),
    snp = c(0, 0, 0, 0, 0, 0, 0),
    coverage = c(1, 18, 24, 24, 33, 33, 34),
    genotype = c("C/C", "G/G", "C/C", "G/G", "C/C", "G/G", "C/C"),
    gt_p_score = c(0, 99, 99, 99, 99, 99, 99),
    gt_conf_score = c(0, 6, 66, 6, 87, 12, 87)
  )

  colnames(df)[1] <- "#chr"


  # Temp output file
  output_file <- tempfile(fileext = ".tsv")

  # Run function 
  save_as_methylkit(df, output_file)

  # Read and check result
  result <- readr::read_tsv(output_file, col_types = readr::cols())

  expect_equal(colnames(result), c("chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"))
  expect_equal(result$chrBase, paste0("J02459.1:", c(3,4,6,7,12,13,14)))
  expect_equal(result$chr, rep("J02459.1", 7))
  expect_equal(result$base, c(3,4,6,7,12,13,14))
  expect_equal(result$strand, c("F","R","F","R","F","R","F"))
  expect_equal(result$coverage, c(1,18,24,24,33,33,34))
  
  # freqC = 100 * mod / (mod + unmod)
  expected_freqC <- 100 * df$mod / (df$mod + df$unmod)
  expected_freqT <- 100 * df$unmod / (df$mod + df$unmod)

  expect_equal(result$freqC, expected_freqC)
  expect_equal(result$freqT, expected_freqT)

  unlink(output_file)
})
