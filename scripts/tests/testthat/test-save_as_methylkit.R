options(run.main=FALSE)

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

  input_file <- tempfile(fileext = ".tsv")
  # write test data to file
  write.table(df, input_file, sep="\t", row.names = FALSE, quote=FALSE)

  # Temp output file
  output_file <- tempfile(fileext = ".tsv")

  system(paste("../../rastair_call_to_methylkit.sh", input_file, ">", output_file))

  # Read and check result
  result <- read.csv(output_file, sep="\t", quote="", header=TRUE)

  unlink(input_file)
  unlink(output_file)

  expect_equal(colnames(result), c("chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"))
  expect_equal(result$chrBase, paste0("J02459.1:", c(4, 5, 7, 8, 13, 14, 15)))
  expect_equal(result$chr, rep("J02459.1", 7))
  expect_equal(result$base, c(4, 5, 7, 8, 13, 14, 15))
  expect_equal(result$strand, c("F","R","F","R","F","R","F"))
  expect_equal(result$coverage, df$unmod + df$mod)

  # freqC = 100 * mod / (mod + unmod)
  expected_freqC <- round(100 * df$mod / (df$mod + df$unmod), 2)
  expected_freqT <- round(100 * df$unmod / (df$mod + df$unmod), 2)

  expect_equal(result$freqC, expected_freqC)
  expect_equal(result$freqT, expected_freqT)
})
