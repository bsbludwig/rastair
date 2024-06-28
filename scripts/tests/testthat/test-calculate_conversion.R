options(run.main=FALSE)
source("../../calculate_conversion.R")

generate_test_table <- function() {
  methylation <- data.frame(X.chr = rep("chr1", 100), start = 0:99, end=1:100, name=NA, strand=c("+","-"), unmod=rbinom(100, 200, 0.03), mod=rbinom(100, 200, 0.97), snp=0, nosnp=round(rnorm(100, 200, 50)))
  methylation[sample.int(100, 5),c("mod", "unmod")] <- 0
  methylation$beta_est <- methylation$mod/(methylation$mod+methylation$unmod)
  methylation[is.nan(methylation$beta_est),"beta_est"] <- NA
  return(methylation)
}

test_that("beta calculation works", {
  test_table <- generate_test_table()
  summary_table <- calculate_summary(test_table)
  expect_equal(length(summary_table), 9)
  expect_true(summary_table["mean"] == mean(test_table$beta_est, na.rm=TRUE))
  expect_true(summary_table["median"] == median(test_table$beta_est, na.rm=TRUE))
  expect_true(summary_table["sum_mean"] == sum(test_table$mod)/sum(test_table$mod+test_table$unmod))

  # try without NAs
  test_table <- subset(test_table, !is.na(beta_est))
  summary_table <- calculate_summary(test_table)
  expect_equal(length(summary_table), 9)
  expect_true(summary_table["mean"] == mean(test_table$beta_est, na.rm=TRUE))
  expect_true(summary_table["median"] == median(test_table$beta_est, na.rm=TRUE))
  expect_true(summary_table["sum_mean"] == sum(test_table$mod)/sum(test_table$mod+test_table$unmod))
})

test_that("context annotation works", {
  test_table <- read.delim("../data/lambda_calls.bed.gz", stringsAsFactors=TRUE, na.strings = c("NA","NaN",".",""))
  annotated_table <- annotate_with_context(test_table, "../../../test_data/test.fasta.gz")
  expect_equal(dim(annotated_table)[1], 3052)
  expect_gte(dim(annotated_table)[2], 5)
  expect_true(all(substr(annotated_table$sequence, 2, 3) == "CG"))
})

test_that("plot generation works", {
  test_table <- read.delim("../data/lambda_calls.bed.gz", stringsAsFactors=TRUE, na.strings = c("NA","NaN",".",""))
  annotated_table <- annotate_with_context(test_table, "../../../test_data/test.fasta.gz")

  tf <- tempfile(, fileext = ".png")
  plot_context_stats(annotated_table, c("-o", tf))
  expect_true(file.exists(tf))
  expect_gt(file.size(tf), 0)
})
