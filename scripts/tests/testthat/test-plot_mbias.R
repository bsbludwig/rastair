options(run.main=FALSE)
source("../../plot_mbias.R")

generate_test_table <- function() {
  mbias_table <- data.frame(X.pos = rep(1:50, 4), type=as.factor(c(rep("OT/1", 50), rep("OT/2", 50),rep("OB/1", 50),rep("OB/2", 50))), unmod=rbinom(50*4, 500, 0.05), mod=rbinom(50*4, 500, 0.95))
  mbias_table <- mbias_table[order(mbias_table$type),]
  # Add in some bias in R2s
  mbias_table[mbias_table$type=="OT/2" & mbias_table$X.pos > 45, "unmod"] <- rbinom(5, 500, 0.3)
  mbias_table[mbias_table$type=="OT/2" & mbias_table$X.pos > 45, "mod"] <- rbinom(5, 500, 0.7)

  mbias_table[mbias_table$type=="OB/2" & mbias_table$X.pos < 7, "unmod"] <- rbinom(6, 500, 0.35)
  mbias_table[mbias_table$type=="OB/2" & mbias_table$X.pos < 7, "mod"] <- rbinom(6, 500, 0.65)

  mbias_table$beta <- round(mbias_table$mod/(mbias_table$mod+mbias_table$unmod), 5)
  return(mbias_table)
}

test_that("post processing works", {
  test_table <- generate_test_table()
  expect_equal(dim(test_table)[2], 5)
  ob_1_front_pre <- mean(test_table[test_table$type=="OB/1" & test_table$X.pos < 7,"beta"])
  test_table <- post_process_mbias_table(test_table)
  expect_equal(dim(test_table)[2], 7)
  expect_contains(names(test_table), c("strand", "read_pair"))
  ob_1_front_post <- mean(test_table[test_table$type=="OB/1" & test_table$X.pos < 7,"beta"])
  expect_false(ob_1_front_pre==ob_1_front_post)
  ob_1_end_post <- mean(test_table[test_table$type=="OB/1" & test_table$X.pos > 44,"beta"])
  expect_equal(ob_1_front_pre, ob_1_end_post)
})

test_that("CI calculation works", {
  test_table <- generate_test_table()
  test_table <- post_process_mbias_table(test_table)
  test_table <- add_CI(test_table)
  expect_contains(names(test_table), c("upper_ci", "lower_ci"))

  expect_true(all(test_table$upper_ci >= test_table$lower_ci))
  expect_true(all(test_table$upper_ci >= test_table$beta & test_table$beta >= test_table$lower_ci))
  expect_true(all(test_table$upper_ci >= 0 & test_table$upper_ci <= 1))
  expect_true(all(test_table$lower_ci >= 0 & test_table$lower_ci <= 1))
})

test_that("cutoff calculation works", {
  test_table <- generate_test_table()
  test_table <- post_process_mbias_table(test_table)
  test_table <- add_CI(test_table)

  cutoffs <- rbind(data.frame(strand="OT", read_pair="1", as.list(calculate_cutoff(subset(test_table, strand=="OT" & read_pair == "1")))),
                   data.frame(strand="OT", read_pair="2", as.list(calculate_cutoff(subset(test_table, strand=="OT" & read_pair == "2")))),
                   data.frame(strand="OB", read_pair="1", as.list(calculate_cutoff(subset(test_table, strand=="OB" & read_pair == "1")))),
                   data.frame(strand="OB", read_pair="2", as.list(calculate_cutoff(subset(test_table, strand=="OB" & read_pair == "2")))))
  cutoffs$strand <- as.factor(cutoffs$strand)
  cutoffs$read_pair <- as.factor(cutoffs$read_pair)

  expect_gte(subset(cutoffs, strand=="OT" & read_pair == "2")$left, 5)
  expect_gte(subset(cutoffs, strand=="OB" & read_pair == "2")$left, 6)
  expect_equal(subset(cutoffs, strand=="OB" & read_pair == "1")$left, 0)
  expect_equal(subset(cutoffs, strand=="OT" & read_pair == "1")$left, 0)
  expect_equal(subset(cutoffs, strand=="OB" & read_pair == "1")$right, 0)
  expect_equal(subset(cutoffs, strand=="OT" & read_pair == "1")$right, 0)
  expect_equal(subset(cutoffs, strand=="OB" & read_pair == "2")$right, 0)
  expect_equal(subset(cutoffs, strand=="OT" & read_pair == "2")$right, 0)
})
