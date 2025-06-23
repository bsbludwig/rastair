options(run.main=FALSE)
source("../../parse_mbias.R")

test_that("parse_mbias correctly transposes and reshapes example input", {

  # Create test input file with your exact example
  test_input <- tempfile(fileext = ".txt")
  test_output <- tempfile(fileext = ".csv")

  example_data <- "strand\tread_pair\tleft\tright
                    OT\t1\t2\t0
                    OT\t2\t7\t0
                    OB\t1\t5\t0
                    OB\t2\t5\t0"
  
  writeLines(example_data, test_input)

  # Run the function
  parse_mbias(test_input, test_output)

  # Read and check output
  result <- read.csv(test_output, header = FALSE)

  # Should be 2 rows, 4 columns
  expect_equal(dim(result), c(2, 4))

  # Check values row by row
  expect_equal(as.numeric(result[1, ]), c(2, 0, 7, 0))
  expect_equal(as.numeric(result[2, ]), c(5, 0, 5, 0))

  unlink(c(test_input, test_output))
})
