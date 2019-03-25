context("Check inputs")
library(invctr)

test_that("RNGmeasures", {
  expect_output(RNGmeasures(y = round(runif(100,1,9)), minScale = 1, maxScale = 9), "1 variables")
})
