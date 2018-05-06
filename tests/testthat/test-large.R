library(pROC)

context("large data sets")

test_that("roc can deal with 1E5 data points", {
  response <- rbinom(1E5, 1, .5)
  predictor <- rnorm(1E5)
  roc(response, predictor)
})

test_that("roc can deal with 1E6 data points", {
  response <- rbinom(1E6, 1, .5)
  predictor <- rnorm(1E6)
  roc(response, predictor)
})

test_that("roc can deal with 1E7 data points", {
  response <- rbinom(1E7, 1, .5)
  predictor <- rnorm(1E7)
  roc(response, predictor)
})
