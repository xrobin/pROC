library(pROC)
data(aSAH)

context("coords")
source(system.file("extdata", "test-coords-expected.R", package="pROC"))


test_that("coords with thresholds works", {
	obtained <- coords(r.s100b, "all", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords)
})

test_that("coords with percent works", {
	obtained <- coords(r.s100b.percent, "all", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	# Adjust for percent
	obtained[c("specificity", "sensitivity", "npv", "ppv", "1-specificity", "1-sensitivity", "1-npv", "1-ppv"),] <- 	obtained[c("specificity", "sensitivity", "npv", "ppv", "1-specificity", "1-sensitivity", "1-npv", "1-ppv"),] / 100
	expect_equal(obtained, expected.coords)
})