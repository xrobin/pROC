library(pROC)
data(aSAH)

context("roc.utils")

test_that("roc.utils.thr.idx finds correc thresholds with direction=<", {
	obtained <- pROC:::roc.utils.thr.idx(r.s100b, c(-Inf, 0.205, 0.055, Inf))
	expect_equal(obtained, c(1, 18, 4, 51))
})

test_that("roc.utils.thr.idx finds correc thresholds with direction=>", {
	obtained <- pROC:::roc.utils.thr.idx(r.s100b, c(Inf, -Inf, 0.05, 0.055, 0.52, 0.205))
	expect_equal(obtained, c(51, 1, 3, 4, 40, 18))
})

