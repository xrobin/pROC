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

test_that("roc.utils.calc.coords works", {
	obtained <- pROC:::roc.utils.calc.coords(r.s100b, -1:-4, c(1, .5, .1, 0), c(0, .5, .9, 1), c(12, .9))
	expect_equal(obtained, expected_roc.utils.calc.coords)
})

test_that("roc.utils.calc.coords works with percent", {
	obtained <- pROC:::roc.utils.calc.coords(r.s100b.percent, -1:-4, c(100, 50, 10, 0), c(0, 50, 90, 100), c(12, .9))
	expect_equal(obtained, expected_roc.utils.calc.coords.percent)
})

