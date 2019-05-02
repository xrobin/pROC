library(pROC)
data(aSAH)

context("coords")

test_that("coords with thresholds works", {
	obtained <- coords(r.s100b, "all", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords)
})


test_that("coords with percent works", {
	obtained.percent <- coords(r.s100b.percent, "all", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	# Adjust for percent
	obtained.percent[c("specificity", "sensitivity", "accuracy", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"),] <- obtained.percent[c("specificity", "sensitivity", "accuracy", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"),] / 100
	expect_equal(obtained.percent, expected.coords)
})


test_that("coords with local maximas thresholds works", {
	obtained <- coords(r.s100b, "local maximas", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expected.thresholds <- c(-Inf, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.135, 0.155, 0.205, 0.245, 0.29, 0.325, 0.345, 0.395, 0.435, 0.475, 0.485, 0.51)
	expect_equivalent(obtained, expected.coords[,expected.coords["threshold",] %in% expected.thresholds])
})


test_that("coords with best threshold works", {
	obtained <- coords(r.s100b, "best", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords[,expected.coords["threshold",] == 0.205])
})


test_that("coords with arbitrary thresholds works", {
	obtained <- coords(r.s100b, c(0.205, 0.055), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equivalent(obtained, expected.coords[, c(18, 4)])
})

test_that("coords with arbitrary thresholds at exact data point works", {
	expect_equal(sum(aSAH$s100b == 0.05),  3)
	expect_equal(sum(aSAH$s100b == 0.52),  1)
	obtained <- coords(r.s100b, c(0.05, 0.52), input = "threshold", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equivalent(obtained, expected.coords[-1, c(3, 40)])
})

test_that("coords with arbitrary thresholds works with direction=>", {
	obtained <- coords(r.100b.reversed, c(0.05, 0.055, 0.205, 0.52), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equivalent(obtained, expected.coords.reverse)
})


test_that("coords with single arbitrary threshold works", {
	obtained <- coords(r.s100b, c(0.205), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords[, c(18), drop=T])
})


test_that("coords with arbitrary thresholds at exact data point works", {
	expect_equal(sum(aSAH$s100b == 0.05),  3)
	expect_equal(sum(aSAH$s100b == 0.52),  1)
	obtained <- coords(r.s100b, c(0.05), input = "threshold", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords[-1, 3])
	obtained <- coords(r.s100b, c(0.52), input = "threshold", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords[-1, 40])
})


test_that("coords with arbitrary thresholds works with direction=>", {
	obtained <- coords(r.100b.reversed, c(0.05), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords.reverse[, 1])
	obtained <- coords(r.100b.reversed, c(0.055), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords.reverse[, 2])
	obtained <- coords(r.100b.reversed, c(0.205), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords.reverse[, 3])
	obtained <- coords(r.100b.reversed, c(0.52), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"))
	expect_equal(obtained, expected.coords.reverse[, 4])
})


test_that("coords with sensitivity works", {
	obtained <- coords(r.s100b, seq(0, 1, .1), input = "sensitivity", ret = c("threshold", "specificity", "sensitivity"))
	expect_equal(unname(obtained["threshold",]), c(Inf, rep(NA, 9), -Inf))
	expect_equal(unname(obtained["sensitivity",]), seq(0, 1, .1))
	expect_equal(unname(obtained["specificity",]), c(1, 1, 1, 0.972222222222222, 0.888888888888889, 0.833333333333333, 0.805555555555556, 0.56875, 0.447222222222222, 0.230555555555556, 0))
})


test_that("coords with sensitivity works with percent", {
	obtained <- coords(r.s100b.percent, seq(0, 100, 10), input = "sensitivity", ret = c("threshold", "specificity", "sensitivity"))
	expect_equal(unname(obtained["threshold",]), c(Inf, rep(NA, 9), -Inf))
	expect_equal(unname(obtained["sensitivity",]), seq(0, 100, 10))
	expect_equal(unname(obtained["specificity",]), c(1, 1, 1, 0.972222222222222, 0.888888888888889, 0.833333333333333, 0.805555555555556, 0.56875, 0.447222222222222, 0.230555555555556, 0) * 100)
})


test_that("coords with specificity works", {
	obtained <- coords(r.s100b, seq(0, 1, .1), input = "specificity", ret = c("threshold", "specificity", "sensitivity"))
	expect_equal(unname(obtained["threshold",]), c(-Inf, rep(NA, 9), 0.51))
	expect_equal(unname(obtained["specificity",]), seq(0, 1, .1))
	expect_equal(unname(obtained["sensitivity",]), c(1, 0.975609756097561, 0.921951219512195, 0.879674796747967, 0.823693379790941, 0.774390243902439, 0.675609756097561, 0.655284552845528, 0.634146341463415, 0.390243902439024, 0.292682926829268))
})


test_that("coords with specificity works with percent", {
	obtained <- coords(r.s100b.percent, seq(0, 100, 10), input = "specificity", ret = c("threshold", "specificity", "sensitivity"))
	expect_equal(unname(obtained["threshold",]), c(-Inf, rep(NA, 9), 0.51))
	expect_equal(unname(obtained["specificity",]), seq(0, 100, 10))
	expect_equal(unname(obtained["sensitivity",]), c(1, 0.975609756097561, 0.921951219512195, 0.879674796747967, 0.823693379790941, 0.774390243902439, 0.675609756097561, 0.655284552845528, 0.634146341463415, 0.390243902439024, 0.292682926829268) * 100)
})


test_that("coords with specificity works with as.list", {
	obtained <- coords(r.s100b.percent, "best", ret = c("threshold", "specificity", "accuracy"), as.list = TRUE)
	expect_equal(obtained, list(
		threshold = 0.205,
		specificity = expected.coords["specificity", 18] * 100,
		accuracy = expected.coords["accuracy", 18] * 100
	))
})

test_that("coords with specificity works with as.list and drop=FALSE", {
	obtained <- coords(r.s100b.percent, "best", 
					   ret = c("threshold", "specificity", "accuracy"), 
					   as.list = TRUE, drop = FALSE)
	expect_equal(obtained$best, list(
		threshold = 0.205,
		specificity = expected.coords["specificity", 18] * 100,
		accuracy = expected.coords["accuracy", 18] * 100
	))
})


test_that("coords with specificity works with as.list and several thresholds", {
	obtained <- coords(r.s100b.percent, c(0.205, 0.51), 
					   ret = c("threshold", "specificity", "accuracy"), 
					   as.list = TRUE, drop = FALSE)
	expect_equal(names(obtained), c("0.205", "0.51"))
	expect_equal(obtained[[1]], list(
		threshold = 0.205,
		specificity = expected.coords["specificity", 18] * 100,
		accuracy = expected.coords["accuracy", 18] * 100
	))
	expect_equal(obtained[[2]], list(
		threshold = 0.51,
		specificity = expected.coords["specificity", 40] * 100,
		accuracy = expected.coords["accuracy", 40] * 100
	))
})


test_that("drop works", {
	skip("drop doesn't drop over x - doc unclear/inconsistent need to be fixed")
	# First make sure we get matrices with drop = FALSE
	expect_is(coords(r.s100b, 0.51, input = "threshold", ret = c("sensitivity", "specificity"), drop = FALSE), "matrix")
	expect_is(coords(r.s100b, 0.51, input = "threshold", ret = "specificity", drop = FALSE), "matrix")
	expect_is(coords(r.s100b, "local maximas", input = "threshold", ret = "specificity", drop = FALSE), "matrix")	
	expect_is(coords(r.s100b, c(0.51, 0.2), input = "threshold", ret = "specificity", drop = FALSE), "matrix")
	# Look for numeric
	expect_is(coords(r.s100b, 0.51, input = "threshold", ret = c("sensitivity", "specificity"), drop = TRUE), "numeric")
	expect_is(coords(r.s100b, 0.51, input = "threshold", ret = "specificity", drop = TRUE), "numeric")
	expect_is(coords(r.s100b, "local maximas", input = "threshold", ret = "specificity", drop = TRUE), "numeric")	
	expect_is(coords(r.s100b, c(0.51, 0.2), input = "threshold", ret = "specificity", drop = TRUE), "numeric")
})


test_that("coords returns the correct basic values ", {
	obtained <- coords(r.s100b, 0.205, 
					   ret = c("t", "tp", "fp", "tn", "fn",
					   		   "sp", "se", "acc",
					   		   "npv", "ppv", "precision", "recall",
					   		   "tpr", "fpr", "tnr", "fnr", "fdr"))
	
	obtained.percent <- coords(r.s100b.percent, 0.205, 
					   ret = c("t", "tp", "fp", "tn", "fn",
					   		"sp", "se", "acc",
					   		"npv", "ppv", "precision", "recall",
					   		"tpr", "fpr", "tnr", "fnr", "fdr"))
	
	# We assume the following values:
	# tp fp tn fn N
	# 26 14 58 15 113
	
	expected <- c(
		threshold = 0.205,
		tp = 26,
		fp = 14,
		tn = 58,
		fn = 15,
		specificity = 58 / (58 + 14),
		sensitivity = 26 / (26 + 15),
		accuracy = (26 + 58) / 113,
		npv = 58 / (58 + 15),
		ppv = 26 / (26 + 14),
		precision = 26 / (26 + 14),
		recall = 26 / (26 + 15),
		tpr = 26 / (26 + 15),
		fpr = 1 - (58 / (58 + 14)),
		tnr = 58 / (58 + 14),
		fnr = 1 - (26 / (26 + 15)),
		fdr = 14 / (26 + 14)
	)
	
	expect_equal(obtained, expected)
	expect_equal(obtained.percent[1:5], expected[1:5])
	expect_equal(obtained.percent[6:17], expected[6:17]*100)
})