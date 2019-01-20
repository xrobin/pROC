library(pROC)
data(aSAH)

context("multiclass-roc")

test_that("univariate multiclass roc/auc works", {
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b)
	expect_equal(class(uv.mr), "multiclass.roc")
	expect_equal(as.numeric(auc(uv.mr)), 0.6539999352)
	expect_false(uv.mr$percent)
	expect_false(attributes(uv.mr$auc)$partial.auc)
	expect_false(attributes(uv.mr$auc)$percent)
})

test_that("univariate multiclass roc/auc works with percent=TRUE", {
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b, percent=TRUE)
	uv.ma <- auc(uv.mr)
	expect_equal(as.numeric(auc(uv.mr)), 65.39999352)
	expect_true(uv.mr$percent)
	expect_false(attributes(uv.mr$auc)$partial.auc)
	expect_true(attributes(uv.mr$auc)$percent)
})

test_that("univariate multiclass roc/auc works with partial.auc", {
	pauc.spec <- c(1, .9)
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b, partial.auc=pauc.spec)
	uv.ma <- auc(uv.mr)
	expect_equal(as.numeric(uv.mr$auc), 0.0116176879)
	expect_equal(as.numeric(auc(uv.mr)), 0.6539999352)
	expect_equal(as.numeric(auc(uv.mr, partial.auc=pauc.spec)), 0.0116176879)
	expect_false(uv.mr$percent)
	expect_equal(attributes(uv.mr$auc)$partial.auc, pauc.spec)
	expect_false(attributes(uv.mr$auc)$percent)
})


test_that("multivariate multiclass roc/auc works", {
	n <- c(100, 80, 150)
	responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
	# construct prediction matrix: one column per class
	set.seed(42)
	
	# Perfect separation
	preds <- lapply(n, function(x) runif(x, 0.8, 1))
	predictor <- as.matrix(data.frame("X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.3)),
									  "X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0, 0.2)),
									  "X3" = c(runif(n[1] + n[2], 0, 0.5), preds[[3]])))
	mv.mr <- multiclass.roc(responses, predictor)
	expect_equal(class(mv.mr), "mv.multiclass.roc")
	expect_equal(as.numeric(auc(mv.mr)), 1)
	
	# Imperfect separation
	preds <- lapply(n, function(x) runif(x, 0.4, 0.6))
	predictor <- as.matrix(data.frame(
		"X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.7)),
		"X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0.2, 0.8)),
		"X3" = c(runif(n[1] + n[2], 0.3, 0.7), preds[[3]])
	))
	mv.mr <- multiclass.roc(responses, predictor)
	expect_equal(class(mv.mr), "mv.multiclass.roc")
	expect_equal(as.numeric(auc(mv.mr)), 0.6480791667)
	expect_false(mv.mr$percent)
	expect_false(attributes(mv.mr$auc)$partial.auc)
	expect_false(attributes(mv.mr$auc)$percent)
})


test_that("multivariate multiclass roc/auc works with percent=TRUE", {
	n <- c(100, 80, 150)
	responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
	# construct prediction matrix: one column per class
	set.seed(42)
	
	# Perfect separation
	preds <- lapply(n, function(x) runif(x, 0.8, 1))
	predictor <- as.matrix(data.frame("X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.3)),
									  "X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0, 0.2)),
									  "X3" = c(runif(n[1] + n[2], 0, 0.5), preds[[3]])))
	mv.mr <- multiclass.roc(responses, predictor, percent=TRUE)
	expect_equal(as.numeric(auc(mv.mr)), 100)
	
	# Imperfect separation
	preds <- lapply(n, function(x) runif(x, 0.4, 0.6))
	predictor <- as.matrix(data.frame(
		"X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.7)),
		"X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0.2, 0.8)),
		"X3" = c(runif(n[1] + n[2], 0.3, 0.7), preds[[3]])
	))
	mv.mr <- multiclass.roc(responses, predictor, percent=TRUE)
	expect_equal(as.numeric(auc(mv.mr)), 64.80791667)
	expect_true(mv.mr$percent)
	expect_false(attributes(mv.mr$auc)$partial.auc)
	expect_true(attributes(mv.mr$auc)$percent)
})


test_that("univariate multiclass roc/auc works with partial.auc", {
	pauc.spec <- c(1, .9)
	n <- c(100, 80, 150)
	responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
	# construct prediction matrix: one column per class
	set.seed(42)
	
	# Perfect separation
	preds <- lapply(n, function(x) runif(x, 0.8, 1))
	predictor <- as.matrix(data.frame("X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.3)),
									  "X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0, 0.2)),
									  "X3" = c(runif(n[1] + n[2], 0, 0.5), preds[[3]])))
	mv.mr <- multiclass.roc(responses, predictor, partial.auc=c(1, .9))
	expect_equal(as.numeric(mv.mr$auc), .1)
	expect_equal(as.numeric(auc(mv.mr)), 1)
	expect_equal(as.numeric(auc(mv.mr, partial.auc=c(1, .9))), .1)
	
	# Imperfect separation
	preds <- lapply(n, function(x) runif(x, 0.4, 0.6))
	predictor <- as.matrix(data.frame(
		"X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.7)),
		"X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0.2, 0.8)),
		"X3" = c(runif(n[1] + n[2], 0.3, 0.7), preds[[3]])
	))
	mv.mr <- multiclass.roc(responses, predictor, partial.auc=c(1, .9))
	expect_equal(as.numeric(mv.mr$auc), 0.0529250000)
	expect_equal(as.numeric(auc(mv.mr)), 0.6480791667)
	expect_equal(as.numeric(auc(mv.mr, partial.auc=c(1, .9))), 0.0529250000)
	expect_false(mv.mr$percent)
	expect_equal(attributes(mv.mr$auc)$partial.auc, pauc.spec)
	expect_false(attributes(mv.mr$auc)$percent)
})
