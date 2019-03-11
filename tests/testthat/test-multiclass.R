library(pROC)
data(aSAH)

context("multiclass-roc")

test_that("univariate multiclass roc/auc works", {
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b)
	expect_equal(class(uv.mr), "multiclass.roc")
	expect_equal(length(uv.mr$rocs), 6)
	expect_equal(as.numeric(auc(uv.mr)), 0.6539999352)
	expect_false(uv.mr$percent)
	expect_false(attributes(uv.mr$auc)$partial.auc)
	expect_false(attributes(uv.mr$auc)$percent)
})

test_that("univariate multiclass roc works with formula", {
	uv.mr <- multiclass.roc(gos6 ~ s100b, aSAH)
	expect_equal(as.numeric(auc(uv.mr)), 0.6539999352)
	uv.mr <- multiclass.roc(aSAH$gos6 ~ aSAH$s100b)
	expect_equal(as.numeric(auc(uv.mr)), 0.6539999352)
	
	uv.mr <- multiclass.roc(gos6 ~ s100b, aSAH, subset = (gender == "Female"))
	expect_equal(length(uv.mr$response), sum(aSAH$gender == "Female"))
	expect_error(multiclass.roc(gos6 ~ s100b, aSAH, weights=age))
})

test_that("univariate multiclass roc/auc works with percent=TRUE", {
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b, percent=TRUE)
	uv.ma <- auc(uv.mr)
	expect_equal(as.numeric(uv.ma), 65.39999352)
	expect_equal(as.numeric(uv.mr$auc), 65.39999352)
	expect_true(uv.mr$percent)
	expect_true(attributes(uv.mr$auc)$percent)
	expect_true(attributes(uv.ma)$percent)
	expect_false(attributes(uv.mr$auc)$partial.auc)
	expect_false(attributes(uv.ma)$partial.auc)
})

test_that("univariate multiclass roc/auc works with partial.auc", {
	pauc.spec <- c(1, .9)
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b, partial.auc=pauc.spec)
	uv.ma <- auc(uv.mr, partial.auc=pauc.spec)
	expect_equal(as.numeric(uv.mr$auc), 0.0116176879)
	expect_equal(as.numeric(uv.ma), 0.0116176879)
	expect_equal(attributes(uv.mr$auc)$partial.auc, pauc.spec)
	expect_equal(attributes(uv.ma)$partial.auc, pauc.spec)
	# Calling AUC without partial.auc gives a full AUC, even if ROC was called with it
	expect_equal(as.numeric(auc(uv.mr)), 0.6539999352)

	# SE
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b, partial.auc=pauc.spec, partial.auc.focus="se")
	uv.ma <- auc(uv.mr, partial.auc=pauc.spec, partial.auc.focus="se")
	expect_equal(as.numeric(uv.mr$auc), 0.02513286)
	expect_equal(as.numeric(uv.ma), 0.02513286)
	expect_equal(attributes(uv.mr$auc)$partial.auc.focus, "sensitivity")
	expect_equal(attributes(uv.ma)$partial.auc.focus, "sensitivity")
	
})

test_that("univariate multiclass roc/auc works with directions", {
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b, direction="auto")
	expect_equal(sapply(uv.mr$rocs, "[[", "direction"), c("<", ">", ">", ">", ">", ">"))
	expect_equal(as.numeric(uv.mr$auc), 0.6539999352)
	
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b, direction="<")
	expect_equal(sapply(uv.mr$rocs, "[[", "direction"), rep("<", 6))
	expect_equal(as.numeric(uv.mr$auc), 0.3487473175)
	
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b, direction=">")
	expect_equal(sapply(uv.mr$rocs, "[[", "direction"), rep(">", 6))
	expect_equal(as.numeric(uv.mr$auc), 0.6512526825)
})

test_that("univariate multiclass handles missing levels", {
	expect_warning(multiclass.roc(aSAH$gos6, aSAH$s100b), "response level")
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
	expect_equal(length(mv.mr$rocs), 3)
	expect_equal(unname(sapply(mv.mr$rocs, length)), rep(2, 3))
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



test_that("multivariate multiclass with formula works", {
	n <- c(100, 80, 150)
	responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
	# construct prediction matrix: one column per class
	set.seed(42)
	
	# Imperfect separation
	preds <- lapply(n, function(x) runif(x, 0.4, 0.6))
	predictor <- as.matrix(data.frame(
		"X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.7)),
		"X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0.2, 0.8)),
		"X3" = c(runif(n[1] + n[2], 0.3, 0.7), preds[[3]])
	))
	test_data <- cbind(as.data.frame(predictor), "response" = responses)
	mv.mr <- multiclass.roc(responses, predictor)
	mv.mr.f <- multiclass.roc(response ~ ., data=test_data)
	expect_equal(as.numeric(auc(mv.mr.f)), as.numeric(auc(mv.mr)))
	
	mv.mr.f <- multiclass.roc(response ~ X1 + X2 + X3, data=test_data)
	expect_equal(as.numeric(auc(mv.mr)), as.numeric(auc(mv.mr)))
	
	subset <- rbinom(sum(n), 1, .5) == 1
	mv.mr.f <- multiclass.roc(response ~ X1 + X2 + X3, data=test_data, subset = subset)
	expect_equal(length(mv.mr.f$response), sum(subset))
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


test_that("multivariate multiclass roc/auc works with partial.auc", {
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


test_that("multivariate multiclass roc/auc works with direction", {
	n <- c(100, 80, 150)
	responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
	# construct prediction matrix: one column per class
	set.seed(42)
	
	# Perfect separation
	preds <- lapply(n, function(x) runif(x, 0.8, 1))
	predictor <- as.matrix(data.frame("X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.3)),
									  "X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0, 0.2)),
									  "X3" = c(runif(n[1] + n[2], 0, 0.5), preds[[3]])))
	
	expect_error(multiclass.roc(responses, predictor, direction = "auto"))
	
	mr.mv.1 <- multiclass.roc(responses, predictor, direction = "<")
	mr.mv.2 <- multiclass.roc(responses, predictor, direction = ">")
	expect_equal(as.numeric(mr.mv.1$auc), 0)
	expect_equal(as.numeric(mr.mv.2$auc), 1)
	
	for (i in 1:3) {
		for (j in 1:2) {
			expect_equal(mr.mv.1$rocs[[i]][[j]]$direction, "<")
			expect_equal(mr.mv.2$rocs[[i]][[j]]$direction, ">")
		}
	}
})

test_that("multivariate behavior with missing levels/columns", {
	n <- c(10, 10, 10)
	responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
	# construct prediction matrix: one column per class
	set.seed(42)
	
	# Perfect separation
	preds <- lapply(n, function(x) runif(x, 0.8, 1))
	predictor <- as.matrix(data.frame("X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.3)),
									  "X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0, 0.2)),
									  "X3" = c(runif(n[1] + n[2], 0, 0.5), preds[[3]])))
	# Wrong number of predictor rows
	expect_error(multiclass.roc(responses[1:20], predictor), "agree")
	
	# Column in predictor not in response warns:
	expect_warning(multiclass.roc(as.character(responses[1:20]), predictor[1:20,]), "X3")
	
	# Level with no obervation warns:
	expect_warning(multiclass.roc(responses[1:20], predictor[1:20,1:2]), "X3")
	
	# Removed both level and column should be silent
	expect_silent(multiclass.roc(as.character(responses[1:20]), predictor[1:20,1:2]))
	
	# Single column is an error
	expect_error(multiclass.roc(responses, predictor[,1, drop=F]))
	
	# Wrong column names
	pr2 <- predictor
	colnames(pr2) <- c("Y1", "Y2", "Y3")
	expect_error(multiclass.roc(responses, pr2))
	colnames(pr2) <- c("Y1", "X2", "X3")
	expect_warning(multiclass.roc(as.character(responses[11:30]), pr2[11:30,]), "Y1")
})

test_that("Invalid CI functions fail cleanly", {
	uv.mr <- multiclass.roc(aSAH$gos6, aSAH$s100b)
	expect_error(ci.se(uv.mr), "not available for multiclass ROC curves")
	expect_error(ci.se(uv.mr$auc), "not available for multiclass ROC curves")
	expect_error(ci.sp(uv.mr), "not available for multiclass ROC curves")
	expect_error(ci.sp(uv.mr$auc), "not available for multiclass ROC curves")
	expect_error(ci.coords(uv.mr), "not available for multiclass ROC curves")
	expect_error(ci.coords(uv.mr$auc), "not available for multiclass ROC curves")
	expect_error(ci.thresholds(uv.mr), "not available for multiclass ROC curves")
	expect_error(ci.thresholds(uv.mr$auc), "not available for multiclass ROC curves")
})


