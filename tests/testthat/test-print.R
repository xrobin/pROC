library(pROC)
data(aSAH)

context("print")

test_that("print.auc works", {
	expect_output(print(auc(r.wfns)), "^Area under the curve: 0.8237$")
	expect_output(print(auc(r.ndka.percent)), "^Area under the curve: 61.2%$")
	
	expect_output(print(r.ndka.partial$auc), "^Partial area under the curve \\(specificity 1-0.9\\): 0.0107$")
	expect_output(print(r.s100b.percent.partial1$auc), "^Partial area under the curve \\(specificity 99%-90%\\): 2.983%$")
	expect_output(print(r.s100b.partial2$auc), "^Partial area under the curve \\(sensitivity 0.99-0.9\\): 0.01376$")
})

test_that("print.roc works", {
	expect_known_output(print(r.wfns), "print_output/r.wfns")
	expect_known_output(print(r.ndka), "print_output/r.ndka")
	expect_known_output(print(r.s100b), "print_output/r.s100b")
	expect_known_output(print(r.wfns.percent), "print_output/r.wfns.percent")
	expect_known_output(print(r.ndka.percent), "print_output/r.ndka.percent")
	expect_known_output(print(r.s100b.percent), "print_output/r.s100b.percent")
	expect_known_output(print(r.wfns.partial1), "print_output/r.wfns.partial1")
	expect_known_output(print(r.ndka.partial1), "print_output/r.ndka.partial1")
	expect_known_output(print(r.s100b.partial1), "print_output/r.s100b.partial1")
	expect_known_output(print(r.wfns.percent.partial1), "print_output/r.wfns.percent.partial1")
	expect_known_output(print(r.ndka.percent.partial1), "print_output/r.ndka.percent.partial1")
	expect_known_output(print(r.s100b.percent.partial1), "print_output/r.s100b.percent.partial1")
	expect_known_output(print(r.s100b.partial2), "print_output/r.s100b.partial2")
	expect_known_output(print(roc(outcome ~ ndka, aSAH)), "print_output/ndka_formula")
})

test_that("print.multiclass.roc works", {
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka)), "print_output/multiclass")
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, levels=c(3, 4, 5))), "print_output/multiclass_levels")
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, percent=TRUE)), "print_output/multiclass_percent")
})

test_that("print.multiclass.roc works", {
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka)), "print_output/multiclass")
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, levels=c(3, 4, 5))), "print_output/multiclass_levels")
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, percent=TRUE)), "print_output/multiclass_percent")
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, partial.auc=c(1, .9))), "print_output/multiclass_partial")
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, partial.auc=c(1, .9), partial.auc.focus="se")), "print_output/multiclass_partial_se")
	expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$wfns, partial.auc=c(1, .9), partial.auc.correct=TRUE)), "print_output/multiclass_partial_correct")
})

test_that("print.multiclass.roc multivariate works", {
	n <- c(100, 80, 150)
	responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
	set.seed(42)
	preds <- lapply(n, function(x) runif(x, 0.4, 0.6))
	predictor <- as.matrix(data.frame(
		"X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.7)),
		"X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0.2, 0.8)),
		"X3" = c(runif(n[1] + n[2], 0.3, 0.7), preds[[3]])
	))
	expect_known_output(print(multiclass.roc(responses, predictor)), "print_output/mv_multiclass")
	
	expect_known_output(print(multiclass.roc(responses, predictor, levels=c("X2", "X3"))), "print_output/mv_multiclass_levels")
	expect_known_output(print(multiclass.roc(responses, predictor, percent=TRUE)), "print_output/mv_multiclass_percent")
	expect_known_output(print(multiclass.roc(responses, predictor, partial.auc=c(1, .9))), "print_output/mv_multiclass_partial")
	expect_known_output(print(multiclass.roc(responses, predictor, partial.auc=c(1, .9), partial.auc.focus="se")), "print_output/mv_multiclass_partial_se")
	expect_known_output(print(multiclass.roc(responses, predictor, partial.auc=c(1, .9), partial.auc.correct=TRUE)), "print_output/mv_multiclass_partial_correct")
})

test_that("print works with a formula", {
	expect_known_output(print(roc(outcome ~ ndka, aSAH)), "print_output/r.ndka.formula")
	expect_known_output(print(multiclass.roc(gos6 ~ ndka, aSAH)), "print_output/mv_multiclass.ndka.formula")
})

test_that("print works without the auc", {
	expect_known_output(print(roc(outcome ~ ndka, aSAH, auc=FALSE)), "print_output/r.ndka.formula.no_auc")
})

test_that("print works with the CI", {
	skip_slow()
	if (R.version$minor >= "6.0") {
		RNGkind(sample.kind="Rounding")
	}
	set.seed(42) # For reproducible CI
	expect_known_output(print(roc(outcome ~ ndka, aSAH, ci=TRUE)), "print_output/r.ndka.formula.ci")
})

test_that("print.smooth.roc works", {
	expect_known_output(print(smooth(roc(aSAH$outcome, aSAH$ndka))), "print_output/smooth.ndka")
	expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH))), "print_output/smooth.s100b.formula")
	expect_known_output(print(smooth(roc(aSAH$outcome, aSAH$ndka))), "print_output/smooth.wfns")
	expect_known_output(print(smooth(roc(aSAH$outcome, aSAH$ndka), method="binormal")), "print_output/smooth.s100b.binormal")
	expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH), method="fitdistr")), "print_output/smooth.s100b.fitdistr")
	expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH), method="density")), "print_output/smooth.s100b.density")
	
	testthat::skip_if_not_installed("logcondens")
	expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH), method="logcondens")), "print_output/smooth.s100b.logcondens")
	expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH), method="logcondens.smooth")), "print_output/smooth.s100b.logcondens.smooth")
})

test_that("print works with ci.auc", {
	skip_slow()
	if (R.version$minor >= "6.0") {
		RNGkind(sample.kind="Rounding")
	}
	set.seed(42) # For reproducible CI
	expect_known_output(print(ci.auc(r.ndka, method = "bootstrap", boot.n = 3, progress = "none")), "print_output/r.ndka.ci.auc")
})

test_that("print works with ci.coords", {
	skip_slow()
	if (R.version$minor >= "6.0") {
		RNGkind(sample.kind="Rounding")
	}
	set.seed(42) # For reproducible CI
	expect_known_output(print(ci.coords(r.ndka, x = c(0.5, 0.2), boot.n = 3, progress = "none")), "print_output/r.ndka.ci.coords")
})

test_that("print works with ci.thresholds", {
	skip_slow()
	if (R.version$minor >= "6.0") {
		RNGkind(sample.kind="Rounding")
	}
	set.seed(42) # For reproducible CI
	expect_known_output(print(ci.thresholds(r.ndka, thresholds = c(0.5, 0.2), boot.n = 3, progress = "none")), "print_output/r.ndka.ci.thresholds")
})

test_that("print works with ci.se", {
	skip_slow()
	if (R.version$minor >= "6.0") {
		RNGkind(sample.kind="Rounding")
	}
	set.seed(42) # For reproducible CI
	expect_known_output(print(ci.se(r.ndka, boot.n = 3, progress = "none")), "print_output/r.ndka.ci.se")
})

test_that("print works with ci.sp", {
	skip_slow()
	if (R.version$minor >= "6.0") {
		RNGkind(sample.kind="Rounding")
	}
	set.seed(42) # For reproducible CI
	expect_known_output(print(ci.sp(r.ndka, boot.n = 3, progress = "none")), "print_output/r.ndka.ci.sp")
})

