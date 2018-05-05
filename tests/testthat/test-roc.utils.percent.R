library(pROC)
data(aSAH)

context("roc.utils.percent")

expect_equal_ignore_call <- function(x, y, ...) {
	x$call <- NULL
	y$call <- NULL
	attr(x$auc, "roc")$call <- NULL
	attr(y$auc, "roc")$call <- NULL
	attr(attr(x$ci, "auc"), "roc")$call <- NULL
	attr(attr(y$ci, "auc"), "roc")$call <- NULL
	expect_equal(x, y, ...)
}

test_that("roc.utils.topercent works on full AUC", {
	expect_equal_ignore_call(pROC:::roc.utils.topercent.roc(r.wfns), r.wfns.percent)
})

test_that("roc.utils.unpercent works on full AUC", {
	expect_equal_ignore_call(pROC:::roc.utils.unpercent.roc(r.wfns.percent), r.wfns)
})

test_that("roc.utils.topercent works on partial AUC", {
	expect_equal_ignore_call(pROC:::roc.utils.topercent.roc(r.wfns.partial1), r.wfns.percent.partial1)
})

test_that("roc.utils.unpercent works on partial AUC", {
	expect_equal_ignore_call(pROC:::roc.utils.unpercent.roc(r.wfns.percent.partial1), r.wfns.partial1)
})	

test_that("roc.utils.topercent works with CI", {
	r <- roc(aSAH$outcome, aSAH$wfns, ci=TRUE)
	r.percent <- roc(aSAH$outcome, aSAH$wfns, ci=TRUE, percent = TRUE)
	expect_equal_ignore_call(pROC:::roc.utils.topercent.roc(r), r.percent)
})

test_that("roc.utils.unpercent works with CI", {
	r <- roc(aSAH$outcome, aSAH$wfns, ci=TRUE)
	r.percent <- roc(aSAH$outcome, aSAH$wfns, ci=TRUE, percent = TRUE)
	expect_equal_ignore_call(pROC:::roc.utils.unpercent.roc(r.percent), r)
})

test_that("roc.utils.topercent works without AUC", {
	r <- roc(aSAH$outcome, aSAH$wfns, auc=FALSE)
	r.percent <- roc(aSAH$outcome, aSAH$wfns, auc=FALSE, percent = TRUE)
	expect_equal_ignore_call(pROC:::roc.utils.topercent.roc(r), r.percent)
})

test_that("roc.utils.unpercent works without AUC", {
	r <- roc(aSAH$outcome, aSAH$wfns, auc=FALSE)
	r.percent <- roc(aSAH$outcome, aSAH$wfns, auc=FALSE, percent = TRUE)
	expect_equal_ignore_call(pROC:::roc.utils.unpercent.roc(r.percent), r)
})