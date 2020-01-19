library(pROC)
data(aSAH)

context("are.paired")

test_that("are.paired works", {
	# most basic example
	expect_true(are.paired(r.wfns, r.ndka))
	
	# Missing values shouldn't screw up
	aSAH.missing <- aSAH
	aSAH.missing$wfns[1:20] <- NA
	expect_true(are.paired(roc(aSAH.missing$outcome, aSAH.missing$wfns), roc(aSAH.missing$outcome, aSAH.missing$ndka)))
	# Also with different data.frames
	expect_true(are.paired(roc(aSAH.missing$outcome, aSAH.missing$wfns), r.ndka))
	
	# The following should fail though
	expect_false(are.paired(roc(aSAH$outcome[21:113], aSAH$wfns[21:113]), roc(aSAH$outcome, aSAH$ndka)))
	
	# Opposite levels should probably fail
	expect_false(are.paired(roc(aSAH$outcome, aSAH$wfns, levels = c("Good", "Poor")), roc(aSAH$outcome, aSAH$ndka, levels = c("Poor", "Good"))))
})

test_that("are.paired works with formula", {
	r.wfns.f <- roc(outcome ~ wfns, aSAH)
	r.ndka.f <- roc(outcome ~ ndka, aSAH)
	# most basic example
	expect_true(are.paired(r.wfns.f, r.ndka.f))
	
	# Missing values shouldn't screw up
	aSAH.missing <- aSAH
	aSAH.missing$wfns[1:20] <- NA
	expect_true(are.paired(roc(outcome ~ wfns, aSAH.missing), roc(outcome ~ ndka, aSAH.missing)))
	# Also with different data.frames
	expect_true(are.paired(roc(outcome ~ wfns, aSAH.missing), r.ndka.f))
	
	# The following should fail though
	expect_false(are.paired(roc(outcome ~ wfns, aSAH.missing[21:113,]), r.ndka))
	
	# Opposite levels should probably fail
	expect_false(are.paired(roc(outcome ~ wfns, aSAH, levels = c("Good", "Poor")), roc(outcome ~ ndka, aSAH, levels = c("Poor", "Good"))))
})


test_that("are.paired works with auc and mixed roc", {
	expect_true(are.paired(auc(aSAH$outcome, aSAH$wfns), auc(aSAH$outcome, aSAH$ndka)))
	expect_true(are.paired(roc(aSAH$outcome, aSAH$wfns), auc(aSAH$outcome, aSAH$ndka)))
	expect_true(are.paired(auc(aSAH$outcome, aSAH$wfns), roc(aSAH$outcome, aSAH$ndka)))
})

test_that("are.paired return.paired.rocs works", {
	pair <- are.paired(r.wfns, r.ndka, return.paired.rocs = TRUE)
	expect_true(pair)
	expect_identical(attr(pair, "roc1"), r.wfns)
	expect_identical(attr(pair, "roc2"), r.ndka)
})

test_that("are.paired return.paired.rocs works with missing values", {
	aSAH.missing <- aSAH
	aSAH.missing$ndka[1:20] <- NA
	r1 <- roc(aSAH.missing$outcome, aSAH.missing$ndka)
	pair <- are.paired(r1, r.wfns, return.paired.rocs = TRUE)
	expect_true(pair)
	expect_identical(attr(pair, "roc1")$thresholds, roc(aSAH$outcome[21:113], aSAH$ndka[21:113])$thresholds)
	expect_identical(attr(pair, "roc2")$thresholds, roc(aSAH$outcome[21:113], aSAH$wfns[21:113])$thresholds)
})

test_that("are.paired return.paired.rocs doesn't return when unpaired", {
	pair <- are.paired(roc(aSAH$outcome[21:113], aSAH$wfns[21:113]), r.ndka, return.paired.rocs = TRUE)
	expect_null(attributes(pair))
})

test_that("are.paired works with smooth.roc curves", {
	expect_true(are.paired(smooth(r.wfns), smooth(r.ndka)))
	
	# Missing values shouldn't screw up
	aSAH.missing <- aSAH
	aSAH.missing$wfns[1:20] <- NA
	expect_true(are.paired(smooth(roc(aSAH.missing$outcome, aSAH.missing$wfns)), smooth(roc(aSAH.missing$outcome, aSAH.missing$ndka))))
	# Also with different data.frames
	expect_true(are.paired(smooth(roc(aSAH.missing$outcome, aSAH.missing$wfns)), smooth(r.ndka)))
	
	# The following should fail though
	expect_false(are.paired(smooth(roc(aSAH$outcome[21:113], aSAH$wfns[21:113])), smooth(roc(aSAH$outcome, aSAH$ndka))))
	
	# Opposite levels should probably fail
	expect_false(are.paired(smooth(roc(aSAH$outcome, aSAH$wfns, levels = c("Good", "Poor"))), smooth(roc(aSAH$outcome, aSAH$ndka, levels = c("Poor", "Good")))))
})

test_that("are.paired works with auc and mixed roc and smooth", {
	expect_true(are.paired(auc(aSAH$outcome, aSAH$wfns), smooth(roc(aSAH$outcome, aSAH$ndka))))
	expect_true(are.paired(smooth(roc(aSAH$outcome, aSAH$wfns)), auc(aSAH$outcome, aSAH$ndka)))
	expect_true(are.paired(roc(aSAH$outcome, aSAH$wfns), smooth(roc(aSAH$outcome, aSAH$ndka))))
	expect_true(are.paired(smooth(roc(aSAH$outcome, aSAH$wfns)), roc(aSAH$outcome, aSAH$ndka)))
})

test_that("are.paired return.paired.rocs returns smooth curves", {
	aSAH.missing <- aSAH
	aSAH.missing$ndka[1:20] <- NA
	r1 <- roc(aSAH.missing$outcome, aSAH.missing$ndka, smooth=TRUE)
	pair <- are.paired(r1, smooth(r.wfns), return.paired.rocs = TRUE)
	expect_true(pair)
	expect_is(attr(pair, "roc1"), "smooth.roc")
	expect_is(attr(pair, "roc2"), "smooth.roc")
})

test_that("are.paired return.paired.rocs smoothes curves with the right method", {
	skip_slow()
	aSAH.missing <- aSAH
	aSAH.missing$ndka[1:20] <- NA
	smooth.methods <- c("binormal", "density", "fitdistr", "logcondens", "logcondens.smooth")

	for (smooth.method in smooth.methods) {
		r1 <- smooth(roc(aSAH.missing$outcome, aSAH.missing$ndka), method=smooth.method)
		pair <- are.paired(r1, smooth(r.s100b, method=smooth.method), return.paired.rocs = TRUE)
		expect_true(pair)
		expect_identical(attr(pair, "roc1")$smoothing.args$method, smooth.method)
		expect_identical(attr(pair, "roc2")$smoothing.args$method, smooth.method)
	}
})

test_that("are.paired return.paired.rocs doesn't return when unpaired and smooth", {
	pair <- are.paired(smooth(roc(aSAH$outcome[21:113], aSAH$wfns[21:113])), r.ndka, return.paired.rocs = TRUE)
	expect_null(attributes(pair))
	pair <- are.paired(roc(aSAH$outcome[21:113], aSAH$wfns[21:113]), smooth(r.ndka), return.paired.rocs = TRUE)
	expect_null(attributes(pair))
	pair <- are.paired(smooth(roc(aSAH$outcome[21:113], aSAH$wfns[21:113])), smooth(r.ndka), return.paired.rocs = TRUE)
	expect_null(attributes(pair))
})
