library(pROC)
data(aSAH)

test_that("var with delong works", {
	expect_equal(cov(r.wfns, r.ndka), -0.000532967856762438)
	expect_equal(cov(r.ndka, r.s100b), -0.000756164938056579)
	expect_equal(cov(r.s100b, r.wfns), 0.00119615567376754)
})


test_that("var works with auc and mixed roc/auc", {
	expect_equal(cov(auc(r.wfns), auc(r.ndka)), -0.000532967856762438)
	expect_equal(cov(auc(r.ndka), r.s100b), -0.000756164938056579)
	expect_equal(cov(r.s100b, auc(r.wfns)), 0.00119615567376754)
})


test_that("var with delong and percent works", {
	expect_equal(cov(r.wfns.percent, r.ndka.percent), -5.32967856762438)
	expect_equal(cov(r.ndka.percent, r.s100b.percent), -7.56164938056579)
	expect_equal(cov(r.s100b.percent, r.wfns.percent), 11.9615567376754)
})


test_that("var with delong, percent and mixed roc/auc works", {
	expect_equal(cov(auc(r.wfns.percent), r.ndka.percent), -5.32967856762438)
	expect_equal(cov(r.ndka.percent, auc(r.s100b.percent)), -7.56164938056579)
	expect_equal(cov(auc(r.s100b.percent), auc(r.wfns.percent)), 11.9615567376754)
})


test_that("var with different auc specifications warns", {
	expect_warning(cov(r.wfns, r.ndka.percent))
	expect_warning(cov(r.wfns.percent, r.ndka))
	# Also mixing auc/roc
	expect_warning(cov(auc(r.wfns), r.ndka.percent))
	expect_warning(cov(r.wfns, auc(r.ndka.percent)))
	expect_warning(cov(r.wfns, auc(r.ndka.percent)))
})
