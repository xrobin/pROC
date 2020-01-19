library(pROC)
data(aSAH)

test_that("cov with delong works", {
	expect_equal(cov(r.wfns, r.ndka), -0.000532967856762438)
	expect_equal(cov(r.ndka, r.s100b), -0.000756164938056579)
	expect_equal(cov(r.s100b, r.wfns), 0.00119615567376754)
})


test_that("cov with obuchowski works", {
	expect_equal(cov(r.wfns, r.ndka, method = "obuchowski"), -3.917223e-06)
	expect_equal(cov(r.ndka, r.s100b, method = "obuchowski"), 0.0007945308)
	expect_equal(cov(r.s100b, r.wfns, method = "obuchowski"), 0.0008560803)
})


test_that("cov works with auc and mixed roc/auc", {
	expect_equal(cov(auc(r.wfns), auc(r.ndka)), -0.000532967856762438)
	expect_equal(cov(auc(r.ndka), r.s100b), -0.000756164938056579)
	expect_equal(cov(r.s100b, auc(r.wfns)), 0.00119615567376754)
})


test_that("cov with delong and percent works", {
	expect_equal(cov(r.wfns.percent, r.ndka.percent), -5.32967856762438)
	expect_equal(cov(r.ndka.percent, r.s100b.percent), -7.56164938056579)
	expect_equal(cov(r.s100b.percent, r.wfns.percent), 11.9615567376754)
})


test_that("cov with delong, percent and mixed roc/auc works", {
	expect_equal(cov(auc(r.wfns.percent), r.ndka.percent), -5.32967856762438)
	expect_equal(cov(r.ndka.percent, auc(r.s100b.percent)), -7.56164938056579)
	expect_equal(cov(auc(r.s100b.percent), auc(r.wfns.percent)), 11.9615567376754)
})


test_that("cov with obuchowski, percent and mixed roc/auc works", {
	expect_equal(cov(auc(r.wfns.percent), r.ndka.percent, method = "obuchowski"), -0.03917223)
	expect_equal(cov(r.ndka.percent, auc(r.s100b.percent), method = "obuchowski"), 7.9453082)
	expect_equal(cov(auc(r.s100b.percent), auc(r.wfns.percent), method = "obuchowski"), 8.560803)
})


test_that("cov with different auc specifications warns", {
	expect_warning(cov(r.wfns, r.ndka.percent))
	expect_warning(cov(r.wfns.percent, r.ndka))
	# Also mixing auc/roc
	expect_warning(cov(auc(r.wfns), r.ndka.percent))
	expect_warning(cov(r.wfns, auc(r.ndka.percent)))
	expect_warning(cov(r.wfns, auc(r.ndka.percent)))
})


test_that("cov with delong, percent and direction = >", {
	expect_equal(cov(r.ndka.percent, r.s100b.percent), -7.56164938056579)
})


test_that("cov with delong, percent, direction = > and mixed roc/auc", {
	r1 <- roc(aSAH$outcome, -aSAH$ndka, percent=TRUE)
	r2 <- roc(aSAH$outcome, -aSAH$s100b, percent=TRUE)
	expect_equal(cov(r1, r2), -7.56164938056579)
	expect_equal(cov(auc(r1), auc(r2)), -7.56164938056579)
	expect_equal(cov(auc(r1), r2), -7.56164938056579)
	expect_equal(cov(r1, auc(r2)), -7.56164938056579)
})


test_that("cov with bootstrap works", {
	skip_slow()
	if (paste0(R.version$major, ".", R.version$minor) >= "3.6.0") {
		RNGkind(sample.kind="Rounding")
	}
	set.seed(42)
	expect_equal(cov(r.wfns, r.ndka, method = "bootstrap", boot.n = 100), -0.0004581385)
	expect_equal(cov(r.ndka.percent, r.s100b.percent, method = "bootstrap", boot.n = 100), -6.312029126)
	expect_equal(cov(r.s100b.partial1, r.wfns.partial1, method = "bootstrap", boot.n = 100), 2.899627e-05)
	expect_equal(cov(r.wfns, r.ndka, method = "bootstrap", boot.n = 100, boot.stratified = FALSE), -0.000419791)
})
