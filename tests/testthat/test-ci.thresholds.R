library(pROC)
data(aSAH)

context("ci.thresholds")

# Only test whether ci.thresholds runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.

# Silence progress bars
options(pROCProgress = list(name = "none"))


for (stratified in c(TRUE, FALSE)) {
	test_that("ci.threshold accepts thresholds=best", {
		n <- round(runif(1, 3, 9)) # keep boot.n small
		obtained <- ci.thresholds(r.wfns, thresholds="best", boot.n = n, 
								  boot.stratified = stratified, conf.level = .91)
		expect_is(obtained, "ci.thresholds")
		expect_is(obtained, "ci")
		expect_equal(names(obtained), c("specificity", "sensitivity"))
		expect_equal(dim(obtained$specificity), c(1, 3))
		expect_equal(dim(obtained$sensitivity), c(1, 3))
		expect_equal(attr(obtained, "conf.level"), .91)
		expect_equal(attr(obtained, "boot.n"), n)
		expect_equal(colnames(obtained$specificity), c("4.5%", "50%", "95.5%"))
		expect_equal(colnames(obtained$sensitivity), c("4.5%", "50%", "95.5%"))
		expect_equal(attr(obtained, "boot.stratified"), stratified)
	})
	
	test_that("ci.threshold accepts thresholds=best", {
		n <- round(runif(1, 3, 9)) # keep boot.n small
		obtained <- ci.thresholds(r.ndka, thresholds = "local maximas", boot.n = n, 
								  boot.stratified = stratified, conf.level = .91)
		expected.thresholds <- coords(r.ndka, x = "l", ret = "t", transpose = FALSE)$threshold
		expect_is(obtained, "ci.thresholds")
		expect_is(obtained, "ci")
		expect_equal(names(obtained), c("specificity", "sensitivity"))
		expect_equal(dim(obtained$specificity), c(length(expected.thresholds), 3))
		expect_equal(dim(obtained$sensitivity), c(length(expected.thresholds), 3))
		expect_equal(attr(obtained, "conf.level"), .91)
		expect_equal(attr(obtained, "boot.n"), n)
		expect_equal(colnames(obtained$specificity), c("4.5%", "50%", "95.5%"))
		expect_equal(colnames(obtained$sensitivity), c("4.5%", "50%", "95.5%"))
		expect_equal(attr(obtained, "boot.stratified"), stratified)
	})
	
	test_that("ci.threshold accepts numeric thresholds", {
		n <- round(runif(1, 3, 9)) # keep boot.n small
		obtained <- ci.thresholds(r.ndka, thresholds = c(0.5, 0.2), boot.n = n, 
								  boot.stratified = stratified, conf.level = .91)
		expected.thresholds <- coords(r.ndka, x = "l", ret = "t", transpose = FALSE)$threshold
		expect_is(obtained, "ci.thresholds")
		expect_is(obtained, "ci")
		expect_equal(names(obtained), c("specificity", "sensitivity"))
		expect_equal(dim(obtained$specificity), c(2, 3))
		expect_equal(dim(obtained$sensitivity), c(2, 3))
		expect_equal(attr(obtained, "conf.level"), .91)
		expect_equal(attr(obtained, "boot.n"), n)
		expect_equal(colnames(obtained$specificity), c("4.5%", "50%", "95.5%"))
		expect_equal(colnames(obtained$sensitivity), c("4.5%", "50%", "95.5%"))
		expect_equal(attr(obtained, "boot.stratified"), stratified)
	})
}

