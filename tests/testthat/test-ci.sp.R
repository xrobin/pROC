library(pROC)
data(aSAH)

context("ci.sp")

# Only test whether ci.sp runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.

# Silence progress bars
options(pROCProgress = list(name = "none"))


for (stratified in c(TRUE, FALSE)) {
	for (test.roc in list(r.s100b, smooth(r.s100b))) {
		test_that("ci.sp with default sensitivities", {
			n <- round(runif(1, 3, 9)) # keep boot.n small
			obtained <- ci.sp(test.roc, boot.n = n, 
							  boot.stratified = stratified, conf.level = .91)
			expect_is(obtained, "ci.sp")
			expect_is(obtained, "ci")
			expect_equal(dim(obtained), c(11, 3))
			expect_equal(attr(obtained, "conf.level"), .91)
			expect_equal(attr(obtained, "boot.n"), n)
			expect_equal(colnames(obtained), c("4.5%", "50%", "95.5%"))
			expect_equal(attr(obtained, "boot.stratified"), stratified)
		})
		
		test_that("ci.sp accepts one sensitivity", {
			n <- round(runif(1, 3, 9)) # keep boot.n small
			obtained <- ci.sp(test.roc, sensitivities = 0.9, boot.n = n, 
							  boot.stratified = stratified, conf.level = .91)
			expect_is(obtained, "ci.sp")
			expect_is(obtained, "ci")
			expect_equal(dim(obtained), c(1, 3))
			expect_equal(attr(obtained, "conf.level"), .91)
			expect_equal(attr(obtained, "boot.n"), n)
			expect_equal(colnames(obtained), c("4.5%", "50%", "95.5%"))
			expect_equal(attr(obtained, "boot.stratified"), stratified)
		})
	}
}
