library(pROC)
data(aSAH)

context("ci.coords")

# Silence progress bars
options(pROCProgress = list(name = "none"))

test_that("ci.coords accepts threshold output with x=best", {
	expect_error(ci.coords(r.wfns, x="best", input="specificity", ret=c("threshold", "specificity", "sensitivity"), boot.n = 1), NA)
})

test_that("ci.coords rejects threshold output except with x=best", {
	expect_error(ci.coords(r.wfns, x=0.9, input="specificity", ret=c("threshold", "specificity", "sensitivity"), boot.n = 1))
})

test_that("ci.coords accepts threshold output with x=best or if input was threshold", {
	expect_s3_class(ci.coords(r.wfns, x=2, input="threshold", ret=c("threshold", "specificity", "sensitivity"), boot.n = 1), "ci.coords")
	expect_s3_class(ci.coords(r.wfns, x="best", ret=c("threshold", "specificity", "sensitivity"), boot.n = 1), "ci.coords")
})

# Only test whether ci.coords runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.
for (stratified in c(TRUE, FALSE)) {
	for (test.roc in list(r.s100b, smooth(r.s100b))) {
		test_that("ci.coords accepts one x and one ret", {
			obtained <- ci.coords(test.roc, x=0.8, input = "sensitivity", ret="sp", 
								  boot.n=3, conf.level = .91, boot.stratified = stratified)
			expect_equal(attr(obtained, "ret"), "specificity")
			expect_equal(names(obtained), attr(obtained, "ret"))
			for (ci.mat in obtained) {
				expect_equal(dim(ci.mat), c(1, 3))
				expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
			}
		})
		
		test_that("ci.coords accepts one x and multiple ret", {
			obtained <- ci.coords(test.roc, x=0.8, input = "sensitivity", ret=c("sp", "ppv", "tp", "1-sensitivity"), 
								  boot.n=3, conf.level = .91, boot.stratified = stratified)
			expect_equal(attr(obtained, "ret"), c("specificity", "ppv", "tp", "1-sensitivity"))
			expect_equal(names(obtained), attr(obtained, "ret"))
			for (ci.mat in obtained) {
				expect_equal(dim(ci.mat), c(1, 3))
				expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
			}
		})
		
		test_that("ci.coords accepts multiple x and one ret", {
			obtained <- ci.coords(test.roc, x=c(0.8, 0.9), input = "sensitivity", ret="sp", 
								  boot.n=3, conf.level = .91, boot.stratified = stratified)
			expect_equal(attr(obtained, "ret"), "specificity")
			expect_equal(names(obtained), attr(obtained, "ret"))
			for (ci.mat in obtained) {
				expect_equal(dim(ci.mat), c(2, 3))
				expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
			}
		})
		
		test_that("ci.coords accepts multiple x and ret", {
			obtained <- ci.coords(test.roc, x=c(0.9, 0.8), input = "sensitivity", ret=c("sp", "ppv", "tp", "1-se"), 
								  boot.n=3, conf.level = .91, boot.stratified = stratified)
			expect_equal(attr(obtained, "ret"), c("specificity", "ppv", "tp", "1-sensitivity"))
			expect_equal(names(obtained), attr(obtained, "ret"))
			for (ci.mat in obtained) {
				expect_equal(dim(ci.mat), c(2, 3))
				expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
			}
		})
	}
}
