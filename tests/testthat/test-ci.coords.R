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
