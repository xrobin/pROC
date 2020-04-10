library(pROC)

context("onLoad")

test_that("Progress bar is set after pROC is loaded", {
	progress.opt <- getOption("pROCProgress")
	expect_is(progress.opt, "list")
	expect_true("name" %in% names(progress.opt))
})

test_that("Progress bar is set after pROC is loaded", {
	skip("Breaks everything for some reason...")
	options("pROCProgress"=NULL)
	expect_null(getOption("pROCProgress"))
	detach("package:pROC", unload = TRUE)
	library(pROC)
	progress.opt <- getOption("pROCProgress")
	expect_is(progress.opt, "list")
	expect_true("name" %in% names(progress.opt))
})

test_that("Progress bar is set by .onLoad", {
	options("pROCProgress"=NULL)
	expect_null(getOption("pROCProgress"))
	pROC:::.onLoad()
	progress.opt <- getOption("pROCProgress")
	expect_is(progress.opt, "list")
	expect_true("name" %in% names(progress.opt))
})

test_that(".onLoad doesn't override user setting", {
	old.progress.opt <- getOption("pROCProgress")
	options("pROCProgress"=list(dummy=TRUE))
	expect_false("name" %in% names(getOption("pROCProgress")))
	pROC:::.onLoad()
	expect_false("name" %in% names(getOption("pROCProgress")))
	# Restore
	options("pROCProgress"=old.progress.opt)
})

test_that(".parseRcppVersion works", {
	expect_equal(pROC:::.parseRcppVersion("65538"), "1.0.2")
	expect_equal(pROC:::.parseRcppVersion("1"), "0.0.1")
})

test_that("We're running the right Rcpp version", {
	skip_if_not(exists("run_slow_tests") && run_slow_tests, message = "Skipping error-prone Rcpp version check")
	skip_if(Rcpp:::getRcppVersion() == '1.0.3', "RCPP_VERSION broken in 1.0.3")
	
	# This check will often fail, RCPP_VERSION is regularly out of sync,
	# for instance Rcpp 1.0.4.6 has RCPP_VERSION 1.0.4. We can't expect
	# it to be silent, however we still want it to execute without error
	# expect_silent(pROC:::.checkRcppVersion())
	pROC:::.checkRcppVersion()
	
	# Replace the actual RcppVersion with a dummy function that returns 1
	# (= 0.0.1) so we actually see a warning
	saved.RcppVersion <- pROC:::RcppVersion
	assignInNamespace("RcppVersion", function() {return("1")}, "pROC")
	expect_warning(pROC:::.checkRcppVersion())
	# Restore
	assignInNamespace("RcppVersion", saved.RcppVersion, "pROC")
	
})