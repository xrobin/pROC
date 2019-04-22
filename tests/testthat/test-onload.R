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