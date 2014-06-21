context("Predictor")

test_that("Predictor (base) works", {

	# With no tie, regular
	notie.controls <- c(-1, -0.5, 0, 0.5, 1)
	notie.cases <- notie.controls + 0.25
	
	# Basic operation: return correct cases / controls
	getcontrols <- pROC:::runit_Predictor_getControls(notie.controls, notie.cases)
	expect_that(getcontrols, is_identical_to(notie.controls))
	
	getcases <- pROC:::runit_Predictor_getCases(notie.controls, notie.cases)
	expect_that(getcases, is_identical_to(notie.cases))
	
	# isCases / isControl
	expect_that(pROC:::runit_Predictor_isControl(notie.controls, notie.cases, 4), is_true())
	expect_that(pROC:::runit_Predictor_isControl(notie.controls, notie.cases, 5), is_false())
	expect_that(pROC:::runit_Predictor_isCase(notie.controls, notie.cases, 4), is_false())
	expect_that(pROC:::runit_Predictor_isCase(notie.controls, notie.cases, 5), is_true())
	
	
	# Bracket operator works
	notie <- pROC:::runit_Predictor_bracketOperatorVector(notie.controls, notie.cases)
	expect_that(notie, is_identical_to(c(notie.controls, notie.cases)))
	
	# getOrder works with no ties
	getorder <- pROC:::runit_Predictor_getOrder(notie.controls, notie.cases)
	expect_that(getorder, is_identical_to(as.integer(c(0, 5, 1, 6, 2, 7, 3, 8, 4, 9))))
	
	# Ties should be broken as controls first
	tie.controls <-  seq(-1, 1, 1/4)
	tie.cases <- tie.controls + 0.25
	getorder <- pROC:::runit_Predictor_getOrder(tie.controls, tie.cases)
	expect_that(getorder, is_identical_to(as.integer(c(0, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15, 8, 16, 17))))
	
	
})