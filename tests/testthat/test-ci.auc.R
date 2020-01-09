library(pROC)
data(aSAH)

context("ci.auc")

expected.ci.auc <- c(0.501244999271703, 0.611957994579946, 0.722670989888189)

test_that("ci.auc with delong works", {
	test.ci <- ci.auc(r.ndka)
	expect_is(test.ci, "ci.auc")
	expect_equal(as.numeric(test.ci), expected.ci.auc)
})


test_that("ci.auc with delong and percent works", {
	expect_equal(as.numeric(ci.auc(r.ndka.percent)), expected.ci.auc * 100)
})


test_that("ci.auc works with an auc", {
	expect_equal(as.numeric(ci.auc(auc(r.ndka))), expected.ci.auc)
})


test_that("ci.auc works with a formula", {
	expect_equal(as.numeric(ci.auc(outcome ~ ndka, data = aSAH)), expected.ci.auc)
	expect_equal(as.numeric(ci.auc(outcome ~ ndka, data = aSAH, subset = (gender == "Female"))), 
				 c(0.5261398281, 0.6671428571, 0.8081458862))
})


test_that("ci.auc works with a response, predictor", {
	expect_equal(as.numeric(ci.auc(aSAH$outcome, aSAH$ndka)), expected.ci.auc)
})


test_that("ci.auc works with a direction = >", {
	expect_equal(as.numeric(ci.auc(aSAH$outcome, -aSAH$ndka)), expected.ci.auc)
})


test_that("ci.auc works with a direction = > and percent", {
	expect_equal(as.numeric(ci.auc(aSAH$outcome, -aSAH$ndka, percent = TRUE)), expected.ci.auc * 100)
})


test_that("ci.auc.auc works with a partial AUC from a roc with full AUC", {
	ci.s100b <- ci.auc(r.s100b)
	expect_equal(attr(ci.s100b, "method"), "delong")
	pauc.s100b <- auc(r.s100b, partial.auc = c(1, .9), partial.auc.focus = "se", partial.auc.correct = TRUE)
	ci.pauc.s100b <- ci.auc(pauc.s100b, boot.n = 10, progress = "none")
	expect_equal(attr(ci.pauc.s100b, "method"), "bootstrap")
	expect_equal(attr(attr(ci.pauc.s100b, "auc"), "partial.auc"), c(1, .9))
	expect_equal(attr(attr(ci.pauc.s100b, "auc"), "partial.auc.focus"), "sensitivity")
	expect_equal(attr(attr(ci.pauc.s100b, "auc"), "partial.auc.correct"), TRUE)
})
