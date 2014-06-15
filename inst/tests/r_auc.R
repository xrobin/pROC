context("auc")

pauc_correct <- function(auc, min, max) {
	(1+(auc-min)/(max-min))/2
}

test_roc_aucs <- function(response, predictor, expected, info, direction="<", percent=1) {
	# Full
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1)$auc), is_equivalent_to(expected$full * percent), info = info)
	# Partial.uncorrected
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .99) * percent)$auc), is_equivalent_to(expected$sp.1.99 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .99) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.1.99 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .99) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.1.99 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .95) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.1.95 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .95) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.1.95 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .90) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.1.90 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .90) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.1.90 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.9, .8) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.9.8 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.9, .8) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.9.8 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.5, 0) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.5.0 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.5, 0) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.5.0 * percent), info = info)
	
	# Partial.corrected
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .99) * percent, partial.auc.correct = TRUE)$auc), is_equivalent_to(pauc_correct(expected$sp.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .99) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .99) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .95) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.1.95, 0.00125, .05) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .95) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.1.95, 0.00125, .05) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .90) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.1.90, 0.005, .1) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(1, .90) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.1.90, 0.005, .1) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.9, .8) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.9.8, 0.015, .1) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.9, .8) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.9.8, 0.015, .1) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.5, 0) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.5.0, 0.375, .5) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.5, 0) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.5.0, 0.375, .5) * percent), info = info)
	
	# Reverse partial.auc
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.99, 1) * percent, partial.auc.correct = TRUE)$auc), is_equivalent_to(pauc_correct(expected$sp.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.99, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.99, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.95, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.1.95, 0.00125, .05) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.95, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.1.95, 0.00125, .05) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.90, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.1.90, 0.005, .1) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.90, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.1.90, 0.005, .1) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.8, .9) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.9.8, 0.015, .1) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.8, .9) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.9.8, 0.015, .1) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(0, .5) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")$auc), is_equivalent_to(pauc_correct(expected$sp.5.0, 0.375, .5) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(0, .5) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")$auc), is_equivalent_to(pauc_correct(expected$se.5.0, 0.375, .5) * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.99, 1) * percent)$auc), is_equivalent_to(expected$sp.1.99 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.99, 1) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.1.99 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.99, 1) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.1.99 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.95, 1) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.1.95 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.95, 1) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.1.95 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.90, 1) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.1.90 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.90, 1) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.1.90 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.8, .9) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.9.8 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(.8, .9) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.9.8 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(0, .5) * percent, partial.auc.focus = "sp")$auc), is_equivalent_to(expected$sp.5.0 * percent), info = info)
	expect_that(as.numeric(roc(response, predictor, direction=direction, percent = percent > 1, partial.auc = c(0, .5) * percent, partial.auc.focus = "se")$auc), is_equivalent_to(expected$se.5.0 * percent), info = info)
}

test_aucs <- function(roc, expected, info, percent = 1) {
	# Full
	expect_that(as.numeric(auc(roc)), is_equivalent_to(expected$full * percent), info = info)
	# Partial.uncorrected
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .99) * percent)), is_equivalent_to(expected$sp.1.99 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .99) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.1.99 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .99) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.1.99 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .95) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.1.95 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .95) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.1.95 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .90) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.1.90 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .90) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.1.90 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.9, .8) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.9.8 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.9, .8) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.9.8 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.5, 0) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.5.0 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.5, 0) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.5.0 * percent), info = info)
	
	# Partial.corrected
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .99) * percent, partial.auc.correct = TRUE)), is_equivalent_to(pauc_correct(expected$sp.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .99) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .99) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .95) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.1.95, 0.00125, .05) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .95) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.1.95, 0.00125, .05) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .90) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.1.90, 0.005, .1) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(1, .90) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.1.90, 0.005, .1) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.9, .8) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.9.8, 0.015, .1) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.9, .8) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.9.8, 0.015, .1) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.5, 0) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.5.0, 0.375, .5) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.5, 0) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.5.0, 0.375, .5) * percent), info = info)

	# Reverse partial.auc
	expect_that(as.numeric(auc(roc, partial.auc = c(.99, 1) * percent, partial.auc.correct = TRUE)), is_equivalent_to(pauc_correct(expected$sp.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.99, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.99, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.1.99, 5e-05, .01) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.95, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.1.95, 0.00125, .05) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.95, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.1.95, 0.00125, .05) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.90, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.1.90, 0.005, .1) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.90, 1) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.1.90, 0.005, .1) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.8, .9) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.9.8, 0.015, .1) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.8, .9) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.9.8, 0.015, .1) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(0, .5) * percent, partial.auc.correct = TRUE, partial.auc.focus = "sp")), is_equivalent_to(pauc_correct(expected$sp.5.0, 0.375, .5) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(0, .5) * percent, partial.auc.correct = TRUE, partial.auc.focus = "se")), is_equivalent_to(pauc_correct(expected$se.5.0, 0.375, .5) * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.99, 1) * percent)), is_equivalent_to(expected$sp.1.99 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.99, 1) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.1.99 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.99, 1) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.1.99 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.95, 1) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.1.95 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.95, 1) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.1.95 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.90, 1) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.1.90 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.90, 1) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.1.90 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.8, .9) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.9.8 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(.8, .9) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.9.8 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(0, .5) * percent, partial.auc.focus = "sp")), is_equivalent_to(expected$sp.5.0 * percent), info = info)
	expect_that(as.numeric(auc(roc, partial.auc = c(0, .5) * percent, partial.auc.focus = "se")), is_equivalent_to(expected$se.5.0 * percent), info = info)
}

test_aucs_aSAH <- function(marker, expected) {
	test_aucs(roc(aSAH$outcome, aSAH[[marker]]), expected, info = marker)
	test_roc_aucs(aSAH$outcome, aSAH[[marker]], expected, info = marker)
	test_aucs(roc(aSAH$outcome, aSAH[[marker]], percent = TRUE), expected, info = marker, percent = 100)
	test_roc_aucs(aSAH$outcome, aSAH[[marker]], expected, info = marker, percent = 100)
}

test_that("roc works", {
	expected.perfect <- list(
		full = 1,
		sp.1.99 = .01,
		sp.1.95 = .05,
		sp.1.90 = .1,
		sp.9.8 = .1,
		sp.5.0 = .5,
		se.1.99 = .01,
		se.1.95 = .05,
		se.1.90 = .1,
		se.9.8 = .1,
		se.5.0 = .5
		)
	test_aucs(roc(c(0, 1), c(0, 1)), expected.perfect, info = "perfect")
	test_roc_aucs(c(0, 1), c(0, 1), expected.perfect, info = "perfect")
	test_aucs(roc(c(0, 1), c(0, 1), percent = TRUE), expected.perfect, info = "perfect", percent = 100)
	test_roc_aucs(c(0, 1), c(0, 1), expected.perfect, info = "perfect", percent = 100)
	
	expected.null <- list(
		full = .5,
		sp.1.99 = 5e-05,
		sp.1.95 = 0.00125,
		sp.1.90 = 0.005,
		sp.9.8 = 0.015,
		sp.5.0 = 0.375,
		se.1.99 = 5e-05,
		se.1.95 = 0.00125,
		se.1.90 = 0.005,
		se.9.8 = 0.015,
		se.5.0 = 0.375
		)
	test_aucs(roc(c(0, 0, 1, 1), c(0, 1, 0, 1)), expected.null, info = "null")
	test_roc_aucs(c(0, 0, 1, 1), c(0, 1, 0, 1), expected.null, info = "null")
	test_aucs(roc(c(0, 0, 1, 1), c(0, 1, 0, 1), percent = TRUE), expected.null, info = "null", percent = 100)
	test_roc_aucs(c(0, 0, 1, 1), c(0, 1, 0, 1), expected.null, info = "null", percent = 100)

	expected.opposite <- list(
		full = 0,
		sp.1.99 = 0,
		sp.1.95 = 0,
		sp.1.90 = 0,
		sp.9.8 = 0,
		sp.5.0 = 0,
		se.1.99 = 0,
		se.1.95 = 0,
		se.1.90 = 0,
		se.9.8 = 0,
		se.5.0 = 0
		)
	test_aucs(roc(c(0, 1), c(0, 1), direction=">"), expected.opposite, info = "opposite")
	test_roc_aucs(c(0, 1), c(0, 1), expected.opposite, info = "opposite", direction=">")
	test_aucs(roc(c(0, 1), c(0, 1), direction=">", percent = TRUE), expected.opposite, info = "opposite", percent = 100)
	test_roc_aucs(c(0, 1), c(0, 1), expected.opposite, info = "opposite", direction=">", percent = 100)
	
	expected.wfns <- list(
		full = 0.823678861788618,
		sp.1.99 = 0.000395121951219513,
		sp.1.95 = 0.00987804878048782,
		sp.1.90 = 0.0334417344173442,
		sp.9.8 = 0.0598373983739837,
		sp.5.0 = 0.488134475939354,
		se.1.99 = 0.000526736111111112,
		se.1.95 = 0.0131612748419151,
		se.1.90 = 0.0400999322493225,
		se.9.8 = 0.0609953703703703,
		se.5.0 = 0.483358739837398
	)
	test_aucs_aSAH("wfns", expected.wfns)
	
	expected.s100b <- list(
		full = 0.731368563685637,
		sp.1.99 = 0.00292682926829269,
		sp.1.95 = 0.0155487804878049,
		sp.1.90 = 0.0327574525745257,
		sp.9.8 = 0.0478319783197832,
		sp.5.0 = 0.448128387533875,
		se.1.99 = 0,
		se.1.95 = 0.00393038617886179,
		se.1.90 = 0.0137635501355013,
		se.9.8 = 0.0350575880758807,
		se.5.0 = 0.478489159891599
	)
	test_aucs_aSAH("s100b", expected.s100b)
	
	expected.ndka <- list(
		full = 0.611957994579946,
		sp.1.99 = 0.00024390243902439,
		sp.1.95 = 0.00284552845528456,
		sp.1.90 = 0.0107046070460705,
		sp.9.8 = 0.0277777777777778,
		sp.5.0 = 0.416836043360434,
		se.1.99 = 0.000138888888888889,
		se.1.95 = 0.00140582655826558,
		se.1.90 = 0.0037940379403794,
		se.9.8 = 0.0242547425474255,
		se.5.0 = 0.428523035230352
	)
	test_aucs_aSAH("ndka", expected.ndka)

})