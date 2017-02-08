library(pROC)
data(aSAH)

r.wfns <- roc(aSAH$outcome, aSAH$wfns)
r.ndka <- roc(aSAH$outcome, aSAH$ndka)

context("auc")

test_that("full auc works", {
	a.wfns <- auc(r.wfns)
	expect_equal(as.numeric(a.wfns), 0.823678861788618)
	
	a.ndka <- auc(r.ndka)
	expect_equal(as.numeric(a.ndka), 0.611957994579946)
})


test_that("partial auc works on arbitrary intervals", {
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(1, .9))), 0.0334417344173442)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.9, 1))), 0.0334417344173442)
	
	# Arbitrary intervals
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.9, .8))), 0.0598373983739837)
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.5, 0))), 0.488134475939354)
	
	# NDKA
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(1, .9))), 0.0107046070460705)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.9, 1))), 0.0107046070460705)
	
	# Arbitrary intervals
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.9, .8))), 0.0277777777777778)
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.5, 0))), 0.416836043360434)
	
	# Full interval == full auc
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(1, 0))), 0.823678861788618)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(1, 0))), 0.611957994579946)
})


test_that("partial auc works with focus on SE", {
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(1, .9), partial.auc.focus = "se")), 0.0400999322493225)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.9, 1), partial.auc.focus = "se")), 0.0400999322493225)
	
	# Arbitrary intervals
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.9, .8), partial.auc.focus = "se")), 0.0609953703703703)
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.5, 0), partial.auc.focus = "se")), 0.483358739837398)
	
	# NDKA
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(1, .9), partial.auc.focus = "se")), 0.0037940379403794)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.9, 1), partial.auc.focus = "se")), 0.0037940379403794)
	
	# Arbitrary intervals
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.9, .8), partial.auc.focus = "se")), 0.0242547425474255)
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.5, 0), partial.auc.focus = "se")), 0.428523035230352)
	
	# Full interval == full auc
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(1, 0), partial.auc.focus = "se")), 0.823678861788618)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(1, 0), partial.auc.focus = "se")), 0.611957994579946)
})


test_that("partial auc works with correction enabled", {
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(1, .9), partial.auc.correct = TRUE)), 0.649693339038653)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.9, 1), partial.auc.correct = TRUE)), 0.649693339038653)
	
	# Arbitrary intervals
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.9, .8), partial.auc.correct = TRUE)), 0.763749402199904)
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.5, 0), partial.auc.correct = TRUE)), 0.952537903757416)
	
	# NDKA
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(1, .9), partial.auc.correct = TRUE)), 0.530024247610897)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.9, 1), partial.auc.correct = TRUE)), 0.530024247610897)
	
	# Arbitrary intervals
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.9, .8), partial.auc.correct = TRUE)), 0.575163398692811)
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.5, 0), partial.auc.correct = TRUE)), 0.667344173441734)
	
	# Full interval == full auc
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(1, 0), partial.auc.correct = TRUE)), 0.823678861788618)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(1, 0), partial.auc.correct = TRUE)), 0.611957994579946)
})


test_that("partial auc works with focus on SE and correction enabled", {
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(1, .9), partial.auc.focus = "se", partial.auc.correct = TRUE)), 0.68473648552275)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.9, 1), partial.auc.focus = "se", partial.auc.correct = TRUE)), 0.68473648552275)
	
	# Arbitrary intervals
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.9, .8), partial.auc.focus = "se", partial.auc.correct = TRUE)), 0.770561002178649)
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(.5, 0), partial.auc.focus = "se", partial.auc.correct = TRUE)), 0.933434959349593)
	
	# NDKA
	expect_warning(auc(r.ndka, partial.auc = c(1, .9), partial.auc.focus = "se", partial.auc.correct = TRUE))
	expect_identical(as.numeric(auc(r.ndka, partial.auc = c(1, .9), partial.auc.focus = "se", partial.auc.correct = TRUE)), NA_real_)
	# direction is unspecified
	expect_warning(auc(r.ndka, partial.auc = c(.9, 1), partial.auc.focus = "se", partial.auc.correct = TRUE))
	expect_identical(as.numeric(auc(r.ndka, partial.auc = c(.9, 1), partial.auc.focus = "se", partial.auc.correct = TRUE)), NA_real_)
	
	# Arbitrary intervals
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.9, .8), partial.auc.focus = "se", partial.auc.correct = TRUE)), 0.554439662043679)
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(.5, 0), partial.auc.focus = "se", partial.auc.correct = TRUE)), 0.714092140921409)
	
	# Full interval == full auc
	expect_equal(as.numeric(auc(r.wfns, partial.auc = c(1, 0), partial.auc.focus = "se", partial.auc.correct = TRUE)), 0.823678861788618)
	# direction is unspecified
	expect_equal(as.numeric(auc(r.ndka, partial.auc = c(1, 0), partial.auc.focus = "se", partial.auc.correct = TRUE)), 0.611957994579946)
})


test_that("auc can create a roc curve", {
	expect_equal(as.numeric(auc(aSAH$outcome, aSAH$wfns)), as.numeric(auc(r.wfns)))
	expect_equal(as.numeric(auc(aSAH$outcome, aSAH$ndka)), as.numeric(auc(r.ndka)))
	# With formula
	expect_equal(as.numeric(auc(outcome ~ wfns, aSAH)), as.numeric(auc(r.wfns)))
	expect_equal(as.numeric(auc(outcome ~ ndka, aSAH)), as.numeric(auc(r.ndka)))
})
