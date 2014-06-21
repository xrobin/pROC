context("ci.auc")

test_that("ci.auc works", {
	
	# DeLong
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns))), is_equivalent_to(c(0.748534887819453, 0.823678861788618, 0.898822835757783)))
	
	# Bootstrap
	set.seed(42)
	wfns.full.bootstrap.ci <- c(0.741361788617886, 0.825372628726287, 0.893805047425474)
	wfns.partial.bootstrap.ci <- c(0.0166804878048781, 0.0340497967479675, 0.0549866945197229)
	wfns.corrected.bootstrap.ci <- c(0.561476251604621, 0.652893667094566, 0.763087865893279)
	set.seed(42)
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns), method="bootstrap")), is_equivalent_to(wfns.full.bootstrap.ci))
	set.seed(42);
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns, percent=TRUE), method="bootstrap")), is_equivalent_to(wfns.full.bootstrap.ci * 100))
	# Multiple invocations should yield of the same results
	set.seed(42)
	expect_that(as.numeric(roc(aSAH$outcome, aSAH$wfns, ci.method="bootstrap", ci=TRUE)$ci), is_equivalent_to(wfns.full.bootstrap.ci))
	set.seed(42)
	expect_that(as.numeric(roc(aSAH$outcome, aSAH$wfns, ci.method="bootstrap", ci=TRUE, percent=TRUE)$ci), is_equivalent_to(wfns.full.bootstrap.ci * 100))
	
	# partial.auc
	set.seed(42); 
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns, partial.auc=c(1, .9)), method="bootstrap")), is_equivalent_to(wfns.partial.bootstrap.ci))
	set.seed(42); # method = "bootstrap" is implicit
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns, partial.auc=c(1, .9)))), is_equivalent_to(wfns.partial.bootstrap.ci))
	set.seed(42); 
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns), method="bootstrap", partial.auc=c(1, .9), reuse.auc=FALSE)), is_equivalent_to(wfns.partial.bootstrap.ci))
	set.seed(42); 
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns, partial.auc=c(100, 90), percent=TRUE), method="bootstrap")), is_equivalent_to(wfns.partial.bootstrap.ci * 100))
	set.seed(42); 
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns, percent=TRUE), method="bootstrap", partial.auc=c(100, 90), reuse.auc=FALSE)), is_equivalent_to(wfns.partial.bootstrap.ci * 100))

	# .artial.auc.correct
	set.seed(42); 
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns, partial.auc=c(1, .9), partial.auc.correct = TRUE), method="bootstrap")), is_equivalent_to(wfns.corrected.bootstrap.ci))
	set.seed(42); 
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns), method="bootstrap", partial.auc=c(1, .9), partial.auc.correct = TRUE, reuse.auc=FALSE)), is_equivalent_to(wfns.corrected.bootstrap.ci))
	set.seed(42); 
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns, partial.auc=c(100, 90), partial.auc.correct = TRUE, percent=TRUE), method="bootstrap")), is_equivalent_to(wfns.corrected.bootstrap.ci * 100))
	set.seed(42); 
	expect_that(as.numeric(ci(roc(aSAH$outcome, aSAH$wfns, percent=TRUE), method="bootstrap", partial.auc=c(100, 90), partial.auc.correct = TRUE, reuse.auc=FALSE)), is_equivalent_to(wfns.corrected.bootstrap.ci * 100))
	
})	