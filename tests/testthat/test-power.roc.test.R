library(pROC)
data(aSAH)

context("power.roc.test")

test_that("power.roc.test basic function", {
  res <- power.roc.test(r.s100b)
  expect_equal(as.numeric(res$auc), as.numeric(r.s100b$auc))
  expect_equal(res$ncases, length(r.s100b$cases))
  expect_equal(res$ncontrols, length(r.s100b$controls))
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.9904833, tolerance = 0.000001)
})

test_that("power.roc.test with percent works", {
	res <- power.roc.test(r.s100b.percent)
	expect_equal(as.numeric(res$auc), as.numeric(r.s100b$auc))
	expect_equal(res$ncases, length(r.s100b$cases))
	expect_equal(res$ncontrols, length(r.s100b$controls))
	expect_equal(res$sig.level, 0.05)
	expect_equal(res$power, 0.9904833, tolerance = 0.000001)
})

test_that("power.roc.test with given auc function", {
  res <- power.roc.test(ncases=41, ncontrols=72, auc=0.73, sig.level=0.05)
  expect_equal(as.numeric(res$auc), 0.73)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.9897453, tolerance = 0.000001)
})


test_that("power.roc.test sig.level can be omitted", {
  res <- power.roc.test(ncases=41, ncontrols=72, auc=0.73)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.9897453, tolerance = 0.000001)
})

test_that("power.roc.test can determine ncases & ncontrols", {
  res <- power.roc.test(auc=r.s100b$auc, sig.level=0.05, power=0.95, kappa=1.7)
  expect_equal(as.numeric(res$auc), as.numeric(r.s100b$auc))
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.95)
  expect_equal(res$ncases, 29.29764, tolerance = 0.000001)
  expect_equal(res$ncontrols, 49.806, tolerance = 0.000001)
})
  
test_that("power.roc.test can determine sig.level", {
  res <- power.roc.test(ncases=41, ncontrols=72, auc=0.73, power=0.95, sig.level=NULL)
  expect_equal(as.numeric(res$auc), 0.73)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(res$power, 0.95)
  expect_equal(res$sig.level, 0.009238584, tolerance = 0.000001)
})

test_that("power.roc.test can determine AUC", {
  res <- power.roc.test(ncases=41, ncontrols=72, sig.level=0.05, power=0.95)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(res$power, 0.95)
  expect_equal(res$sig.level, 0.05)
  expect_equal(as.numeric(res$auc), 0.6961054, tolerance = 0.000001)
})

test_that("power.roc.test can take 2 ROC curves with DeLong variance", {
  res <- power.roc.test(r.ndka, r.wfns)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.7131594, tolerance = 0.000001)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})


test_that("power.roc.test can take 2 percent ROC curves with DeLong variance", {
	res <- power.roc.test(r.ndka.percent, r.wfns.percent)
	expect_equal(res$ncases, 41)
	expect_equal(res$ncontrols, 72)
	expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
	expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
	expect_equal(res$power, 0.7131594, tolerance = 0.000001)
	expect_equal(res$sig.level, 0.05)
	expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test can take 2 ROC curves with Obuchowski variance", {
  res <- power.roc.test(r.ndka, r.wfns, method="obuchowski")
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.8061004, tolerance = 0.000001)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test ncases/ncontrols can take 2 ROC curves with DeLong variance", {
  res <- power.roc.test(r.ndka, r.wfns, power=0.9)
  expect_equal(res$ncases, 64.77777, tolerance = 0.000001)
  expect_equal(res$ncontrols, 113.7561, tolerance = 0.000001)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test ncases/ncontrols can take 2 ROC curves with Obuchowski variance", {
  res <- power.roc.test(r.ndka, r.wfns, power=0.9, method="obuchowski")
  expect_equal(res$ncases, 53.23685, tolerance = 0.000001)
  expect_equal(res$ncontrols, 93.48911, tolerance = 0.000001)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test sig.level can take 2 ROC curves with DeLong variance", {
  res <- power.roc.test(r.ndka, r.wfns, power=0.9, sig.level=NULL)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.1836639, tolerance = 0.000001)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test sig.level can take 2 ROC curves with Obuchowski variance", {
  res <- power.roc.test(r.ndka, r.wfns, power=0.9, sig.level=NULL, method="obuchowski")
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.1150686, tolerance = 0.000001)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test works with partial AUC", {
  skip_slow()
  skip("Bootstrap cannot be tested yet")
  r.wfns.partial <<- roc(aSAH$outcome, aSAH$wfns, quiet = TRUE, partial.auc=c(1, 0.9))
  r.ndka.partial <<- roc(aSAH$outcome, aSAH$ndka, quiet = TRUE, partial.auc=c(1, 0.9))
  power.roc.test(r.wfns.partial, r.ndka.partial, power=0.9)
})


test_that("power.roc.test works with partial AUC", {
  r.wfns.partial <<- roc(aSAH$outcome, aSAH$wfns, quiet = TRUE, partial.auc=c(1, 0.9))
  r.ndka.partial <<- roc(aSAH$outcome, aSAH$ndka, quiet = TRUE, partial.auc=c(1, 0.9))
  res <- power.roc.test(r.wfns.partial, r.ndka.partial, power=0.9, method="obuchowski")
  
  expect_equal(res$ncases, 0.5061498, tolerance = 0.000001)
  expect_equal(res$ncontrols, 0.8888484, tolerance = 0.000001)
  expect_equal(as.numeric(res$auc1), as.numeric(r.wfns.partial$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.ndka.partial$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test works with binormal parameters", {
  ob.params <- list(A1=2.6, B1=1, A2=1.9, B2=1, rn=0.6, ra=0.6, FPR11=0,
                    FPR12=0.2, FPR21=0, FPR22=0.2, delta=0.037) 
  
  res1 <- power.roc.test(ob.params, power=0.8, sig.level=0.05)
  expect_equal(res1$ncases, 107.0238, tolerance = 0.000001)
  expect_equal(res1$ncontrols, 107.0238, tolerance = 0.000001)
  expect_equal(res1$power, 0.8)
  expect_equal(res1$sig.level, 0.05)
  
  res2 <- power.roc.test(ob.params, power=0.8, sig.level=NULL, ncases=107)
  expect_equal(res2$ncases, 107)
  expect_equal(res2$ncontrols, 107)
  expect_equal(res2$power, 0.8)
  expect_equal(res2$sig.level, 0.05004012, tolerance = 0.000001)
  
  res3 <- power.roc.test(ob.params, power=NULL, sig.level=0.05, ncases=107)
  expect_equal(res3$ncases, 107)
  expect_equal(res3$ncontrols, 107)
  expect_equal(res3$sig.level, 0.05)
  expect_equal(res3$power, 0.7999286, tolerance = 0.000001)
})

## With only binormal parameters given
# From example 2 of Obuchowski and McClish, 1997.


test_that("power.roc.test returns correct results from litterature", {
	context("Check results in Obuchowski 2004 Table 4")
	# Note: the table reports at least 10 in each cell, and adapts
	# the complement value to match kappa. So in 0.25/0.95 we have 10/40
	# although both values are < 10.
	# Note2: some values don't match exactly, specifically
	# expected.ncases[0.5, 0.6] and expected.ncontrols[4, 0.6]
	# are off by 1 (< 1%).
	kappas <- c(0.25, 0.5, 1, 2, 4)
	thetas <- c(0.6, 0.7, 0.8, 0.9, 0.95)
	expected.ncontrols <- matrix(c(
		84, 	20,	10,	10,	10,
		101,	25,	10,	10,	10,
		135,	33,	14,	10,	10,
		203,	50,	21,	20,	20,
		#339,	84,	40,	40,	40
		340,	84,	40,	40,	40 # Fixed
	), nrow = 5, byrow = TRUE,
	dimnames = list(kappas, thetas))
	expected.ncases <- matrix(c(
		334,	80,	40,	40,	40,
		#201,	49,	20,	20,	20,
		202,	49,	20,	20,	20, # Fixed
		135,	33,	14,	10,	10,
		102,	25,	11,	10,	10,
		85, 	21,	10,	10,	10
	), nrow = 5, byrow = TRUE, 
	dimnames = list(kappas, thetas))
	
	for (kappa in kappas) {
		for (theta in thetas) {
			context(sprintf("kappa: %s, theta: %s", kappa, theta))
			pr <- power.roc.test(auc=theta, sig.level=0.05, power=0.9, kappa=kappa, alternative="one.sided")
			expect_equal(max(10, ifelse(ceiling(pr$ncases) < 10, 10, 0) * kappa, ceiling(pr$ncontrols)), expected.ncontrols[as.character(kappa), as.character(theta)])
			expect_equal(max(10, ifelse(ceiling(pr$ncontrols) < 10, 10, 0) / kappa, ceiling(pr$ncases)), expected.ncases[as.character(kappa), as.character(theta)])
		}
	}
})


test_that("kappa works with a single ROC curve", {
	# kappa from data
	res <- power.roc.test(r.s100b, sig.level = 0.05, power = 0.9)
	expect_equal(res$ncases, 23.5598674)
	expect_equal(res$ncontrols, 41.3734257)
	expect_equal(res$ncases/res$ncontrols, length(r.s100b$cases)/length(r.s100b$controls))
	# set kappa
	res <- power.roc.test(r.s100b, sig.level = 0.05, power = 0.9, kappa = 1)
	expect_equal(res$ncases, 29.5697422)
	expect_equal(res$ncontrols, 29.5697422)
})


test_that("kappa works with two ROC curves", {
	# kappa from data
	res <- power.roc.test(r.s100b, r.ndka, sig.level = 0.05, power = 0.9)
	expect_equal(res$ncases, 213.117677)
	expect_equal(res$ncontrols, 374.255432)
	expect_equal(res$ncases/res$ncontrols, length(r.s100b$cases)/length(r.s100b$controls))
	# set kappa
	res <- power.roc.test(r.s100b, r.ndka, sig.level = 0.05, power = 0.9, kappa = 1)
	expect_equal(res$ncases, 213.117677)
	expect_equal(res$ncases, 213.117677)
	# ...
})
