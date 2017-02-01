library(pROC)
data(aSAH)
r_ndka <- roc(aSAH$outcome, aSAH$ndka)
r_wfns <- roc(aSAH$outcome, aSAH$wfns)
r_s100b <- roc(aSAH$outcome, aSAH$s100b)

context("DeLong Placements C++ code works")

for (percent in c(FALSE, TRUE)) {
	for (marker in c("ndka", "wfns", "s100b")) {
		desc <- sprintf("delongPlacementsCpp runs with %s (percent = %s)", marker, percent)
		r <- roc(aSAH$outcome, aSAH[[marker]], percent = percent)
		test_that(desc, {
			placements <- pROC:::delongPlacementsCpp(r)
			expect_identical(names(placements), c("theta", "X", "Y"))
			expect_equal(length(placements$theta), 1)
			expect_equal(length(placements$X), sum(aSAH$outcome == "Poor"))
			expect_equal(length(placements$Y), sum(aSAH$outcome == "Good"))
		})
	}
	
	for (marker in c("ndka", "wfns", "s100b")) {
		desc <- sprintf("delongPlacementsCpp runs with reversed levels and %s (percent = %s)", marker, percent)
		r <- roc(aSAH$outcome, aSAH[[marker]], levels = c("Poor", "Good"), percent = percent)
		test_that(desc, {
			placements <- pROC:::delongPlacementsCpp(r)
			expect_identical(names(placements), c("theta", "X", "Y"))
			expect_equal(length(placements$theta), 1)
			expect_equal(length(placements$X), sum(aSAH$outcome == "Good"))
			expect_equal(length(placements$Y), sum(aSAH$outcome == "Poor"))
		})
	}
	
	for (marker in c("ndka", "wfns", "s100b")) {
		desc <- sprintf("delongPlacementsCpp runs with reversed direction and %s (percent = %s)", marker, percent)
		r <- roc(aSAH$outcome, aSAH[[marker]], direction = ">", percent = percent)
		test_that(desc, {
			placements <- pROC:::delongPlacementsCpp(r)
			expect_identical(names(placements), c("theta", "X", "Y"))
			expect_equal(length(placements$theta), 1)
			expect_equal(length(placements$X), sum(aSAH$outcome == "Poor"))
			expect_equal(length(placements$Y), sum(aSAH$outcome == "Good"))
		})
	}
	
	for (marker in c("ndka", "wfns", "s100b")) {
		desc <- sprintf("delongPlacementsCpp runs with reversed levels reversed direction and %s (percent = %s)", marker, percent)
		r <- roc(aSAH$outcome, aSAH[[marker]], levels = c("Poor", "Good"), direction = ">", percent = percent)
		test_that(desc, {
			placements <- pROC:::delongPlacementsCpp(r)
			expect_identical(names(placements), c("theta", "X", "Y"))
			expect_equal(length(placements$theta), 1)
			expect_equal(length(placements$X), sum(aSAH$outcome == "Poor"))
			expect_equal(length(placements$Y), sum(aSAH$outcome == "Good"))
		})
	}
}