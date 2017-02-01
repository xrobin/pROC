library(pROC)
data(aSAH)

# Get static, correct output of the function
source(system.file("extdata", "test-deLongPlacementsCpp-expected.R", package="pROC"))

context("DeLong Placements C++ code works")

for (percent in c(FALSE, TRUE)) {
	for (marker in c("ndka", "wfns", "s100b")) {
		desc <- sprintf("delongPlacementsCpp runs with %s (percent = %s)", marker, percent)
		r <- roc(aSAH$outcome, aSAH[[marker]], percent = percent)
		test_that(desc, {
			placements <- pROC:::delongPlacementsCpp(r)
			expect_equal(placements, expected.placements[[marker]][["forward"]])
		})
	}
	
	for (marker in c("ndka", "wfns", "s100b")) {
		desc <- sprintf("delongPlacementsCpp runs with reversed levels and %s (percent = %s)", marker, percent)
		r <- roc(aSAH$outcome, aSAH[[marker]], levels = c("Poor", "Good"), percent = percent)
		test_that(desc, {
			placements <- pROC:::delongPlacementsCpp(r)
			expect_identical(names(placements), c("theta", "X", "Y"))
		})
	}
	
	for (marker in c("ndka", "wfns", "s100b")) {
		desc <- sprintf("delongPlacementsCpp runs with reversed direction and %s (percent = %s)", marker, percent)
		r <- roc(aSAH$outcome, aSAH[[marker]], direction = ">", percent = percent)
		test_that(desc, {
			placements <- pROC:::delongPlacementsCpp(r)
			expect_identical(names(placements), c("theta", "X", "Y"))
		})
	}
	
	for (marker in c("ndka", "wfns", "s100b")) {
		desc <- sprintf("delongPlacementsCpp runs with reversed levels reversed direction and %s (percent = %s)", marker, percent)
		r <- roc(aSAH$outcome, aSAH[[marker]], levels = c("Poor", "Good"), direction = ">", percent = percent)
		test_that(desc, {
			placements <- pROC:::delongPlacementsCpp(r)
			expect_identical(names(placements), c("theta", "X", "Y"))
		})
	}
}