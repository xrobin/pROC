library(testthat)
library(pROC)
data(aSAH)


# Set environment variable RUN_SLOW_TESTS to run the slower tests
run_slow_tests <- identical(Sys.getenv("RUN_SLOW_TESTS"), "true")


test_that("We can build basic ROC curves", {
	# These curves will be re-used throughout the tests
	r.wfns <<- roc(aSAH$outcome, aSAH$wfns, quiet = TRUE)
	r.ndka <<- roc(aSAH$outcome, aSAH$ndka, quiet = TRUE)
	r.s100b <<- roc(aSAH$outcome, aSAH$s100b, quiet = TRUE)
	r.wfns.percent <<- roc(aSAH$outcome, aSAH$wfns, percent = TRUE, quiet = TRUE)
	r.ndka.percent <<- roc(aSAH$outcome, aSAH$ndka, percent = TRUE, quiet = TRUE)
	r.s100b.percent <<- roc(aSAH$outcome, aSAH$s100b, percent = TRUE, quiet = TRUE)
})

test_check("pROC")
