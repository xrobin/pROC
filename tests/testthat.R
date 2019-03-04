library(testthat)
library(pROC)
data(aSAH)

# Set environment variable RUN_SLOW_TESTS to run the slower tests
run_slow_tests <- identical(Sys.getenv("RUN_SLOW_TESTS"), "true")

test_check("pROC")
