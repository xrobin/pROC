library(pROC)
data(aSAH)

# Make sure the value looks like a p value.
expect_p_value <- function(p.value) {
	expect_is(p.value, "numeric")
	expect_lte(p.value, 1)
	expect_gte(p.value, 0)
}

# Make sure we got a htest
expect_htest <- function(ht) {
	expect_is(ht, "htest")
	expect_p_value(ht$p.value)
}

expect_venkatraman_htest <- function(ht) {
	expect_htest(ht)
	expect_equal(unname(ht$null.value), 0)
	expect_named(ht$null.value, "difference in AUC")
	expect_is(ht$statistic, c("numeric", "integer")) # Can be either?
	expect_named(ht$statistic, "E")
	expect_is(ht$parameter, "numeric")
	expect_named(ht$parameter, "boot.n")
}

test_that("paired venkatraman works as expected", {
	skip_if_not(exists("run_slow_tests") && run_slow_tests, message = "Slow test skipped")
	ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12)
	expect_venkatraman_htest(ht)
	expect_equal(ht$alternative, "two.sided")
	expect_equal(ht$method, "Venkatraman's test for two paired ROC curves")
	expect_equal(unname(ht$parameter), 12)
	# Test output
	ht$statistic <- c(E = 42)
	ht$p.value <- 0
	expect_known_output(print(ht), "print_output/roc.test-venkatraman.paired")
})

test_that("unpaired venkatraman works as expected", {
	skip_if_not(exists("run_slow_tests") && run_slow_tests, message = "Slow test skipped")
	expect_warning(ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, paired = FALSE), "paired")
	expect_venkatraman_htest(ht)
	expect_equal(ht$alternative, "two.sided")
	expect_equal(ht$method, "Venkatraman's test for two unpaired ROC curves")
	expect_equal(unname(ht$parameter), 12)
	# Test output
	ht$statistic <- c(E = 41)
	ht$p.value <- 0.548347196932
	expect_known_output(print(ht), "print_output/roc.test-venkatraman.unpaired")
})

test_that("non stratified venkatraman works as expected", {
	skip_if_not(exists("run_slow_tests") && run_slow_tests, message = "Slow test skipped")
	ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, boot.stratified = FALSE)
	expect_venkatraman_htest(ht)
	expect_equal(ht$alternative, "two.sided")
	expect_equal(ht$method, "Venkatraman's test for two paired ROC curves")
	expect_equal(unname(ht$parameter), 12)
	# Test output
	ht$statistic <- c(E = 43)
	ht$p.value <- 0.05
	expect_known_output(print(ht), "print_output/roc.test-venkatraman.unstratified")
})

test_that("non stratified, unpaired venkatraman works as expected", {
	skip_if_not(exists("run_slow_tests") && run_slow_tests, message = "Slow test skipped")
	expect_warning(ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, boot.stratified = FALSE, paired = FALSE), "paired")
	expect_venkatraman_htest(ht)
	expect_equal(ht$alternative, "two.sided")
	expect_equal(ht$method, "Venkatraman's test for two unpaired ROC curves")
	expect_equal(unname(ht$parameter), 12)
	# Test output
	ht$statistic <- c(E = 43)
	ht$p.value <- 0.05
	expect_known_output(print(ht), "print_output/roc.test-venkatraman.unpaired.unstratified")
})


