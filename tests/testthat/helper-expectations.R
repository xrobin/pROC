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

# Make sure we got a venkatraman test
expect_venkatraman_htest <- function(ht) {
	expect_htest(ht)
	expect_equal(unname(ht$null.value), 0)
	expect_named(ht$null.value, "difference in AUC")
	expect_is(ht$statistic, c("numeric", "integer")) # Can be either?
	expect_named(ht$statistic, "E")
	expect_is(ht$parameter, "numeric")
	expect_named(ht$parameter, "boot.n")
}

# Make sure we got a boostrap test
expect_bootstrap_htest <- function(ht) {
	expect_htest(ht)
	expect_equal(unname(ht$null.value), 0)
	expect_named(ht$null.value, "difference in AUC")
	expect_is(ht$statistic, c("numeric", "integer")) # Can be either?
	expect_named(ht$statistic, "D")
	expect_is(ht$parameter, "numeric")
	expect_named(ht$parameter, c("boot.n", "boot.stratified"))
}
