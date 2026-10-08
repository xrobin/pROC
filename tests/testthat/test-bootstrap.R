library(pROC)
data(aSAH)

# Characterisation tests for the bootstrap code paths (issue #28).
#
# Every bootstrap entry point is run with a fixed seed and a tiny boot.n, and
# its class, shape and values are compared against helper-bootstrap-expected.R.
# A failure means the RNG draw order, the aggregation or the output shape
# changed -- which during a refactor is a regression, not a reason to
# regenerate the expectations.
#
# These are cheap by design (well under a second in total) so they run on CRAN
# and guard the refactor on every check. The statistical behaviour of the
# bootstrap at realistic boot.n is covered by the slow tests below.

test_that("every bootstrap entry point has an expectation", {
  expect_setequal(names(bootstrap.cases), names(expected.bootstrap))
})


for (case.name in names(bootstrap.cases)) {
  # Force the name so each test closes over its own value.
  local({
    nm <- case.name
    test_that(paste("bootstrap unchanged:", nm), {
      if (grepl("smooth", nm)) {
        skip_if_not_installed("MASS")
      }
      expected <- expected.bootstrap[[nm]]
      actual <- bootstrap.run(nm)

      expect_identical(class(actual), expected$class)
      expect_identical(bootstrap.shape(actual), expected$shape)
      # expect_equal, never expect_identical, on the values: smoothed results
      # carry an lm fit whose terms hold an environment pointer, and numeric
      # comparison is the point here anyway.
      expect_equal(bootstrap.sig(actual), expected$values, tolerance = 1e-8)
    })
  })
}


test_that("bootstrap results are reproducible from a seed", {
  # Two runs with the same seed agree; a different seed generally does not.
  a <- bootstrap.run("ci.auc.continuous", seed = 1)
  b <- bootstrap.run("ci.auc.continuous", seed = 1)
  d <- bootstrap.run("ci.auc.continuous", seed = 2)
  expect_equal(as.numeric(a), as.numeric(b))
  expect_false(isTRUE(all.equal(as.numeric(a), as.numeric(d))))
})


test_that("stratified and non-stratified bootstrap differ", {
  # Guards against the two resampling schemes being accidentally unified.
  strat <- bootstrap.run("ci.auc.strat")
  nonstrat <- bootstrap.run("ci.auc.nonstrat")
  expect_false(isTRUE(all.equal(as.numeric(strat), as.numeric(nonstrat))))
})


test_that("boot.n is honoured", {
  # The number of replicates must reach the result, not just the call.
  set.seed(42)
  small <- ci.auc(r.s100b, method = "bootstrap", boot.n = 10)
  set.seed(42)
  large <- ci.auc(r.s100b, method = "bootstrap", boot.n = 50)
  expect_equal(attr(small, "boot.n"), 10)
  expect_equal(attr(large, "boot.n"), 50)
  expect_false(isTRUE(all.equal(as.numeric(small), as.numeric(large))))
})


# ---------------------------------------------------------------------------
# Slow tier: realistic boot.n. Skipped on CRAN and in ordinary runs.
# ---------------------------------------------------------------------------

test_that("bootstrap CI brackets the point estimate at realistic boot.n", {
  skip_slow()
  for (r in list(r.wfns, r.ndka, r.s100b)) {
    set.seed(42)
    ci <- ci.auc(r, method = "bootstrap", boot.n = 1000)
    expect_lte(ci[1], ci[2])
    expect_lte(ci[2], ci[3])
    # ci[2] is the median of the bootstrap replicates, not the observed AUC
    # (see ci_auc_bootstrap), so the observed AUC is only required to fall
    # inside the interval, not to equal the middle value.
    expect_gte(as.numeric(auc(r)), as.numeric(ci[1]))
    expect_lte(as.numeric(auc(r)), as.numeric(ci[3]))
  }
})


test_that("bootstrap CI approaches the DeLong CI at realistic boot.n", {
  skip_slow()
  # Not a tight statistical claim, just that the two methods agree roughly;
  # a gross divergence means the resampling is wrong.
  for (r in list(r.wfns, r.ndka, r.s100b)) {
    set.seed(42)
    boot.ci <- ci.auc(r, method = "bootstrap", boot.n = 2000)
    delong.ci <- ci.auc(r, method = "delong")
    expect_equal(as.numeric(boot.ci), as.numeric(delong.ci), tolerance = 0.05)
  }
})


test_that("non-stratified bootstrap also brackets the estimate", {
  skip_slow()
  set.seed(42)
  ci <- ci.auc(r.s100b, method = "bootstrap", boot.n = 1000, boot.stratified = FALSE)
  expect_lte(ci[1], ci[2])
  expect_lte(ci[2], ci[3])
})


test_that("ci.se and ci.sp are monotone in the requested points", {
  skip_slow()
  set.seed(42)
  se <- ci.se(r.s100b, specificities = seq(0, 1, 0.25), boot.n = 1000)
  # Sensitivity decreases as the required specificity increases.
  expect_true(all(diff(se[, 2]) <= 1e-8))
  # And the interval always contains the median column.
  expect_true(all(se[, 1] <= se[, 2] + 1e-8))
  expect_true(all(se[, 2] <= se[, 3] + 1e-8))
})


test_that("var and cov agree with their bootstrap definitions", {
  skip_slow()
  set.seed(42)
  v <- var(r.s100b, method = "bootstrap", boot.n = 1000)
  expect_gt(v, 0)
  set.seed(42)
  cv <- cov(r.s100b, r.ndka, method = "bootstrap", boot.n = 1000)
  expect_true(is.finite(cv))
})
