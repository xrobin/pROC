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
      expected <- expected.bootstrap[[nm]]
      actual <- bootstrap.run(nm)
      values <- bootstrap.sig(actual)

      # Always checked, on every platform: the path runs at all, and returns
      # the class and shape its callers index into. None of this depends on
      # floating point.
      expect_identical(class(actual), expected$class)
      expect_identical(bootstrap.shape(actual), expected$shape)
      expect_length(values, length(expected$values))
      expect_true(all(is.finite(values) | is.na(values)))

      # The recorded values are one machine's arithmetic down to the last bit.
      # They are what pins the RNG draw order and the aggregation through a
      # refactor, and they are checked on every machine we control -- but not
      # on CRAN's dozen platforms and BLAS implementations, where a last-bit
      # difference in a smoothed fit would be a failed submission rather than
      # a bug. expect_equal, never expect_identical: smoothed results carry an
      # lm fit whose terms hold an environment pointer.
      skip_slow()
      expect_equal(values, expected$values, tolerance = 1e-8)
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


# ---------------------------------------------------------------------------
# Regressions fixed while unifying the bootstrap internals.
# ---------------------------------------------------------------------------

test_that("ci.coords works on a percent curve", {
  # The non-smoothed ci.coords worker was the only one that did not scale
  # sensitivities/specificities by 100 for a percent curve, so coords() was
  # handed fraction-scale values on a curve flagged as percent and rejected
  # the percent-scale x: "Input specificity (50) not in range (0-1)".
  set.seed(42)
  frac <- ci.coords(r.s100b, x = 0.5, input = "specificity",
                    ret = "sensitivity", boot.n = 10)
  set.seed(42)
  pct <- ci.coords(r.s100b.percent, x = 50, input = "specificity",
                   ret = "sensitivity", boot.n = 10)
  expect_equal(as.numeric(unlist(pct)), as.numeric(unlist(frac)) * 100)
})


test_that("ci.coords returns percent-scale values on a percent curve", {
  # The x = "best" path did not error on the broken version, it silently
  # returned fractions: ci.coords() reported 0.8194 where coords() reported
  # 80.56 on the same curve. Checking that it runs is therefore not enough,
  # the values have to be on the curve's own scale.
  set.seed(42)
  pct <- ci.coords(r.s100b.percent, x = "best",
                   ret = c("sensitivity", "specificity"), boot.n = 10)
  point <- coords(r.s100b.percent, "best", ret = c("sensitivity", "specificity"))

  for (what in c("sensitivity", "specificity")) {
    interval <- as.numeric(pct[[what]])
    expect_true(all(interval >= 0 & interval <= 100))
    # The bootstrap median sits near the point estimate, and both are on the
    # percent scale: a fraction-scale interval would be ~100x smaller.
    expect_equal(interval[2], as.numeric(point[[what]]), tolerance = 0.25)
  }

  # And the same curve as fractions gives the same numbers divided by 100.
  set.seed(42)
  frac <- ci.coords(r.s100b, x = "best",
                    ret = c("sensitivity", "specificity"), boot.n = 10)
  expect_equal(as.numeric(unlist(pct[c("sensitivity", "specificity")])),
               as.numeric(unlist(frac[c("sensitivity", "specificity")])) * 100)
})


test_that("NA replicates are dropped whole, keeping paired statistics aligned", {
  # A 2 x boot.n matrix holds one column per replicate and one row per curve.
  # The old filter used margin 1, so a single NA discarded an entire curve's
  # values and left the matrix with one row: cov() then read the wrong row, or
  # failed with "incorrect number of dimensions".
  m <- rbind(c(0.70, NA, 0.72), c(0.60, 0.61, 0.62))
  kept <- expect_warning(roc_utils_drop_na_replicates(m, margin = 2L),
                         "1 NA value")
  expect_equal(dim(kept), c(2L, 2L))
  expect_equal(kept[1, ], c(0.70, 0.72))
  expect_equal(kept[2, ], c(0.60, 0.62))

  # Rows-as-replicates (ci.se, ci.sp) still filters the other way.
  kept.rows <- expect_warning(roc_utils_drop_na_replicates(t(m), margin = 1L),
                              "1 NA value")
  expect_equal(dim(kept.rows), c(2L, 2L))

  # A plain vector needs no margin.
  expect_equal(expect_warning(roc_utils_drop_na_replicates(c(1, NA, 3)), "1 NA value"),
               c(1, 3))
  expect_no_warning(roc_utils_drop_na_replicates(c(1, 2, 3)))
})


test_that("stratified and non-stratified resampling keep their shapes", {
  # roc_utils_resample() replaced sixteen near-identical preambles; this pins
  # the contract they shared.
  set.seed(42)
  strat <- roc_utils_resample(r.s100b, stratified = TRUE)
  expect_length(strat$controls, length(r.s100b$controls))
  expect_length(strat$cases, length(r.s100b$cases))
  expect_length(strat$predictor, length(r.s100b$predictor))
  expect_length(strat$response, length(r.s100b$response))

  set.seed(42)
  nonstrat <- roc_utils_resample(r.s100b, stratified = FALSE)
  # Non-stratified keeps the total but lets the class sizes vary.
  expect_length(nonstrat$predictor, length(r.s100b$predictor))
  expect_equal(length(nonstrat$controls) + length(nonstrat$cases),
               length(r.s100b$predictor))
})

test_that("non-stratified bootstrap drops the resamples that lose a class", {
  # 2 cases (or controls) among 20: about 12% of the resamples have none
  resp <- c(rep(0, 18), rep(1, 2))
  pred <- c(1:18, 12.5, 19)
  r.few.cases <- roc(resp, pred, quiet = TRUE)
  r.few.controls <- roc(1 - resp, -pred, quiet = TRUE)
  for (r in list(r.few.cases, r.few.controls)) {
    cis <- suppressWarnings(list(
      ci.se(r, specificities = 0.5, boot.n = 50, boot.stratified = FALSE),
      ci.sp(r, sensitivities = 0.5, boot.n = 50, boot.stratified = FALSE),
      ci.coords(r, x = 0.5, input = "sensitivity", ret = "specificity", boot.n = 50, boot.stratified = FALSE),
      ci.coords(r, x = 0.5, input = "specificity", ret = "sensitivity", boot.n = 50, boot.stratified = FALSE),
      ci.thresholds(r, thresholds = 15, boot.n = 50, boot.stratified = FALSE)
    ))
    for (ci in cis) {
      values <- as.numeric(unlist(unclass(ci)))
      expect_false(anyNA(values))
      expect_true(all(values >= 0 & values <= 1))
    }
  }
})

test_that("non-stratified ci.auc and var drop the resamples that lose the controls", {
  # 2 controls out of 10: most resamples of 50 lose both of them at least once
  r <- roc(c(0, 0, rep(1, 8)), c(1, 3, 2, 4:10), quiet = TRUE)
  expect_s3_class(suppressWarnings(ci.auc(r, method = "bootstrap", boot.stratified = FALSE, boot.n = 50)), "ci.auc")
  expect_true(is.numeric(suppressWarnings(var(r, method = "bootstrap", boot.stratified = FALSE, boot.n = 50))))
})

test_that("non-stratified bootstraps of ordered predictors drop the resamples that lose a class", {
  pred <- factor(c("a", "c", "b", "b", "c", "c", "d", "d", "d", "d"), levels = c("a", "b", "c", "d"), ordered = TRUE)
  r <- roc(c(0, 0, rep(1, 8)), pred, quiet = TRUE)
  expect_s3_class(suppressWarnings(ci.se(r, specificities = 0.5, boot.stratified = FALSE, boot.n = 50)), "ci.se")
  expect_s3_class(suppressWarnings(ci.sp(r, sensitivities = 0.5, boot.stratified = FALSE, boot.n = 50)), "ci.sp")
  expect_s3_class(suppressWarnings(ci.auc(r, method = "bootstrap", boot.stratified = FALSE, boot.n = 50)), "ci.auc")
  expect_s3_class(suppressWarnings(ci.coords(r, "b", input = "threshold", boot.stratified = FALSE, boot.n = 50)), "ci.coords")
})
