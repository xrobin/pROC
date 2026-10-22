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
  set.seed(42)
  # 2 cases (or controls) among 20: about 12% of the resamples have none,
  # so some of 200 resamples always do (all keep both classes with p < 1e-11)
  resp <- c(rep(0, 18), rep(1, 2))
  pred <- c(1:18, 12.5, 19)
  r.few.cases <- roc(resp, pred, quiet = TRUE)
  r.few.controls <- roc(1 - resp, -pred, quiet = TRUE)
  na.warning <- "NA value\\(s\\) produced during bootstrap were ignored"
  for (r in list(r.few.cases, r.few.controls)) {
    cis <- list()
    expect_warning(cis$se <- ci.se(r, specificities = 0.5, boot.n = 200, boot.stratified = FALSE), na.warning)
    expect_warning(cis$sp <- ci.sp(r, sensitivities = 0.5, boot.n = 200, boot.stratified = FALSE), na.warning)
    expect_warning(cis$sp.se <- ci.coords(r, x = 0.5, input = "sensitivity", ret = "specificity", boot.n = 200, boot.stratified = FALSE), na.warning)
    expect_warning(cis$se.sp <- ci.coords(r, x = 0.5, input = "specificity", ret = "sensitivity", boot.n = 200, boot.stratified = FALSE), na.warning)
    expect_warning(cis$thr <- ci.thresholds(r, thresholds = 15, boot.n = 200, boot.stratified = FALSE), na.warning)
    for (ci in cis) {
      values <- as.numeric(unlist(unclass(ci)))
      expect_false(anyNA(values))
      expect_true(all(values >= 0 & values <= 1))
    }
  }
})

test_that("non-stratified ci.auc and var drop the resamples that lose the controls", {
  set.seed(42)
  # 2 controls out of 10: about 11% of the resamples lose both of them, so
  # some of 200 resamples always do (none does with p < 1e-9)
  r <- roc(c(0, 0, rep(1, 8)), c(1, 3, 2, 4:10), quiet = TRUE)
  na.warning <- "NA value\\(s\\) produced during bootstrap were ignored"
  expect_warning(ci <- ci.auc(r, method = "bootstrap", boot.stratified = FALSE, boot.n = 200), na.warning)
  expect_s3_class(ci, "ci.auc")
  expect_warning(v <- var(r, method = "bootstrap", boot.stratified = FALSE, boot.n = 200), na.warning)
  expect_true(is.numeric(v))
})

test_that("non-stratified bootstraps of ordered predictors drop the resamples that lose a class", {
  set.seed(42)
  pred <- factor(c("a", "c", "b", "b", "c", "c", "d", "d", "d", "d"), levels = c("a", "b", "c", "d"), ordered = TRUE)
  r <- roc(c(0, 0, rep(1, 8)), pred, quiet = TRUE)
  # 2 controls out of 10: some of 200 resamples always lose both (none does with p < 1e-9)
  na.warning <- "NA value\\(s\\) produced during bootstrap were ignored"
  expect_warning(ci <- ci.se(r, specificities = 0.5, boot.stratified = FALSE, boot.n = 200), na.warning)
  expect_s3_class(ci, "ci.se")
  expect_warning(ci <- ci.sp(r, sensitivities = 0.5, boot.stratified = FALSE, boot.n = 200), na.warning)
  expect_s3_class(ci, "ci.sp")
  expect_warning(ci <- ci.auc(r, method = "bootstrap", boot.stratified = FALSE, boot.n = 200), na.warning)
  expect_s3_class(ci, "ci.auc")
  expect_warning(ci <- ci.coords(r, "b", input = "threshold", boot.stratified = FALSE, boot.n = 200), na.warning)
  expect_s3_class(ci, "ci.coords")
})

test_that("boot.n must be a non-negative number", {
  set.seed(42)
  msg <- "'boot.n' must be a non-negative number."
  for (boot.n in list(-1, NA, Inf, "a", numeric(0))) {
    expect_error(ci.auc(r.s100b, method = "bootstrap", boot.n = boot.n), msg, fixed = TRUE)
    expect_error(ci.auc(smooth(r.s100b), boot.n = boot.n), msg, fixed = TRUE)
    expect_error(var(r.s100b, method = "bootstrap", boot.n = boot.n), msg, fixed = TRUE)
    expect_error(cov(r.s100b, r.ndka, method = "bootstrap", boot.n = boot.n), msg, fixed = TRUE)
    expect_error(roc.test(r.s100b, r.ndka, method = "bootstrap", boot.n = boot.n), msg, fixed = TRUE)
    expect_error(roc.test(r.s100b, r.ndka, method = "venkatraman", boot.n = boot.n), msg, fixed = TRUE)
    expect_error(ci.se(r.s100b, boot.n = boot.n), msg, fixed = TRUE)
    expect_error(ci.sp(r.s100b, boot.n = boot.n), msg, fixed = TRUE)
    expect_error(ci.thresholds(r.s100b, boot.n = boot.n), msg, fixed = TRUE)
    expect_error(ci.coords(r.s100b, "best", boot.n = boot.n), msg, fixed = TRUE)
  }
  # An integer-valued double or an integer is fine
  expect_equal(attr(ci.auc(r.s100b, method = "bootstrap", boot.n = 2L), "boot.n"), 2L)
})

test_that("boot.n is read as seq_len() does", {
  # A fraction is truncated
  set.seed(42)
  a <- ci.auc(r.s100b, method = "bootstrap", boot.n = 2.5)
  set.seed(42)
  expect_identical(as.numeric(a), as.numeric(ci.auc(r.s100b, method = "bootstrap", boot.n = 2)))
  set.seed(42)
  t <- roc.test(r.s100b, r.ndka, method = "venkatraman", boot.n = 2.5)
  set.seed(42)
  t2 <- roc.test(r.s100b, r.ndka, method = "venkatraman", boot.n = 2)
  expect_identical(t$parameter, t2$parameter)
  expect_identical(t$p.value, t2$p.value)
  # A vector uses its first element, with a warning
  set.seed(42)
  expect_warning(a <- ci.auc(r.s100b, method = "bootstrap", boot.n = c(2, 10)), "first element used of 'boot.n' argument", fixed = TRUE)
  set.seed(42)
  expect_identical(as.numeric(a), as.numeric(ci.auc(r.s100b, method = "bootstrap", boot.n = 2)))
  # boot.n = 0 runs no replicate: a NA CI or p-value, as before
  expect_identical(as.numeric(ci.auc(r.s100b, method = "bootstrap", boot.n = 0)), rep(NA_real_, 3))
  expect_identical(roc.test(r.s100b, r.ndka, method = "venkatraman", boot.n = 0)$p.value, NaN)
})

# Replicates that fail (cluster F decision on PR #148, rows F1-F15): a
# resample that loses a class, cannot be smoothed, or has no "best" point in
# the partial AUC range gives an NA replicate, dropped with one warning.
# set.seed(1): the seed of the table in the decision. With F9's data, most
# other seeds hit an error of coords() on a degenerate binormal fit (pending,
# not a replicate failure).
f.resp <- c(rep(0, 27), rep(1, 3)) # 3 cases in 30
f.ri <- roc(f.resp, c(1:27, 20.5, 25.5, 28), quiet = TRUE)
f.rc <- roc(c(0, 0, rep(1, 8)), c(1, 3, 2, 4:10), quiet = TRUE) # 2 controls in 10

test_that("F1: ci.auc drops the non-stratified resamples that lose a class", {
  set.seed(1)
  expect_warning(ci <- ci.auc(f.ri, method = "bootstrap", boot.n = 200, boot.stratified = FALSE), "NA value")
  expect_false(anyNA(ci))
})

test_that("F2: ci.sp drops the non-stratified resamples that lose a class", {
  set.seed(1)
  expect_warning(ci <- ci.sp(f.ri, sensitivities = 0.5, boot.n = 200, boot.stratified = FALSE), "NA value")
  expect_false(anyNA(ci))
})

test_that("F3: ci.coords drops the non-stratified resamples that lose a class", {
  set.seed(1)
  expect_warning(ci <- ci.coords(f.ri, x = 0.5, input = "sensitivity", ret = "specificity", boot.n = 200, boot.stratified = FALSE), "NA value")
  expect_false(anyNA(ci$specificity))
})

test_that("F4: ci.thresholds drops the non-stratified resamples that lose a class", {
  set.seed(1)
  expect_warning(ci <- ci.thresholds(f.ri, thresholds = 22, boot.n = 200, boot.stratified = FALSE), "NA value")
  expect_false(anyNA(ci$sensitivity))
  expect_false(anyNA(ci$specificity))
})

test_that("F5: roc.test drops the non-stratified resamples that lose a class", {
  set.seed(1)
  r2 <- roc(f.resp, c(1:27, 26, 2, 28), quiet = TRUE)
  expect_warning(t <- roc.test(f.ri, r2, method = "bootstrap", boot.n = 200, boot.stratified = FALSE), "NA value")
  expect_true(is.finite(t$statistic))
})

test_that("F6: ci.auc drops the resamples that lose all controls", {
  set.seed(1)
  expect_warning(ci <- ci.auc(f.rc, method = "bootstrap", boot.stratified = FALSE, boot.n = 50), "NA value")
  expect_false(anyNA(ci))
})

test_that("F7: var drops the resamples that lose a class", {
  set.seed(1)
  expect_warning(v <- var(f.rc, method = "bootstrap", boot.stratified = FALSE, boot.n = 50), "NA value")
  expect_true(is.finite(v))
})

test_that("F8: ci.se of an ordered predictor drops the resamples that lose a class", {
  ro <- roc(c(0, 0, rep(1, 8)), factor(c("a", "c", "b", "b", "c", "c", "d", "d", "d", "d"), ordered = TRUE), quiet = TRUE)
  set.seed(1)
  expect_warning(ci <- ci.se(ro, specificities = 0.5, boot.stratified = FALSE, boot.n = 50), "NA value")
  expect_false(anyNA(ci))
})

test_that("F9: ci.coords of a smoothed curve keeps the replicates aligned when smoothing fails", {
  s9 <- smooth(roc(rep(0:1, each = 6), c(1, 1, 1, 2, 2, 3, 1, 2, 2, 3, 3, 3), quiet = TRUE))
  set.seed(1)
  expect_warning(ci <- ci.coords(s9, x = c(0.5, 0.9), input = "specificity", ret = c("sensitivity", "specificity"), boot.n = 100), "NA value")
  # Each replicate returns the requested specificities: an aligned CI is the input
  expect_equal(unname(ci$specificity), matrix(c(0.5, 0.9), 2, 3))
})

test_that("F10: ci.se of a smoothed curve with no smoothable resample is NA", {
  s10 <- smooth(roc(c(0, 0, 0, 1, 1, 1, 1, 1), c(1.71, 0.10, 0.80, 0.63, 3.56, 2.54, 1.72, 0.58), quiet = TRUE))
  set.seed(1)
  expect_warning(ci <- ci.se(s10, specificities = c(0.5, 0.9), boot.n = 5), "NA value")
  expect_equal(dim(ci), c(2, 3))
  expect_true(all(is.na(ci)))
})

test_that("F11: ci.coords drops the resamples with no best point in the partial AUC range", {
  r11 <- roc(rep(0:1, each = 5), c(1, 2, 3, 4, 6, 5, 7, 8, 9, 10), partial.auc = c(0.9, 0.7), quiet = TRUE)
  set.seed(1)
  expect_warning(ci <- ci.coords(r11, "best", boot.n = 50), "NA value")
  expect_false(anyNA(ci$specificity))
})

test_that("F12: ci.auc of a smoothed percent curve with a corrected partial AUC", {
  s12 <- smooth(roc(aSAH$outcome, aSAH$ndka, partial.auc = c(100, 80), partial.auc.correct = TRUE, percent = TRUE, quiet = TRUE))
  set.seed(1)
  expect_no_warning(ci <- ci.auc(s12, boot.n = 20))
  expect_false(anyNA(ci))
})

test_that("F13: var of a smoothed curve with a corrected partial AUC does not warn per replicate", {
  s13 <- smooth(roc(rep(0:1, each = 8), c(0.1, 0.5, 0.9, 1.3, 1.6, 2.0, 2.4, 2.9, 0.3, 0.8, 1.2, 1.7, 2.2, 2.6, 3.1, 3.4),
    partial.auc = c(1, 0.9), partial.auc.correct = TRUE, quiet = TRUE
  ))
  set.seed(1)
  expect_no_warning(v <- var(s13, boot.n = 100))
  expect_true(is.finite(v))
})

test_that("F14: roc.test of a smoothed curve with an extra smoothing argument writes no file", {
  wd <- tempfile("pROC-F14-")
  dir.create(wd)
  owd <- setwd(wd)
  on.exit({
    setwd(owd)
    unlink(wd, recursive = TRUE)
  })
  set.seed(1)
  expect_warning(
    t <- roc.test(smooth(r.s100b, method = "density", percent = FALSE), r.ndka, method = "bootstrap", boot.n = 20),
    "Different smoothing parameters"
  )
  expect_true(is.finite(t$statistic))
  expect_false(file.exists("pROC_bug.RData"))
})

test_that("F15: the keyword/level warning of coords repeats once per replicate, as in develop", {
  rb <- roc(c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1),
    factor(c("a", "best", "c", "a", "best", "c", "c", "best", "c", "c"), levels = c("a", "best", "c"), ordered = TRUE),
    quiet = TRUE
  )
  # One warning for x = "best" itself, and one per replicate (the only
  # warnings of this call)
  n.warnings <- 0
  set.seed(1)
  withCallingHandlers(
    ci.coords(rb, "best", ret = c("specificity", "sensitivity"), boot.n = 12),
    warning = function(w) {
      n.warnings <<- n.warnings + 1
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(n.warnings, 13)
})
