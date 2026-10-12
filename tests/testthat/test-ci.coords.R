library(pROC)
data(aSAH)

context("ci.coords")

test_that("ci.coords accepts threshold output with x=best", {
  expect_warning(
    expect_error(ci.coords(r.wfns, x = "best", input = "specificity", ret = c("threshold", "specificity", "sensitivity"), boot.n = 1), NA),
    "not available for ordered"
  )
})

test_that("ci.coords rejects threshold output except with x=best", {
  expect_error(ci.coords(r.s100b, x = 0.9, input = "specificity", ret = c("threshold", "specificity", "sensitivity"), boot.n = 1))
})

test_that("ci.coords accepts threshold output with x=best or if input was threshold", {
  expect_warning(
    expect_s3_class(ci.coords(r.wfns, x = "2", input = "threshold", ret = c("threshold", "specificity", "sensitivity"), boot.n = 1), "ci.coords"),
    "not available for ordered"
  )
  expect_warning(
    expect_s3_class(ci.coords(r.wfns, x = "best", ret = c("threshold", "specificity", "sensitivity"), boot.n = 1), "ci.coords"),
    "not available for ordered"
  )
})

test_that("ci.coords merges identical best points from empty ordered levels", {
  # Level 2 is unused: thresholds 2 and 3 are the same point in every replicate,
  # and clearly the best one, so that "stop" fails deterministically
  predictor <- ordered(c(rep("1", 20), rep("4", 2), rep("3", 20), rep("1", 2)), levels = c("1", "2", "3", "4"))
  response <- rep(c(0, 1), each = 22)
  r <- roc(response, predictor, quiet = TRUE)
  expect_s3_class(ci.coords(r, x = "best", boot.n = 3), "ci.coords")
  expect_warning(
    expect_s3_class(ci.coords(r, x = "best", ret = c("threshold", "specificity", "sensitivity"), boot.n = 3), "ci.coords"),
    "not available for ordered"
  )
  expect_error(ci.coords(r, x = "best", boot.n = 3, best.policy = "stop"), "More than one")
})

test_that("best.policy unique.stop still stops on distinct best points", {
  res <- data.frame(threshold = c(1.5, 2.5), specificity = c(0.5, 0.5), sensitivity = c(0.8, 0.8))
  expect_error(enforce.best.policy(res, "unique.stop"), "More than one distinct")
  res <- data.frame(specificity = c(0.5, 0.6), sensitivity = c(0.8, 0.7))
  expect_error(enforce.best.policy(res, "unique.stop"), "More than one distinct")
  res <- data.frame(specificity = c(0.5, 0.5), sensitivity = c(0.8, 0.8))
  expect_equal(enforce.best.policy(res, "unique.stop"), res[1, , drop = FALSE])
})

# Only test whether ci.coords runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.
valid_coords_input <- coord.is.monotone <- c(
  "threshold", "sensitivity", "specificity", "tn", "tp", "fn", "fp", "tpr",
  "tnr", "fpr", "fnr", "1-specificity", "1-sensitivity", "recall"
)
for (input in valid_coords_input) {
  for (stratified in c(TRUE, FALSE)) {
    for (test.roc in list(r.s100b, smooth(r.s100b))) {
      # A smoothed curve has no thresholds, so "threshold" is not a valid
      # input for it: coords() rejects that combination and ci.coords() now
      # agrees instead of silently working from absent thresholds.
      if (input == "threshold" && methods::is(test.roc, "smooth.roc")) next
      context(sprintf("input: %s, stratified: %s, class: %s", input, stratified, class(test.roc)))
      test_that("ci.coords accepts one x and one ret", {
        skip_slow()
        obtained <- ci.coords(test.roc,
          x = 0.8, input = input, ret = "sp",
          boot.n = 3, conf.level = .91, boot.stratified = stratified
        )
        expect_equal(attr(obtained, "ret"), "specificity")
        expect_equal(names(obtained), attr(obtained, "ret"))
        for (ci.mat in obtained) {
          expect_equal(dim(ci.mat), c(1, 3))
          expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
        }
      })

      test_that("ci.coords accepts one x and multiple ret", {
        skip_slow()
        obtained <- ci.coords(test.roc,
          x = 0.8, input = input, ret = c("sp", "ppv", "tp", "1-sensitivity"),
          boot.n = 3, conf.level = .91, boot.stratified = stratified
        )
        expect_equal(attr(obtained, "ret"), c("specificity", "ppv", "tp", "1-sensitivity"))
        expect_equal(names(obtained), attr(obtained, "ret"))
        for (ci.mat in obtained) {
          expect_equal(dim(ci.mat), c(1, 3))
          expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
        }
      })

      test_that("ci.coords accepts multiple x and one ret", {
        skip_slow()
        obtained <- ci.coords(test.roc,
          x = c(0.8, 0.9), input = input, ret = "sp",
          boot.n = 3, conf.level = .91, boot.stratified = stratified
        )
        expect_equal(attr(obtained, "ret"), "specificity")
        expect_equal(names(obtained), attr(obtained, "ret"))
        for (ci.mat in obtained) {
          expect_equal(dim(ci.mat), c(2, 3))
          expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
        }
      })

      test_that("ci.coords accepts multiple x and ret", {
        skip_slow()
        obtained <- ci.coords(test.roc,
          x = c(0.9, 0.8), input = input, ret = c("sp", "ppv", "tp", "1-se"),
          boot.n = 3, conf.level = .91, boot.stratified = stratified
        )
        expect_equal(attr(obtained, "ret"), c("specificity", "ppv", "tp", "1-sensitivity"))
        expect_equal(names(obtained), attr(obtained, "ret"))
        for (ci.mat in obtained) {
          expect_equal(dim(ci.mat), c(2, 3))
          expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
        }
      })
    }
  }
}


test_that("ci.coords works on a smoothed curve with x = 'best'", {
  # ci.coords() resolved 'input' from a two-element default through
  # match.arg(several.ok = FALSE), so it failed for every x on a smoothed
  # curve unless 'input' was given explicitly:
  #   Error in match.arg(x, valid.args, several.ok = FALSE):
  #     'arg' must be of length 1
  # and once past that, x = "best" hit a second problem: the bootstrap worker
  # called coords.roc() on a smooth.roc, bypassing the smooth method that
  # fills the absent thresholds with NA.
  s <- smooth(r.s100b)

  set.seed(1)
  best <- ci.coords(s, "best", ret = "sensitivity", boot.n = 20)
  expect_s3_class(best, "ci.coords")
  # The bootstrap median should sit near the point estimate.
  point <- as.numeric(coords(s, "best", ret = "sensitivity"))
  expect_equal(as.numeric(best$sensitivity)[2], point, tolerance = 0.05)

  # A numeric x works without naming 'input', which defaults to specificity.
  set.seed(1)
  by.default <- ci.coords(s, 0.5, ret = "sensitivity", boot.n = 20)
  set.seed(1)
  explicit <- ci.coords(s, 0.5, input = "specificity", ret = "sensitivity", boot.n = 20)
  expect_equal(as.numeric(unlist(by.default)), as.numeric(unlist(explicit)))
})


test_that("ci.coords rejects input = 'threshold' on a smoothed curve", {
  # As coords() does: a smoothed curve has no thresholds.
  s <- smooth(r.s100b)
  expect_error(ci.coords(s, 0.5, input = "threshold", ret = "sp", boot.n = 3))
  expect_error(coords(s, 0.5, input = "threshold", ret = "sp"))
})

test_that("ci.coords on a smoothed curve keeps replicates aligned when smoothing fails", {
  set.seed(42)
  # A replicate that fails to smooth has the shape of a successful one
  failing <- function(roc) stop("cannot smooth")
  res <- pROC:::bootstrap.smooth.coords(1, attr(smooth(r.s100b), "roc"), TRUE,
    x = c(0.5, 0.9), input = "specificity", ret = c("sensitivity", "specificity"),
    best.method = "youden", best.weights = c(1, 0.5), resmooth = failing,
    best.policy = "stop"
  )
  expect_equal(dim(res), c(2, 2))
  expect_named(res, c("sensitivity", "specificity"))
  expect_true(all(is.na(res)))

  # Every other replicate fails to smooth
  n.calls <- 0
  flaky.density <- function(x, ...) {
    n.calls <<- n.calls + 1
    if (n.calls %% 2 == 0) {
      stop("flaky density")
    }
    density(x, ...)
  }
  s <- smooth(r.s100b, method = "density", density.controls = flaky.density)
  expect_warning(
    ci <- ci.coords(s,
      x = c(0.5, 0.9), input = "specificity",
      ret = c("sensitivity", "specificity"), boot.n = 10
    ),
    "NA value\\(s\\) produced during bootstrap were ignored"
  )
  # Every successful replicate returns the input specificity
  for (i in 1:3) {
    expect_equal(as.numeric(ci$specificity[, i]), c(0.5, 0.9))
  }
})

test_that("ci.coords 'best' ignores resamples with no point in the partial AUC range", {
  set.seed(42)
  resp <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
  pred <- c(1, 2, 3, 4, 6, 5, 7, 8, 9, 10)
  # With 5 controls the only specificity in the range is 0.8, missing in many resamples
  r <- roc(resp, pred, partial.auc = c(0.9, 0.7), quiet = TRUE)
  expect_warning(
    ci <- ci.coords(r, "best", ret = c("specificity", "sensitivity"), boot.n = 50),
    "NA value\\(s\\) produced during bootstrap were ignored"
  )
  expect_s3_class(ci, "ci.coords")
  expect_equal(as.numeric(ci$specificity), c(0.8, 0.8, 0.8))
})

test_that("ci.coords 'best' on smoothed curves ignores resamples with no point in the partial AUC range", {
  set.seed(42)
  s <- smooth(roc(aSAH$outcome, aSAH$s100b, quiet = TRUE, partial.auc = c(0.9, 0.899)))
  # About 12% of the resamples have no point in the range: some of 150 always do
  expect_warning(
    ci <- ci.coords(s, "best", boot.n = 150, progress = "none"),
    "NA value\\(s\\) produced during bootstrap were ignored"
  )
  expect_s3_class(ci, "ci.coords")
  # Every resample is restricted to the (tiny) range, or NA if it has no point in it
  sp <- as.numeric(ci$specificity)
  expect_true(all(is.na(sp) | (sp >= 0.899 - 1e-6 & sp <= 0.9 + 1e-6)))
})

test_that("ci.coords of a smoothed curve with a corrected partial AUC does not warn for each resample", {
  set.seed(42)
  # smooth() recomputed the corrected pAUC of every resample, warning when below the diagonal
  # Each value twice (once shifted by 0.05), so that nearly all the resamples can be smoothed
  controls <- c(0.1, 0.5, 0.9, 1.3, 1.6, 2.0, 2.4, 2.9)
  cases <- c(0.3, 0.8, 1.2, 1.7, 2.2, 2.6, 3.1, 3.4)
  resp <- rep(0:1, each = 16)
  pred <- c(controls, controls + 0.05, cases, cases + 0.05)
  s <- smooth(roc(resp, pred, partial.auc = c(1, 0.9), partial.auc.correct = TRUE, quiet = TRUE))
  expect_no_warning(ci.coords(s, x = 0.95, boot.n = 50), message = "Partial AUC correction not defined")
})

test_that("ci.coords warns about unknown or misspelt 'ret' elements", {
  set.seed(42)
  expect_warning(
    res <- ci.coords(r.s100b, 0.9, input = "specificity", ret = c("sensitivity", "Youden"), boot.n = 3),
    "Youden"
  )
  expect_equal(names(res), "sensitivity")
  w <- tryCatch(ci.coords(r.s100b, 0.9, input = "specificity", ret = c("sensitivity", "Youden"), boot.n = 3), warning = identity)
  expect_identical(conditionMessage(w), "Invalid coordinate(s) in 'ret' ignored: 'Youden'.")
  expect_identical(w$call[[1]], as.name("ci.coords.roc"))
  # nothing matches: a clear error
  expect_error(ci.coords(r.s100b, 0.9, input = "specificity", ret = "bogus", boot.n = 3), "Invalid coordinate(s) in 'ret': 'bogus'.", fixed = TRUE)
  expect_error(ci.coords(smooth(r.s100b), "best", ret = "threshold", boot.n = 3), "Invalid coordinate(s) in 'ret': 'threshold'.", fixed = TRUE)
  # the default ret does not warn when the threshold is dropped
  expect_no_warning(
    ci.coords(r.s100b, 0.9, input = "specificity", boot.n = 3),
    message = "'ret'"
  )
})

test_that("ci.coords on a smoothed curve has no threshold column, as coords", {
  s <- smooth(r.s100b)
  # ret = "all" expands without threshold
  res <- ci.coords(s, 0.9, input = "specificity", ret = "all", boot.n = 5)
  expect_equal(names(res), names(coords(s, 0.9, input = "specificity", ret = "all")))
  expect_false("threshold" %in% names(res))
  expect_true(all(res$sensitivity >= 0 & res$sensitivity <= 1))
  expect_true(all(res$tp >= 0 & res$tp <= length(r.s100b$cases)))
  # "threshold" is ignored with a warning, the other columns are unchanged
  set.seed(42)
  expect_warning(
    res.thr <- ci.coords(s, 0.9, input = "specificity", ret = c("threshold", "sensitivity", "youden"), boot.n = 5),
    "threshold"
  )
  set.seed(42)
  res.nothr <- ci.coords(s, 0.9, input = "specificity", ret = c("sensitivity", "youden"), boot.n = 5)
  expect_equal(res.thr, res.nothr)
  # only "threshold": an error, as coords
  expect_error(ci.coords(s, 0.9, input = "specificity", ret = "threshold", boot.n = 5))
  expect_error(coords(s, 0.9, input = "specificity", ret = "threshold"))
})

test_that("ci.coords of a smoothed curve with numeric x uses best.weights", {
  s <- smooth(r.s100b)
  set.seed(42)
  res <- ci.coords(s, 0.9, input = "specificity", ret = c("sensitivity", "youden"), best.weights = c(2, 0.2), boot.n = 5)
  # youden = se + r * sp with r = 2 and sp = 0.9 in every replicate
  expect_equal(res$youden, res$sensitivity + 2 * 0.9, ignore_attr = TRUE)
})

test_that("ci.coords with NA in a numeric x gives NA rows", {
  set.seed(42)
  ci <- ci.coords(r.s100b, x = c(0.5, NA), input = "specificity", ret = "sensitivity", boot.n = 10)
  set.seed(42)
  ci1 <- ci.coords(r.s100b, x = 0.5, input = "specificity", ret = "sensitivity", boot.n = 10)
  expect_identical(ci$sensitivity[1, ], ci1$sensitivity[1, ])
  expect_identical(unname(ci$sensitivity[2, ]), rep(NA_real_, 3))
  expect_output(print(ci), "NA")
})

test_that("ci.coords 'best' on smoothed curves searches every resample within the partial AUC range", {
  set.seed(42)
  # The resamples lost the partial AUC range and returned the global best
  s <- smooth(r.s100b.partial)
  s2 <- roc(aSAH$outcome, aSAH$s100b, partial.auc = c(1, 0.9), smooth = TRUE, quiet = TRUE)
  for (curve in list(s, s2)) {
    ci <- ci.coords(curve, "best", boot.n = 20, progress = "none")
    sp <- as.numeric(ci$specificity)
    expect_true(all(sp >= 0.9 - 1e-6 & sp <= 1 + 1e-6))
  }
  # Focus on sensitivity, in percent
  s3 <- smooth(r.s100b.percent.partial2)
  ci <- ci.coords(s3, "best", boot.n = 20, progress = "none")
  se <- as.numeric(ci$sensitivity)
  expect_true(all(se >= 90 - 1e-4 & se <= 99 + 1e-4))
})
