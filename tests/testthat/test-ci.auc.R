library(pROC)
data(aSAH)

context("ci.auc")

expected.ci.auc <- c(0.497330565621015, 0.611957994579946, 0.715404233940336)

test_that("ci.auc with delong works", {
  test.ci <- ci.auc(r.ndka)
  expect_is(test.ci, "ci.auc")
  expect_equal(as.numeric(test.ci), expected.ci.auc)
})

test_that("ci.auc with delong stays within (0, 1) on the logit scale", {
  r <- roc(controls = c(1, 2, 3, 4, 5), cases = c(4.5, 6, 7, 8, 9), quiet = TRUE)
  expect_equal(as.numeric(ci.auc(r)), c(0.572180762635470, 0.96, 0.997683441755321))
})

test_that("ci.auc with delong handles direction and extreme conf.level", {
  r.rev <- roc(aSAH$outcome, aSAH$ndka, direction = ">", quiet = TRUE)
  expect_equal(as.numeric(ci.auc(r.rev)), 1 - rev(expected.ci.auc))
  expect_equal(as.numeric(ci.auc(r.ndka, conf.level = 1)), c(0, as.numeric(auc(r.ndka)), 1))
  r <- roc(controls = 1:5, cases = c(1.5, 2.5), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, conf.level = 0))
  expect_identical(ci, rep(ci[2], 3))
})

test_that("ci.auc with delong of a perfect curve is a null interval", {
  r <- roc(controls = c(1, 2, 3), cases = c(4, 5, 6), quiet = TRUE)
  expect_warning(ci <- ci.auc(r), "ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.", fixed = TRUE)
  expect_identical(as.numeric(ci), c(1, 1, 1))
  r0 <- roc(controls = c(4, 5, 6), cases = c(1, 2, 3), direction = "<", quiet = TRUE)
  expect_identical(as.numeric(ci.auc(r0)), c(0, 0, 0))
  rt <- roc(controls = rep(1, 3), cases = rep(1, 3), quiet = TRUE)
  expect_identical(as.numeric(ci.auc(rt, conf.level = 1)), c(0.5, 0.5, 0.5))
})


test_that("ci.auc is one-sided with sides", {
  # The bound of a one-sided 95% CI is that of the two-sided 90% CI
  two.sided <- as.numeric(ci.auc(r.ndka, conf.level = 0.9))
  lower <- ci.auc(r.ndka, sides = "lower")
  expect_equal(as.numeric(lower), c(two.sided[1:2], 1))
  expect_identical(attr(lower, "sides"), "lower")
  expect_identical(attr(lower, "conf.level"), 0.95)
  expect_equal(as.numeric(ci.auc(r.ndka, sides = "upper")), c(0, two.sided[2:3]))
  expect_identical(attr(ci.auc(r.ndka), "sides"), "two.sided")
  # Partial matching
  expect_identical(ci.auc(r.ndka, side = "l"), lower)
  expect_output(print(lower), "95% one-sided CI: 0.5161-1 (DeLong)", fixed = TRUE)
  # The default label of the plots shows the bound only
  expect_identical(pROC:::ggroc_auc_label(r.ndka$auc, lower), "AUC: 0.612 (>0.516)")
  expect_identical(pROC:::ggroc_auc_label(r.ndka.percent$auc, ci.auc(r.ndka.percent, sides = "upper")), "AUC: 61.2% (<70.0%)")
  expect_identical(pROC:::ggroc_auc_label(r.ndka$auc, lower, "%.2f %.2f %.2f"), "0.61 0.52 1.00")
  # percent, auc and ci entry points, smoothed curves (ml)
  expect_equal(as.numeric(ci.auc(r.ndka.percent, sides = "lower")), c(two.sided[1:2], 1) * 100)
  expect_equal(as.numeric(ci(auc(r.ndka), sides = "upper")), c(0, two.sided[2:3]))
  s <- smooth(r.ndka)
  expect_equal(as.numeric(ci.auc(s, sides = "lower")), c(as.numeric(ci.auc(s, conf.level = 0.9))[1:2], 1))
  # A null variance leaves no interval on the bounded side
  r <- roc(c(0, 0, 1, 1), c(1, 1, 1, 1), levels = c(0, 1), direction = "<")
  expect_identical(as.numeric(ci.auc(r, sides = "lower")), c(0.5, 0.5, 1))
  expect_identical(as.numeric(ci.auc(r, sides = "upper")), c(0, 0.5, 0.5))
  # The formula method and roc(ci = TRUE)
  expect_identical(ci.auc(outcome ~ ndka, data = aSAH, sides = "lower", quiet = TRUE)[1:3], lower[1:3])
  expect_identical(roc(aSAH$outcome, aSAH$ndka, ci = TRUE, sides = "lower", quiet = TRUE)$ci[1:3], lower[1:3])
  # smooth() and are.paired() replay the sides of a CI they re-use
  r <- r.ndka
  r$ci <- lower
  expect_identical(attr(smooth(r, reuse.ci = TRUE)$ci, "sides"), "lower")
  paired <- are.paired(r, r.s100b, return.paired.rocs = TRUE, reuse.ci = TRUE)
  expect_identical(attr(attr(paired, "roc1")$ci, "sides"), "lower")
  expect_error(ci.auc(r.ndka, sides = "one.sided"))
})

test_that("ci.auc with bootstrap is one-sided with sides", {
  boot.ci <- function(roc, ...) {
    set.seed(42) # the same replicates for each side
    ci.auc(roc, method = "bootstrap", boot.n = 20, ...)
  }
  two.sided <- as.numeric(boot.ci(r.ndka, conf.level = 0.9))
  lower <- boot.ci(r.ndka, sides = "lower")
  expect_equal(as.numeric(lower), c(two.sided[1:2], 1))
  expect_named(lower, c("5%", "50%", "100%"))
  # No replicate: no interval
  expect_true(all(is.na(pROC:::roc_utils_ci_quantile(numeric(0), 0.95, "lower", c(0, 1)))))
  expect_equal(as.numeric(boot.ci(r.ndka, sides = "upper")), c(0, two.sided[2:3]))
  # The other bound is the limit of the partial AUC, corrected or not
  two.sided <- as.numeric(boot.ci(r.s100b.partial1, conf.level = 0.9))
  expect_equal(as.numeric(boot.ci(r.s100b.partial1, sides = "lower")), c(two.sided[1:2], 0.09))
  expect_equal(as.numeric(boot.ci(r.s100b.partial1, sides = "upper")), c(0, two.sided[2:3]))
  pauc <- auc(r.s100b, partial.auc = c(1, 0.9), partial.auc.correct = TRUE)
  expect_equal(boot.ci(pauc, sides = "lower")[[3]], 1)
  expect_equal(boot.ci(pauc, sides = "upper")[[1]], (1 - 0.005 / 0.095) / 2)
  # percent, smoothed curves
  expect_equal(boot.ci(r.s100b.percent.partial1, sides = "lower")[[3]], 9)
  s <- smooth(r.ndka, method = "density")
  expect_equal(as.numeric(boot.ci(s, sides = "upper")), c(0, as.numeric(boot.ci(s, conf.level = 0.9))[2:3]))
})

test_that("ci.auc with delong and percent works", {
  expect_equal(as.numeric(ci.auc(r.ndka.percent)), expected.ci.auc * 100)
})


test_that("ci.auc works with an auc", {
  expect_equal(as.numeric(ci.auc(auc(r.ndka))), expected.ci.auc)
})


test_that("ci.auc works with a formula", {
  expect_equal(as.numeric(ci.auc(outcome ~ ndka, data = aSAH)), expected.ci.auc)
  expect_equal(
    as.numeric(ci.auc(outcome ~ ndka, data = aSAH, subset = (gender == "Female"))),
    c(0.515076172852731, 0.667142857142857, 0.790883412202345)
  )
})


test_that("ci.auc works with a response, predictor", {
  expect_equal(as.numeric(ci.auc(aSAH$outcome, aSAH$ndka)), expected.ci.auc)
})


test_that("ci.auc works with a direction = >", {
  expect_equal(as.numeric(ci.auc(aSAH$outcome, -aSAH$ndka)), expected.ci.auc)
})


test_that("ci.auc works with a direction = > and percent", {
  expect_equal(as.numeric(ci.auc(aSAH$outcome, -aSAH$ndka, percent = TRUE)), expected.ci.auc * 100)
})


test_that("ci.auc.auc works with a partial AUC from a roc with full AUC", {
  ci.s100b <- ci.auc(r.s100b)
  expect_equal(attr(ci.s100b, "method"), "delong")
  pauc.s100b <- auc(r.s100b, partial.auc = c(1, .9), partial.auc.focus = "se", partial.auc.correct = TRUE)
  ci.pauc.s100b <- ci.auc(pauc.s100b, boot.n = 10)
  expect_equal(attr(ci.pauc.s100b, "method"), "bootstrap")
  expect_equal(attr(attr(ci.pauc.s100b, "auc"), "partial.auc"), c(1, .9))
  expect_equal(attr(attr(ci.pauc.s100b, "auc"), "partial.auc.focus"), "sensitivity")
  expect_equal(attr(attr(ci.pauc.s100b, "auc"), "partial.auc.correct"), TRUE)
})


# Only test whether ci.auc runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.
for (stratified in c(TRUE, FALSE)) {
  for (test.roc in list(r.s100b, smooth(r.s100b), auc(r.s100b), r.s100b.partial1, r.s100b.partial2$auc)) {
    test_that("ci.auc with bootstrap works", {
      n <- round(runif(1, 3, 9)) # keep boot.n small
      test.ci <- ci.auc(test.roc, method = "bootstrap", boot.n = n, conf.level = .91, boot.stratified = stratified)
      expect_is(test.ci, "ci.auc")
      expect_equal(attr(test.ci, "conf.level"), .91)
      expect_equal(attr(test.ci, "boot.n"), n)
      expect_equal(attr(test.ci, "names"), c("4.5%", "50%", "95.5%"))
      expect_equal(attr(test.ci, "boot.stratified"), stratified)
      expect_equal(attr(test.ci, "method"), "bootstrap")
    })
  }
}

test_that("ci.auc of a smoothed percent curve with a corrected partial AUC is 100x the fraction one", {
  s <- smooth(roc(aSAH$outcome, aSAH$ndka, partial.auc = c(1, 0.8), partial.auc.correct = TRUE, quiet = TRUE))
  s.percent <- smooth(roc(aSAH$outcome, aSAH$ndka, partial.auc = c(100, 80), partial.auc.correct = TRUE, percent = TRUE, quiet = TRUE))
  set.seed(42)
  ci.fraction <- ci.auc(s, boot.n = 10)
  set.seed(42)
  ci.percent <- ci.auc(s.percent, boot.n = 10)
  expect_false(anyNA(ci.percent))
  expect_equal(as.numeric(ci.percent), 100 * as.numeric(ci.fraction))
})

test_that("ci.auc and ci formula/default methods with smooth = TRUE can recompute the AUC", {
  s <- smooth(roc(aSAH$outcome, aSAH$s100b, auc = FALSE, quiet = TRUE))
  set.seed(42)
  expected <- as.numeric(ci.auc(s, boot.n = 5, reuse.auc = FALSE, progress = "none"))
  set.seed(42)
  expect_equal(as.numeric(ci.auc(outcome ~ s100b, data = aSAH, smooth = TRUE, boot.n = 5, reuse.auc = FALSE, progress = "none", quiet = TRUE)), expected)
  set.seed(42)
  expect_equal(as.numeric(ci.auc(aSAH$outcome, aSAH$s100b, smooth = TRUE, boot.n = 5, auc = FALSE, progress = "none", quiet = TRUE)), expected)
  set.seed(42)
  expect_equal(as.numeric(ci(outcome ~ s100b, data = aSAH, smooth = TRUE, boot.n = 5, reuse.auc = FALSE, progress = "none", quiet = TRUE)), expected)
  set.seed(42)
  expect_equal(as.numeric(ci(aSAH$outcome, aSAH$s100b, smooth = TRUE, boot.n = 5, auc = FALSE, progress = "none", quiet = TRUE)), expected)
})

test_that("ci functions reject a missing or vector conf.level", {
  set.seed(42)
  for (conf.level in list(NA, NA_real_, c(0.9, 0.95), "0.95")) {
    expect_error(ci.auc(r.s100b, conf.level = conf.level), "conf.level must be a single number within the interval [0,1].", fixed = TRUE)
    expect_error(ci.auc(smooth(r.s100b), conf.level = conf.level), "conf.level must be a single number within the interval [0,1].", fixed = TRUE)
    expect_error(ci.se(r.s100b, conf.level = conf.level), "'conf.level' must be a single number within the interval [0,1].", fixed = TRUE)
    expect_error(ci.se(smooth(r.s100b), conf.level = conf.level), "'conf.level' must be a single number within the interval [0,1].", fixed = TRUE)
    expect_error(ci.sp(r.s100b, conf.level = conf.level), "'conf.level' must be a single number within the interval [0,1].", fixed = TRUE)
    expect_error(ci.sp(smooth(r.s100b), conf.level = conf.level), "'conf.level' must be a single number within the interval [0,1].", fixed = TRUE)
    expect_error(ci.thresholds(r.s100b, conf.level = conf.level), "'conf.level' must be a single number within the interval [0,1].", fixed = TRUE)
    expect_error(ci.coords(r.s100b, "best", conf.level = conf.level), "'conf.level' must be a single number within the interval [0,1].", fixed = TRUE)
    expect_error(ci.coords(smooth(r.s100b), "best", conf.level = conf.level), "'conf.level' must be a single number within the interval [0,1].", fixed = TRUE)
  }
  # A logical is a number, as in base R (t.test(conf.level = TRUE))
  expect_identical(as.numeric(ci.auc(r.s100b, conf.level = TRUE)), as.numeric(ci.auc(r.s100b, conf.level = 1)))
})

test_that("ci.auc and ci formula/default methods with smooth = TRUE take method", {
  s <- smooth(r.s100b)
  set.seed(42)
  expected <- as.numeric(ci.auc(s, method = "bootstrap", boot.n = 5, progress = "none"))
  set.seed(42)
  expect_equal(as.numeric(ci.auc(aSAH$outcome, aSAH$s100b, smooth = TRUE, method = "bootstrap", boot.n = 5, progress = "none", quiet = TRUE)), expected)
  set.seed(42)
  expect_equal(as.numeric(ci.auc(outcome ~ s100b, data = aSAH, smooth = TRUE, method = "bootstrap", boot.n = 5, progress = "none", quiet = TRUE)), expected)
  set.seed(42)
  expect_equal(as.numeric(ci(aSAH$outcome, aSAH$s100b, smooth = TRUE, method = "bootstrap", boot.n = 5, progress = "none", quiet = TRUE)), expected)
})

test_that("ci.auc of a smoothed curve takes method and bootstraps DeLong with a warning", {
  s <- smooth(r.s100b)
  set.seed(42)
  expected <- as.numeric(ci.auc(s, method = "bootstrap", boot.n = 5, progress = "none"))
  # DeLong is not supported for smoothed curves: same bootstrap CI, with a warning
  msg <- "DeLong's method is not supported for smoothed ROCs. Using bootstrap instead."
  set.seed(42)
  expect_warning(ci <- ci.auc(s, method = "delong", boot.n = 5, progress = "none"), msg, fixed = TRUE)
  expect_equal(as.numeric(ci), expected)
  set.seed(42)
  expect_warning(ci <- ci.auc(aSAH$outcome, aSAH$s100b, smooth = TRUE, method = "delong", boot.n = 5, progress = "none", quiet = TRUE), msg, fixed = TRUE)
  expect_equal(as.numeric(ci), expected)
  # A CI of the unsmoothed curve computed with DeLong is reused with the
  # default method of the smoothed curve, without warning
  r <- r.s100b
  r$ci <- ci.auc(r, method = "delong")
  expect_no_warning(s.ci <- smooth(r, reuse.ci = TRUE)$ci)
  expect_equal(as.numeric(s.ci), as.numeric(ci.auc(s)))
  expect_identical(attr(s.ci, "method"), "ml")
  # An explicit bootstrap is kept
  r$ci <- ci.auc(r, method = "bootstrap", boot.n = 2, progress = "none")
  expect_identical(attr(smooth(r, reuse.ci = TRUE)$ci, "method"), "bootstrap")
})
