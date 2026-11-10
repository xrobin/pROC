library(pROC)
data(aSAH)

context("smooth binormal proper")

skip_if_not_installed("mvtnorm")

# Degenerate data of issue #152: only one control scores among the cases
mtcars.glm <- glm(am ~ wt + hp, data = mtcars, family = "binomial")
r.mtcars <- roc(mtcars$am, predict(mtcars.glm, type = "response"), quiet = TRUE)

# The fit takes about 0.3 s: fit it once for the whole file
proper.s100b <- smooth(r.s100b, method = "proper")

# Slopes of the curve, from FPF = 0 to 1: non-increasing on a proper curve
proper_slopes <- function(s) {
  fpf <- rev(1 - s$specificities)
  tpf <- rev(s$sensitivities)
  keep <- diff(fpf) > 0
  (diff(tpf) / diff(fpf))[keep]
}

test_that("Eq. 36 gives the AUCs of Pesce & Metz (2007), Table 1", {
  # Pesce & Metz 2007, p. 826, Table 1: AUC, da, c
  table1 <- rbind(
    c(0.6, 0.1, -0.154357), c(0.6, 0.358, -0.046907), c(0.6, 0.358287, 0),
    c(0.85, 1.24365, -0.428571), c(0.85, 1.46571, -0.14), c(0.85, 1.46574, 0),
    c(0.95, 0.0999998, -0.853738), c(0.95, 2.324, -0.292369), c(0.95, 2.32618, 0)
  )
  for (i in seq_len(nrow(table1))) {
    expect_equal(pROC:::proper_binormal_auc(table1[i, 2], table1[i, 3]), table1[i, 1], tolerance = 1e-6)
  }
})

test_that("Eq. 36 has the closed forms at c = 0, at da = 0 and under reflection", {
  # c = 0: the conventional equal-variance binormal AUC (Metz & Pan p. 12)
  for (da in c(0, 0.5, 1.7, 4)) {
    expect_identical(pROC:::proper_binormal_auc(da, 0), pnorm(da / sqrt(2)))
  }
  # da = 0: F(0, 0; rho) = 1/4 + arcsin(rho) / (2 pi) (Owen 1956, p. 1078)
  for (c in c(-0.9, -0.3, -1e-6, 1e-4, 0.5)) {
    rho <- -(1 - c^2) / (1 + c^2)
    expect_equal(pROC:::proper_binormal_auc(0, c), 1 + asin(rho) / pi, tolerance = 1e-8)
  }
  # the sign of c only reflects the curve (Pesce & Metz p. 819)
  expect_identical(pROC:::proper_binormal_auc(1.2, -0.4), pROC:::proper_binormal_auc(1.2, 0.4))
})

test_that("Eq. 36 is the area under the curve of Eqs. 33-34", {
  for (s in list(proper.s100b, smooth(r.mtcars, method = "proper"), smooth(r.wfns, method = "proper"))) {
    x <- s$specificities
    y <- s$sensitivities
    expect_equal(as.numeric(s$auc), sum(diff(x) * (y[-1] + y[-length(y)]) / 2), tolerance = 1e-4)
  }
})

test_that("the proper fit smooths the degenerate data of issue #152", {
  s <- smooth(r.mtcars, method = "proper")
  expect_is(s, "smooth.roc")
  expect_equal(s$proper$categories, data.frame(controls = c(18L, 0L, 1L, 0L), cases = c(0L, 4L, 0L, 9L)))
  expect_true(all(diff(proper_slopes(s)) <= 0))
  expect_true(s$auc > 0.98 && s$auc < 0.99)
  # the least-squares fit refuses them
  expect_error(smooth(r.mtcars), "ROC curve not smoothable (degenerate data).", fixed = TRUE)
  # its limit (a vertical line: c = 1 > 2 da) gives the start da = 0, c > 0
  # (Pesce & Metz p. 818, footnote), not the chance line (0, 0)
  start <- pROC:::proper_binormal_start(c(18L, 0L, 1L, 0L), c(0L, 4L, 0L, 9L))
  expect_identical(start[["da"]], 0)
  expect_true(start[["c"]] > 0.5)
})

test_that("proper fits of the aSAH curves", {
  s <- proper.s100b
  expect_equal(as.numeric(s$auc), 0.73339, tolerance = 1e-5)
  expect_equal(s$proper$logLik, -338.29553, tolerance = 1e-7)
  expect_length(s$sensitivities, 514)
  expect_true(all(diff(proper_slopes(s)) <= 0))
  # an ordered predictor: one category per level
  s <- smooth(r.wfns, method = "proper")
  expect_equal(s$proper$categories, data.frame(controls = c(37L, 20L, 3L, 8L, 4L), cases = c(2L, 12L, 1L, 8L, 18L)))
  expect_equal(as.numeric(s$auc), 0.84525, tolerance = 1e-5)
  expect_equal(s$proper$logLik, -142.48017, tolerance = 1e-7)
  expect_equal(as.numeric(smooth(r.wfns.num, method = "proper")$auc), as.numeric(s$auc))
})

test_that("the proper fit is symmetric under reflection of the data", {
  # Reversing the classes and the predictor reflects the curve across the
  # negative diagonal: c changes sign, the AUC is the same (Pesce & Metz p. 819)
  s <- proper.s100b
  reflected <- smooth(roc(controls = -aSAH$s100b[aSAH$outcome == "Poor"], cases = -aSAH$s100b[aSAH$outcome == "Good"], direction = "<", quiet = TRUE), method = "proper")
  expect_identical(reflected$proper$categories$controls, rev(s$proper$categories$cases))
  expect_identical(reflected$proper$categories$cases, rev(s$proper$categories$controls))
  expect_equal(reflected$proper$logLik, s$proper$logLik, tolerance = 1e-8)
  expect_equal(as.numeric(reflected$auc), as.numeric(s$auc), tolerance = 1e-4)
  expect_equal(reflected$proper$c, -s$proper$c, tolerance = 1e-3)
})

test_that("the initial cut-points and the likelihood are continuous as c tends to 0", {
  # The bound of Eq. 35 tends to infinity: the cut-points must not be
  # searched from it
  s <- proper.s100b
  k <- s$proper$categories$controls
  l <- s$proper$categories$cases
  v0 <- pROC:::proper_binormal_start_v(1, 0, k, l)
  at0 <- pROC:::proper_binormal_loglik(v0, 1, 0, k, l)
  expect_true(is.finite(at0))
  for (c in c(-1e-12, -1e-8, -1e-6, 1e-12, 1e-8, 1e-6)) {
    v <- pROC:::proper_binormal_start_v(1, c, k, l)
    expect_equal(v, v0, tolerance = 1e-4)
    expect_equal(pROC:::proper_binormal_loglik(v0, 1, c, k, l), at0, tolerance = 1e-6)
  }
})

test_that("a start far from the data finds its initial cut-points", {
  # Categories of a sample of the highly skewed AUC = 0.95 curve of Pesce &
  # Metz Table 1: from the least-squares start (c < 0), the bound of Eq. 35
  # is far from the cut-points, and the initial search over the whole range
  # missed them: "values must be length 1, but FUN(X[[1]]) result is length 0"
  k <- c(12, 0, 33, 0, 4, 0, 1, 0)
  l <- c(0, 2, 0, 2, 0, 1, 0, 45)
  s <- smooth(roc(controls = rep(1:8, k), cases = rep(1:8, l), quiet = TRUE), method = "proper")
  expect_equal(as.numeric(s$auc), 0.97169, tolerance = 1e-4)
  expect_true(s$proper$c < -0.9)
  expect_true(all(diff(proper_slopes(s)) <= 0))
})

test_that("data without a least-squares start are not fitted with the chance line", {
  # No operating point inside the square, yet no perfect performance: the
  # start was (0, 0), a stationary point (AUC 0.5, logLik -17.88)
  r <- roc(controls = 1:10, cases = c(0, 11:20), quiet = TRUE)
  s <- smooth(r, method = "proper")
  expect_null(pROC:::proper_binormal_start(s$proper$categories$controls, s$proper$categories$cases))
  expect_equal(s$proper$logLik, -7.525608, tolerance = 1e-6)
  expect_true(s$auc > 0.97)
})

test_that("data with an empirical AUC below 0.5 are not fitted with the chance line", {
  # The isopleth of AUC <= 0.5 is the chance line, a stationary point: AUC
  # 0.5 exactly (logLik -59.5738)
  k <- c(0, 1, 0, 2, 0, 1, 0, 3, 1, 3, 0)
  l <- c(4, 0, 2, 1, 1, 0, 2, 0, 1, 0, 4)
  r <- roc(controls = rep(seq_along(k), k), cases = rep(seq_along(l), l), direction = "<", quiet = TRUE)
  expect_true(auc(r) < 0.5)
  s <- smooth(r, method = "proper")
  expect_equal(s$proper$logLik, -59.2819, tolerance = 1e-6)
  expect_true(s$auc > 0.55)
})

test_that("data with an empirical AUC of exactly 0.5 are not fitted with the chance line", {
  # Cases on both sides of the controls: AUC 0.5 exactly (logLik -20.7944)
  r <- roc(controls = rep(2, 10), cases = rep(c(1, 3), each = 5), direction = "<", quiet = TRUE)
  expect_equal(as.numeric(auc(r)), 0.5)
  s <- smooth(r, method = "proper")
  expect_equal(s$proper$logLik, -18.88542, tolerance = 1e-6)
  expect_true(s$auc > 0.7)
})

test_that("a fit started at d_a = 0 can leave it", {
  # d_a = |x|: the central difference in x is 0 at x = 0
  s <- smooth(r.wfns, method = "proper")
  k <- s$proper$categories$controls
  l <- s$proper$categories$cases
  from0 <- pROC:::proper_binormal_fit_from(c(da = 0, c = s$proper$c), k, l)
  expect_equal(from0$logLik, s$proper$logLik, tolerance = 1e-8)
  expect_true(from0$da > 1)
})

test_that("the fit lies inside the range of Eq. 35, with its logLik", {
  # The maximum is about c = 0: a step of optim moved atanh(c) across 0 by
  # about 1e-16, below optim's test for an unchanged parameter, and optim
  # returned that point out of the range (logLik -Inf) with the logLik of
  # another point (-13.18)
  s <- smooth(roc(controls = c(1, 1, 2, 2, 3, 3), cases = c(1, 2, 3, 1, 2, 3), quiet = TRUE), method = "proper")
  p <- s$proper
  loglik <- pROC:::proper_binormal_loglik(p$cutpoints, p$da, p$c, p$categories$controls, p$categories$cases)
  expect_identical(loglik, p$logLik)
  # The same from starts at about c = 0 (the step above depends on the
  # platform; these ones changed the sign of c on Linux x86_64)
  k <- p$categories$controls
  l <- p$categories$cases
  for (c in c(-1.1e-16, 1e-16, 3e-16)) {
    fit <- pROC:::proper_binormal_fit_from(c(da = 1e-16, c = c), k, l)
    expect_identical(pROC:::proper_binormal_loglik(fit$cutpoints, fit$da, fit$c, k, l), fit$logLik)
  }
})

test_that("the skewed starts find a maximum the least-squares start misses", {
  # A sample of 25 + 25 cases of a Pesce & Metz population: the
  # least-squares c is about 0, and the starts at c = 0 and +-c all reached
  # logLik -148.4139; the isopleth grid (Metz & Pan p. 20) finds -148.4079
  k <- c(0, 5, 0, 1, 0, 7, 0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 0)
  l <- c(2, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 3, 0, 1, 0, 1, 0, 1, 0, 3, 0, 2, 0, 3)
  s <- smooth(roc(controls = rep(seq_along(k), k), cases = rep(seq_along(l), l), quiet = TRUE), method = "proper")
  expect_equal(s$proper$logLik, -148.40788, tolerance = 1e-7)
  expect_equal(as.numeric(s$auc), 0.6254, tolerance = 1e-3)
})

test_that("direction and percent", {
  s <- proper.s100b
  r.reversed <- roc(aSAH$outcome, -aSAH$s100b, quiet = TRUE)
  expect_identical(r.reversed$direction, ">")
  expect_equal(smooth(r.reversed, method = "proper")$proper[c("da", "c", "logLik", "categories")], s$proper[c("da", "c", "logLik", "categories")])
  sp <- smooth(r.s100b.percent, method = "proper")
  expect_equal(as.numeric(sp$auc), as.numeric(s$auc) * 100)
  expect_equal(sp$sensitivities, s$sensitivities * 100)
  expect_equal(sp$specificities, s$specificities * 100)
  expect_equal(as.numeric(auc(sp, partial.auc = c(100, 80))), as.numeric(auc(s, partial.auc = c(1, 0.8))) * 100)
})

test_that("roc(smooth = TRUE) passes the method, and it is stored for replays", {
  s <- roc(aSAH$outcome, aSAH$s100b, quiet = TRUE, smooth = TRUE, smooth.method = "proper")
  expect_identical(s$smoothing.args$method, "proper")
  expect_equal(s$proper, proper.s100b$proper)
  expect_identical(pROC:::roc_utils_resmooth(r.s100b, s$smoothing.args)$proper, s$proper)
})

test_that("ordered predictors can also be smoothed with the proper method", {
  expect_s3_class(smooth(r.wfns, method = "proper"), "smooth.roc")
  expect_error(smooth(r.wfns, method = "density"),
    "ROC curves of ordered predictors can be smoothed only with binormal or proper smoothing.",
    fixed = TRUE
  )
})

test_that("Obuchowski methods refuse proper curves", {
  # their A and B parameters are those of the conventional binormal model
  s1 <- proper.s100b
  s2 <- smooth(r.ndka, method = "proper")
  msg <- "Obuchowski method is not supported for smoothed ROCs. Use method=\"bootstrap\" instead."
  expect_error(var(s1, method = "obuchowski"), msg, fixed = TRUE)
  expect_error(cov(s1, s2, method = "obuchowski"), msg, fixed = TRUE)
  expect_error(power.roc.test(s1, method = "obuchowski", power = 0.9))
})

test_that("print shows the proper fit", {
  s <- proper.s100b
  # d_a is near 0, where the curve hardly depends on it (Pesce & Metz p. 820)
  expect_output(print(s), "Smoothing: proper \\(d_a: 0\\.0[0-9]+; c: -0\\.38[0-9]*\\)")
})

test_that("perfect performance gives the limit AUC = 1 with a warning", {
  perfect <- "Perfect performance: the proper binormal fit is the limit of AUC = 1 (d_a = 0, c = 1)."
  # Metz & Pan p. 27: perfect separation
  expect_warning(s <- smooth(roc(controls = c(1, 2, 3), cases = c(4, 5, 6), quiet = TRUE), method = "proper"), perfect, fixed = TRUE)
  expect_identical(s$proper[c("da", "c", "logLik")], list(da = 0, c = 1, logLik = NA_real_))
  expect_identical(as.numeric(auc(s)), 1)
  expect_true(all(s$sensitivities[s$specificities < 1] == 1))
  # Points only on the top and left borders, (7/30, 1) and (0, 27/40):
  # consistent with perfect performance within the tied category. The fit
  # stayed at the chance line (logLik -76.4); before, a row of NAs from the
  # operating points failed with "missing value where TRUE/FALSE needed".
  # PROPROC also returns AUC = 1 on such data (return code 3).
  r <- roc(controls = rep(0:1, c(23, 7)), cases = rep(1:4, c(13, 12, 10, 5)), quiet = TRUE, percent = TRUE)
  expect_warning(s <- smooth(r, method = "proper"), perfect, fixed = TRUE)
  expect_identical(as.numeric(auc(s)), 100)
  expect_identical(range(s$sensitivities), c(0, 100))
  expect_error(smooth(r), "ROC curve not smoothable (not enough points).", fixed = TRUE)
})

test_that("a single category is not smoothable", {
  # all values tied
  expect_error(smooth(roc(controls = c(1, 1, 1), cases = c(1, 1), quiet = TRUE), method = "proper"),
    "ROC curve not smoothable (not enough points).",
    fixed = TRUE
  )
})

test_that("a single interior point gives a proper curve through it", {
  r <- roc(controls = c(1, 1, 1, 2, 2), cases = c(1, 2, 2, 2, 2), quiet = TRUE)
  s <- smooth(r, method = "proper")
  expect_true(all(diff(proper_slopes(s)) <= 0))
  # the fit has the likelihood of the saturated model (the data are the
  # categories' observed proportions)
  k <- s$proper$categories$controls
  l <- s$proper$categories$cases
  expect_equal(s$proper$logLik, sum(k * log(k / sum(k))) + sum(l * log(l / sum(l))), tolerance = 1e-6)
})

test_that("bootstrap CI of a proper fit", {
  skip_slow()
  s <- proper.s100b
  ci <- ci.auc(s, boot.n = 10)
  expect_is(ci, "ci.auc")
  expect_true(ci[1] <= ci[2] && ci[2] <= ci[3])
})
