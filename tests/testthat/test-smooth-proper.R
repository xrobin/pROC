library(pROC)
data(aSAH)

context("smooth binormal proper")

skip_if_not_installed("mvtnorm")

# Degenerate data of issue #152: only one control scores among the cases
mtcars.glm <- glm(am ~ wt + hp, data = mtcars, family = "binomial")
r.mtcars <- roc(mtcars$am, predict(mtcars.glm, type = "response"), quiet = TRUE)

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
  for (s in list(smooth(r.s100b, fit = "proper"), smooth(r.mtcars, fit = "proper"), smooth(r.wfns, fit = "proper"))) {
    x <- s$specificities
    y <- s$sensitivities
    expect_equal(as.numeric(s$auc), sum(diff(x) * (y[-1] + y[-length(y)]) / 2), tolerance = 1e-4)
  }
})

test_that("the proper fit smooths the degenerate data of issue #152", {
  s <- smooth(r.mtcars, fit = "proper")
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
  s <- smooth(r.s100b, fit = "proper")
  expect_equal(as.numeric(s$auc), 0.73339, tolerance = 1e-5)
  expect_equal(s$proper$logLik, -338.29553, tolerance = 1e-7)
  expect_length(s$sensitivities, 514)
  expect_true(all(diff(proper_slopes(s)) <= 0))
  # an ordered predictor: one category per level
  s <- smooth(r.wfns, fit = "proper")
  expect_equal(s$proper$categories, data.frame(controls = c(37L, 20L, 3L, 8L, 4L), cases = c(2L, 12L, 1L, 8L, 18L)))
  expect_equal(as.numeric(s$auc), 0.84525, tolerance = 1e-5)
  expect_equal(s$proper$logLik, -142.48017, tolerance = 1e-7)
  expect_equal(as.numeric(smooth(r.wfns.num, fit = "proper")$auc), as.numeric(s$auc))
})

test_that("the proper fit is symmetric under reflection of the data", {
  # Reversing the classes and the predictor reflects the curve across the
  # negative diagonal: c changes sign, the AUC is the same (Pesce & Metz p. 819)
  s <- smooth(r.s100b, fit = "proper")
  reflected <- smooth(roc(controls = -aSAH$s100b[aSAH$outcome == "Poor"], cases = -aSAH$s100b[aSAH$outcome == "Good"], direction = "<", quiet = TRUE), fit = "proper")
  expect_identical(reflected$proper$categories$controls, rev(s$proper$categories$cases))
  expect_identical(reflected$proper$categories$cases, rev(s$proper$categories$controls))
  expect_equal(reflected$proper$logLik, s$proper$logLik, tolerance = 1e-8)
  expect_equal(as.numeric(reflected$auc), as.numeric(s$auc), tolerance = 1e-4)
  expect_equal(reflected$proper$c, -s$proper$c, tolerance = 1e-3)
})

test_that("the initial cut-points and the likelihood are continuous as c tends to 0", {
  # The bound of Eq. 35 tends to infinity: the cut-points must not be
  # searched from it
  s <- smooth(r.s100b, fit = "proper")
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
  s <- smooth(roc(controls = rep(1:8, k), cases = rep(1:8, l), quiet = TRUE), fit = "proper")
  expect_equal(as.numeric(s$auc), 0.97169, tolerance = 1e-4)
  expect_true(s$proper$c < -0.9)
  expect_true(all(diff(proper_slopes(s)) <= 0))
})

test_that("data without a least-squares start are not fitted with the chance line", {
  # No operating point inside the square, yet no perfect performance: the
  # start was (0, 0), a stationary point (AUC 0.5, logLik -17.88)
  r <- roc(controls = 1:10, cases = c(0, 11:20), quiet = TRUE)
  s <- smooth(r, fit = "proper")
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
  s <- smooth(r, fit = "proper")
  expect_equal(s$proper$logLik, -59.2819, tolerance = 1e-6)
  expect_true(s$auc > 0.55)
})

test_that("data with an empirical AUC of exactly 0.5 are not fitted with the chance line", {
  # Cases on both sides of the controls: AUC 0.5 exactly (logLik -20.7944)
  r <- roc(controls = rep(2, 10), cases = rep(c(1, 3), each = 5), direction = "<", quiet = TRUE)
  expect_equal(as.numeric(auc(r)), 0.5)
  s <- smooth(r, fit = "proper")
  expect_equal(s$proper$logLik, -18.88542, tolerance = 1e-6)
  expect_true(s$auc > 0.7)
})

test_that("a fit started at d_a = 0 can leave it", {
  # d_a = |x|: the central difference in x is 0 at x = 0
  s <- smooth(r.wfns, fit = "proper")
  k <- s$proper$categories$controls
  l <- s$proper$categories$cases
  from0 <- pROC:::proper_binormal_fit_from(c(da = 0, c = s$proper$c), k, l)
  expect_equal(from0$logLik, s$proper$logLik, tolerance = 1e-8)
  expect_true(from0$da > 1)
})

test_that("the skewed starts find a maximum the least-squares start misses", {
  # A sample of 25 + 25 cases of a Pesce & Metz population: the
  # least-squares c is about 0, and the starts at c = 0 and +-c all reached
  # logLik -148.4139; the isopleth grid (Metz & Pan p. 20) finds -148.4079
  k <- c(0, 5, 0, 1, 0, 7, 0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 0)
  l <- c(2, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 3, 0, 1, 0, 1, 0, 1, 0, 3, 0, 2, 0, 3)
  s <- smooth(roc(controls = rep(seq_along(k), k), cases = rep(seq_along(l), l), quiet = TRUE), fit = "proper")
  expect_equal(s$proper$logLik, -148.40788, tolerance = 1e-7)
  expect_equal(as.numeric(s$auc), 0.6254, tolerance = 1e-3)
})

test_that("direction and percent", {
  s <- smooth(r.s100b, fit = "proper")
  r.reversed <- roc(aSAH$outcome, -aSAH$s100b, quiet = TRUE)
  expect_identical(r.reversed$direction, ">")
  expect_equal(smooth(r.reversed, fit = "proper")$proper[c("da", "c", "logLik", "categories")], s$proper[c("da", "c", "logLik", "categories")])
  sp <- smooth(r.s100b.percent, fit = "proper")
  expect_equal(as.numeric(sp$auc), as.numeric(s$auc) * 100)
  expect_equal(sp$sensitivities, s$sensitivities * 100)
  expect_equal(sp$specificities, s$specificities * 100)
  expect_equal(as.numeric(auc(sp, partial.auc = c(100, 80))), as.numeric(auc(s, partial.auc = c(1, 0.8))) * 100)
})

test_that("roc(smooth = TRUE) passes fit, and the fit is stored for replays", {
  s <- roc(aSAH$outcome, aSAH$s100b, quiet = TRUE, smooth = TRUE, fit = "proper")
  expect_identical(s$smoothing.args$fit, "proper")
  expect_equal(s$proper, smooth(r.s100b, fit = "proper")$proper)
  expect_identical(pROC:::roc_utils_resmooth(r.s100b, s$smoothing.args)$proper, s$proper)
  # the least-squares fit stays the default, and its smoothing.args are unchanged
  expect_identical(smooth(r.s100b)$smoothing.args, smooth(r.s100b, fit = "lm")$smoothing.args)
  expect_null(smooth(r.s100b)$smoothing.args$fit)
})

test_that("Obuchowski methods refuse proper curves", {
  # their A and B parameters are those of the conventional binormal model
  s1 <- smooth(r.s100b, fit = "proper")
  s2 <- smooth(r.ndka, fit = "proper")
  msg <- "Obuchowski method is not supported for smoothed ROCs. Use method=\"bootstrap\" instead."
  expect_error(var(s1, method = "obuchowski"), msg, fixed = TRUE)
  expect_error(cov(s1, s2, method = "obuchowski"), msg, fixed = TRUE)
  expect_error(power.roc.test(s1, method = "obuchowski", power = 0.9))
})

test_that("print shows the proper fit", {
  s <- smooth(r.s100b, fit = "proper")
  # d_a is near 0, where the curve hardly depends on it (Pesce & Metz p. 820)
  expect_output(print(s), "Smoothing: proper binormal \\(d_a: 0\\.0[0-9]+; c: -0\\.38[0-9]*\\)")
})

test_that("data without curve parameters to estimate are not smoothable", {
  # Metz & Pan p. 27: perfect separation
  expect_error(smooth(roc(controls = c(1, 2, 3), cases = c(4, 5, 6), quiet = TRUE), fit = "proper"),
    "ROC curve not smoothable (perfect performance).",
    fixed = TRUE
  )
  # Points only on the top and left borders, (7/30, 1) and (0, 27/40):
  # consistent with perfect performance within the tied category. The fit
  # stayed at the chance line (logLik -76.4, the maximum is -41.5 at the
  # limit AUC = 1); before, a row of NAs from the operating points
  # failed with "missing value where TRUE/FALSE needed"
  r <- roc(controls = rep(0:1, c(23, 7)), cases = rep(1:4, c(13, 12, 10, 5)), quiet = TRUE)
  expect_error(smooth(r, fit = "proper"), "ROC curve not smoothable (perfect performance).", fixed = TRUE)
  expect_error(smooth(r), "ROC curve not smoothable (not enough points).", fixed = TRUE)
  # all values tied: a single category
  expect_error(smooth(roc(controls = c(1, 1, 1), cases = c(1, 1), quiet = TRUE), fit = "proper"),
    "ROC curve not smoothable (not enough points).",
    fixed = TRUE
  )
})

test_that("a single interior point gives a proper curve through it", {
  r <- roc(controls = c(1, 1, 1, 2, 2), cases = c(1, 2, 2, 2, 2), quiet = TRUE)
  s <- smooth(r, fit = "proper")
  expect_true(all(diff(proper_slopes(s)) <= 0))
  # the fit has the likelihood of the saturated model (the data are the
  # categories' observed proportions)
  k <- s$proper$categories$controls
  l <- s$proper$categories$cases
  expect_equal(s$proper$logLik, sum(k * log(k / sum(k))) + sum(l * log(l / sum(l))), tolerance = 1e-6)
})

test_that("bootstrap CI of a proper fit", {
  skip_slow()
  s <- smooth(r.s100b, fit = "proper")
  ci <- ci.auc(s, method = "bootstrap", boot.n = 10)
  expect_is(ci, "ci.auc")
  expect_identical(attr(ci, "method"), "bootstrap")
  expect_true(ci[1] <= ci[2] && ci[2] <= ci[3])
})

test_that("the gradient of Eq. 36 is that of the AUC", {
  h <- 1e-6
  for (p in list(c(1.24365, -0.428571), c(0.1, -0.154357), c(2.324, 0.292369), c(0.5, 0.9), c(1.46574, 0), c(0, -0.38))) {
    numerical <- c(
      (pROC:::proper_binormal_auc(p[1] + h, p[2]) - pROC:::proper_binormal_auc(max(0, p[1] - h), p[2])) / (p[1] + h - max(0, p[1] - h)),
      (pROC:::proper_binormal_auc(p[1], p[2] + h) - pROC:::proper_binormal_auc(p[1], p[2] - h)) / (2 * h)
    )
    expect_equal(pROC:::proper_binormal_auc_gradient(p[1], p[2]), numerical, tolerance = 1e-6)
  }
  # the AUC is not differentiable in c at da = c = 0 (1 + arcsin(rho) / pi)
  expect_identical(pROC:::proper_binormal_auc_gradient(0, 0)[2], NaN)
})

test_that("the information matrix is that of the multinomial likelihood", {
  # N0 sum dp dp' / p + N1 sum dq dq' / q, with the derivatives of the
  # probabilities of the categories by finite differences
  for (s in list(smooth(r.s100b, fit = "proper"), smooth(r.mtcars, fit = "proper"), smooth(r.wfns, fit = "proper"))) {
    k <- s$proper$categories$controls
    l <- s$proper$categories$cases
    theta <- c(s$proper$da, s$proper$c, s$proper$cutpoints)
    probs <- function(theta, slope, sign) {
      range <- pROC:::proper_binormal_range(theta[1], theta[2])
      v <- c(range[1], theta[-(1:2)], range[2])
      pROC:::proper_binormal_probs(pROC:::proper_binormal_args(v, theta[1], theta[2], slope(theta[2]), sign))
    }
    information <- 0
    for (group in list(list(n = sum(k), slope = function(c) 1 - c, sign = -1), list(n = sum(l), slope = function(c) 1 + c, sign = 1))) {
      p <- probs(theta, group$slope, group$sign)
      J <- vapply(seq_along(theta), function(i) {
        h <- replace(numeric(length(theta)), i, 1e-6)
        (probs(theta + h, group$slope, group$sign) - probs(theta - h, group$slope, group$sign)) / 2e-6
      }, numeric(length(p)))
      information <- information + group$n * crossprod(J[p > 0, ] / sqrt(p[p > 0]))
    }
    expect_equal(pROC:::proper_binormal_information(s$proper$cutpoints, s$proper$da, s$proper$c, k, l), information, tolerance = 1e-6)
  }
})

test_that("var and ci.auc of a proper fit use the fit by default", {
  s <- smooth(r.s100b, fit = "proper")
  expect_equal(var(s), 0.0016672558, tolerance = 1e-5)
  expect_identical(var(s, method = "binormal"), var(s))
  expect_identical(var(auc(s)), var(s))
  ci <- ci.auc(s)
  expect_identical(attr(ci, "method"), "binormal")
  expect_equal(as.numeric(ci), qnorm(c(0.025, 0.5, 0.975), as.numeric(s$auc), sqrt(var(s))))
  expect_identical(ci(s), ci)
  expect_identical(ci.auc(s, method = "binormal"), ci)
  expect_equal(as.numeric(ci.auc(s, conf.level = 0.9)), qnorm(c(0.05, 0.5, 0.95), as.numeric(s$auc), sqrt(var(s))))
  expect_output(print(ci), "95% CI: 0.6534-0.8134 (binormal)", fixed = TRUE)
  # an ordered predictor
  expect_equal(var(smooth(r.wfns, fit = "proper")), 0.0014777159, tolerance = 1e-5)
  # percent
  sp <- smooth(r.s100b.percent, fit = "proper")
  expect_equal(var(sp), var(s) * 100^2)
  expect_equal(as.numeric(ci(sp)), as.numeric(ci) * 100)
  expect_output(print(ci(sp)), "95% CI: 65.34%-81.34% (binormal)", fixed = TRUE)
  # the all-in-one call
  expect_identical(attr(roc(aSAH$outcome, aSAH$s100b, quiet = TRUE, smooth = TRUE, fit = "proper", ci = TRUE)$ci, "method"), "binormal")
})

test_that("the binormal method needs the full AUC of a proper fit", {
  msg.fit <- "Binormal method is supported only for ROC curves smoothed with method=\"binormal\" and fit=\"proper\". Use method=\"bootstrap\" instead."
  msg.partial <- "Binormal method is not supported for partial AUC. Use method=\"bootstrap\" instead."
  # the least-squares fit, which is not a maximum likelihood fit
  expect_error(var(smooth(r.s100b), method = "binormal"), msg.fit, fixed = TRUE)
  expect_error(ci.auc(smooth(r.s100b), method = "binormal"), msg.fit, fixed = TRUE)
  expect_identical(attr(ci.auc(smooth(r.s100b), boot.n = 2), "method"), "bootstrap")
  # an empirical curve
  expect_error(var(r.s100b, method = "binormal"), msg.fit, fixed = TRUE)
  # a partial AUC
  s <- smooth(r.wfns.partial, fit = "proper")
  expect_error(var(s, method = "binormal"), msg.partial, fixed = TRUE)
  expect_error(ci.auc(s, method = "binormal"), msg.partial, fixed = TRUE)
  expect_identical(attr(ci.auc(s, boot.n = 2), "method"), "bootstrap")
})

test_that("roc.test of unpaired proper fits", {
  # Dorfman & Alf (1969) Eq. 10, Pesce & Metz (2007) p. 822
  women <- aSAH$gender == "Female"
  s1 <- smooth(roc(aSAH$outcome[women], aSAH$s100b[women], quiet = TRUE), fit = "proper")
  s2 <- smooth(roc(aSAH$outcome[!women], aSAH$s100b[!women], quiet = TRUE), fit = "proper")
  t <- roc.test(s1, s2)
  expect_identical(t$method, "Binormal test for two ROC curves")
  z <- (as.numeric(s1$auc) - as.numeric(s2$auc)) / sqrt(var(s1) + var(s2))
  expect_equal(t$statistic, c(Z = z))
  expect_equal(t$p.value, 2 * pnorm(-abs(z)))
  expect_equal(roc.test(s1, s2, alternative = "less")$p.value, pnorm(z))
  expect_equal(roc.test(s1, s2, alternative = "greater")$p.value, pnorm(-z))
  expect_identical(roc.test(s1, s2, method = "binormal")$statistic, t$statistic)
  expect_equal(unname(t$estimate), c(as.numeric(s1$auc), as.numeric(s2$auc)))
  expect_identical(t$roc1$proper, s1$proper)
  expect_identical(t$roc2$proper, s2$proper)
  # percent
  s1p <- smooth(roc(aSAH$outcome[women], aSAH$s100b[women], percent = TRUE, quiet = TRUE), fit = "proper")
  s2p <- smooth(roc(aSAH$outcome[!women], aSAH$s100b[!women], percent = TRUE, quiet = TRUE), fit = "proper")
  expect_equal(roc.test(s1p, s2p)$statistic, t$statistic)
})

test_that("the binormal method refuses paired curves", {
  s1 <- smooth(r.s100b, fit = "proper")
  s2 <- smooth(r.ndka, fit = "proper")
  msg <- "Binormal method is not supported for paired ROC curves. Use method=\"bootstrap\" instead."
  expect_error(roc.test(s1, s2, method = "binormal"), msg, fixed = TRUE)
  # the user can test them as unpaired
  expect_identical(roc.test(s1, s2, paired = FALSE)$method, "Binormal test for two ROC curves")
  expect_warning(
    expect_error(roc.test(smooth(r.s100b), s2, method = "binormal", paired = FALSE),
      "Binormal method is supported only for ROC curves smoothed with method=\"binormal\" and fit=\"proper\". Use method=\"bootstrap\" instead.",
      fixed = TRUE
    ),
    "Different smoothing parameters in the ROC curves.",
    fixed = TRUE
  )
  # paired proper fits are bootstrapped by default
  skip_slow()
  expect_identical(roc.test(s1, s2, boot.n = 2)$method, "Bootstrap test for two correlated ROC curves")
})
