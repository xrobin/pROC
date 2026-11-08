library(pROC)
data(aSAH)

# Expected information n p' p'^T / p of the categories, with the
# derivatives p' of their probabilities 'probs(theta)' by central differences
numerical_information <- function(probs, theta, n) {
  h <- 1e-6
  jacobian <- vapply(seq_along(theta), function(i) {
    e <- replace(numeric(length(theta)), i, h)
    (probs(theta + e) - probs(theta - e)) / (2 * h)
  }, numeric(length(n)))
  crossprod(jacobian * sqrt(n / probs(theta)))
}

s.s100b <- smooth(r.s100b)

test_that("the ML variance of a binormal fit is O&M Eq. 6 with the inverse expected information", {
  fit <- s.s100b$binormal
  k <- fit$categories$controls
  l <- fit$categories$cases
  probs <- function(theta) unlist(pROC:::binormal_probs(theta[1], theta[2], theta[-(1:2)]))
  info <- numerical_information(probs, c(fit$a, fit$b, fit$cutpoints), rep(c(sum(k), sum(l)), each = length(k)))
  # gradient of the AUC pnorm(a / sqrt(1 + b^2)) in a and b
  gradient <- c(
    dnorm(fit$a / sqrt(1 + fit$b^2)) / sqrt(1 + fit$b^2),
    -dnorm(fit$a / sqrt(1 + fit$b^2)) * fit$a * fit$b / (1 + fit$b^2)^1.5,
    numeric(length(fit$cutpoints))
  )
  expect_equal(var(s.s100b), drop(gradient %*% solve(info, gradient)), tolerance = 1e-6)
  expect_equal(var(s.s100b), 0.0026694004298121847)
})

test_that("var and ci.auc of binormal and proper fits default to the ML method", {
  expect_identical(var(s.s100b), var(s.s100b, method = "ml"))
  expect_identical(var(s.s100b), var(auc(s.s100b)))
  ci <- ci.auc(s.s100b)
  expect_identical(attr(ci, "method"), "ml")
  expect_identical(ci, ci(s.s100b))
  expect_identical(ci, ci.auc(s.s100b, method = "ml"))
  # Wald interval around the AUC
  expect_equal(as.numeric(ci), as.numeric(auc(s.s100b)) + c(-1, 0, 1) * qnorm(0.975) * sqrt(var(s.s100b)))
  expect_equal(as.numeric(ci.auc(s.s100b, conf.level = 0.9)), as.numeric(auc(s.s100b)) + c(-1, 0, 1) * qnorm(0.95) * sqrt(var(s.s100b)))
  expect_output(print(ci), "95% CI: 0.632-0.8345 (maximum likelihood)", fixed = TRUE)
  # from the shortcuts
  for (shortcut in list(
    ci.auc(aSAH$outcome, aSAH$s100b, smooth = TRUE, quiet = TRUE),
    roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, ci = TRUE, quiet = TRUE)$ci
  )) {
    expect_identical(as.numeric(shortcut), as.numeric(ci))
    expect_identical(attr(shortcut, "method"), "ml")
  }
})

test_that("ML var and ci.auc of percent curves are rescaled", {
  s <- smooth(r.s100b.percent)
  expect_equal(var(s), var(s.s100b) * 100^2)
  expect_equal(as.numeric(ci.auc(s)), as.numeric(ci.auc(s.s100b)) * 100)
})

test_that("the ML method is refused where it does not apply", {
  msg <- "The ML method is only supported for the full AUC of curves smoothed with method=\"binormal\" or \"proper\". Use method=\"bootstrap\" instead."
  expect_error(var(r.s100b, method = "ml"), msg, fixed = TRUE)
  # partial AUC: bootstrap by default
  s.partial <- smooth(r.s100b.partial)
  expect_error(var(s.partial, method = "ml"), msg, fixed = TRUE)
  expect_error(ci.auc(s.partial, method = "ml"), msg, fixed = TRUE)
  expect_identical(attr(ci.auc(s.partial, boot.n = 2, progress = "none"), "method"), "bootstrap")
  # other smoothings
  s.density <- smooth(r.s100b, method = "density")
  expect_error(var(s.density, method = "ml"), msg, fixed = TRUE)
  expect_error(ci.auc(s.density, method = "ml"), msg, fixed = TRUE)
  expect_identical(attr(ci.auc(s.density, boot.n = 2, progress = "none"), "method"), "bootstrap")
})

test_that("binormal curves smoothed by pROC 1.19 keep the bootstrap", {
  # they hold the least-squares model, not the fit
  legacy <- function(s) {
    s$binormal <- NULL
    s$model <- lm(sp ~ se, data.frame(sp = c(1, 0.5, -0.2), se = c(-1, 0.3, 1)))
    s
  }
  s <- legacy(s.s100b)
  expect_identical(attr(ci.auc(s, boot.n = 2, progress = "none"), "method"), "bootstrap")
  expect_error(var(s, method = "ml"), "The ML method is only supported", fixed = TRUE)
  male <- aSAH$gender == "Male"
  s1 <- legacy(smooth(roc(aSAH$outcome[male], aSAH$s100b[male], quiet = TRUE)))
  s2 <- smooth(roc(aSAH$outcome[!male], aSAH$s100b[!male], quiet = TRUE))
  expect_identical(roc.test(s1, s2, boot.n = 2, progress = "none")$method, "Bootstrap test for two ROC curves")
  expect_error(roc.test(s1, s2, method = "ml"), "The ML method is only supported", fixed = TRUE)
})

test_that("unpaired roc.test of binormal fits defaults to the delta method", {
  male <- aSAH$gender == "Male"
  s1 <- smooth(roc(aSAH$outcome[male], aSAH$s100b[male], quiet = TRUE))
  s2 <- smooth(roc(aSAH$outcome[!male], aSAH$s100b[!male], quiet = TRUE))
  t <- roc.test(s1, s2)
  expect_identical(t$method, "Maximum-likelihood test for two ROC curves")
  z <- (as.numeric(auc(s1)) - as.numeric(auc(s2))) / sqrt(var(s1) + var(s2))
  expect_identical(t$statistic, c(Z = z))
  expect_equal(t$p.value, 2 * pnorm(-abs(z)))
  expect_equal(roc.test(s1, s2, alternative = "greater")$p.value, pnorm(-z))
  expect_equal(roc.test(s1, s2, alternative = "less")$p.value, pnorm(z))
  expect_identical(roc.test(s1, s2, method = "ml")[c("statistic", "p.value", "method")], t[c("statistic", "p.value", "method")])
  expect_identical(t$roc1$binormal, s1$binormal)
  expect_identical(t$roc2$binormal, s2$binormal)
  expect_identical(unname(t$estimate), c(as.numeric(auc(s1)), as.numeric(auc(s2))))
  # Percent curves
  s1.percent <- smooth(roc(aSAH$outcome[male], aSAH$s100b[male], percent = TRUE, quiet = TRUE))
  s2.percent <- smooth(roc(aSAH$outcome[!male], aSAH$s100b[!male], percent = TRUE, quiet = TRUE))
  expect_equal(roc.test(s1.percent, s2.percent)$statistic, t$statistic)
  # One curve not smoothed: bootstrap, and ML is refused
  r2 <- roc(aSAH$outcome[!male], aSAH$s100b[!male], quiet = TRUE)
  expect_warning(
    expect_error(roc.test(s1, r2, method = "ml"), "The ML method is only supported", fixed = TRUE),
    "Different smoothing parameters in the ROC curves.", fixed = TRUE
  )
})

test_that("paired roc.test refuses the ML method", {
  s.ndka <- smooth(r.ndka)
  expect_error(roc.test(s.s100b, s.ndka, method = "ml"),
    "The ML method is not supported for paired ROC curves. Use method=\"bootstrap\" instead.",
    fixed = TRUE
  )
  expect_identical(roc.test(s.s100b, s.ndka, boot.n = 2, progress = "none")$method, "Bootstrap test for two correlated ROC curves")
  # unless the user says the curves are unpaired
  expect_identical(roc.test(s.s100b, s.ndka, paired = FALSE)$method, "Maximum-likelihood test for two ROC curves")
})

# Probabilities of the categories of the proper binormal model at theta = c(da, c, cut-points)
proper_probs <- function(theta) {
  range <- pROC:::proper_binormal_range(theta[1], theta[2])
  v <- c(range[1], theta[-(1:2)], range[2])
  c(
    pROC:::proper_binormal_probs(pROC:::proper_binormal_args(v, theta[1], theta[2], 1 - theta[2], -1)),
    pROC:::proper_binormal_probs(pROC:::proper_binormal_args(v, theta[1], theta[2], 1 + theta[2], 1))
  )
}
proper_auc_gradient <- function(da, c, n) {
  h <- 1e-6
  c(
    (pROC:::proper_binormal_auc(da + h, c) - pROC:::proper_binormal_auc(da - h, c)) / (2 * h),
    (pROC:::proper_binormal_auc(da, c + h) - pROC:::proper_binormal_auc(da, c - h)) / (2 * h),
    numeric(n - 2)
  )
}

test_that("the ML variance of a proper fit uses the inverse expected information", {
  skip_if_not_installed("mvtnorm")
  s <- smooth(r.wfns, method = "proper")
  fit <- s$proper
  k <- fit$categories$controls
  l <- fit$categories$cases
  theta <- c(fit$da, fit$c, fit$cutpoints)
  info <- pROC:::proper_binormal_information(fit$da, fit$c, fit$cutpoints, k, l)
  expect_equal(info, numerical_information(proper_probs, theta, rep(c(sum(k), sum(l)), each = length(k))), tolerance = 1e-6)
  gradient <- proper_auc_gradient(fit$da, fit$c, length(theta))
  expect_equal(var(s), drop(gradient %*% solve(info, gradient)), tolerance = 1e-8)
  expect_identical(attr(ci.auc(s), "method"), "ml")
})

test_that("the proper ML variance of a single point leaves out the singular direction, with a warning", {
  skip_if_not_installed("mvtnorm")
  # two categories, three parameters: the information is singular
  expect_warning(
    s <- smooth(roc(controls = c(1, 1, 1, 2, 2), cases = c(1, 2, 2, 2, 2), quiet = TRUE), method = "proper"),
    "A proper curve through a single point is arbitrary and can be misleading.", fixed = TRUE
  )
  fit <- s$proper
  e <- eigen(pROC:::proper_binormal_information(fit$da, fit$c, fit$cutpoints, fit$categories$controls, fit$categories$cases), symmetric = TRUE)
  expect_lt(e$values[3], 1e-6 * e$values[2])
  gradient <- proper_auc_gradient(fit$da, fit$c, 3)
  msg <- "The ML variance of a proper curve through a single point is unreliable and can be misleading."
  expect_warning(v <- var(s), msg, fixed = TRUE)
  expect_equal(v, sum(crossprod(e$vectors[, 1:2], gradient)^2 / e$values[1:2]))
  expect_true(is.finite(v) && v > 0)
  expect_warning(ci.auc(s), msg, fixed = TRUE)
  # two points determine the curve: no warning
  r2 <- roc(controls = c(1, 1, 1, 2, 2, 3), cases = c(1, 2, 2, 3, 3, 3), quiet = TRUE)
  expect_no_warning(s2 <- smooth(r2, method = "proper"))
  expect_no_warning(var(s2))
})

test_that("the proper ML variance of perfect performance is 0", {
  skip_if_not_installed("mvtnorm")
  expect_warning(s <- smooth(roc(controls = c(1, 2, 3), cases = c(4, 5, 6), quiet = TRUE), method = "proper"), "Perfect performance")
  expect_warning(expect_identical(var(s), 0), "var() of a ROC curve with AUC == 1 is always 0 and can be misleading.", fixed = TRUE)
  expect_warning(ci <- ci.auc(s), "ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.", fixed = TRUE)
  expect_identical(as.numeric(ci), c(1, 1, 1))
})

test_that("unpaired roc.test of two perfect proper fits has Z = 0", {
  skip_if_not_installed("mvtnorm")
  expect_warning(s1 <- smooth(roc(controls = c(1, 2, 3), cases = c(4, 5, 6), quiet = TRUE), method = "proper"), "Perfect performance")
  expect_warning(s2 <- smooth(roc(controls = c(1, 2), cases = c(4, 5, 6, 7), quiet = TRUE), method = "proper"), "Perfect performance")
  expect_warning(t <- roc.test(s1, s2), "roc.test() of two ROC curves with AUC == 1 has always p.value = 1 and can be misleading.", fixed = TRUE)
  expect_identical(t$statistic, c(Z = 0))
  expect_identical(t$p.value, 1)
})

test_that("unpaired roc.test of proper fits defaults to the delta method", {
  skip_if_not_installed("mvtnorm")
  male <- aSAH$gender == "Male"
  s1 <- smooth(roc(aSAH$outcome[male], aSAH$wfns[male], quiet = TRUE), method = "proper")
  s2 <- smooth(roc(aSAH$outcome[!male], aSAH$wfns[!male], quiet = TRUE), method = "proper")
  t <- roc.test(s1, s2)
  expect_identical(t$method, "Maximum-likelihood test for two ROC curves")
  expect_identical(t$statistic, c(Z = (as.numeric(auc(s1)) - as.numeric(auc(s2))) / sqrt(var(s1) + var(s2))))
})

test_that("the proper ML variance warns at the chance line", {
  skip_if_not_installed("mvtnorm")
  # empirical AUC below 0.5: the fit is the chance line
  s <- smooth(roc(controls = c(2, 4, 6, 8), cases = c(1, 3, 5, 7), direction = "<", quiet = TRUE), method = "proper")
  msg <- "The ML variance of a proper curve with AUC == 0.5 is unstable and can be misleading."
  expect_warning(v <- var(s), msg, fixed = TRUE)
  expect_true(is.finite(v))
  expect_warning(ci.auc(s), msg, fixed = TRUE)
  # not at a binormal fit of the same data, nor at a proper fit above the chance line
  expect_no_warning(var(smooth(roc(controls = c(2, 4, 6, 8), cases = c(1, 3, 5, 7), direction = "<", quiet = TRUE))))
  expect_no_warning(var(smooth(r.wfns, method = "proper")))
})
