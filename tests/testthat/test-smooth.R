library(pROC)
data(aSAH)

context("smooth")

# Define some density functions

unif.density <- function(x, n, from, to, bw, kernel, ...) {
  smooth.x <- seq(from = from, to = to, length.out = n)
  smooth.y <- dunif(smooth.x, min = min(x), max = max(x))
  return(smooth.y)
}

norm.density <- function(x, n, from, to, bw, kernel, ...) {
  smooth.x <- seq(from = from, to = to, length.out = n)
  smooth.y <- dnorm(smooth.x, mean = mean(x), sd = sd(x))
  return(smooth.y)
}

lnorm.density <- function(x, n, from, to, bw, kernel, ...) {
  smooth.x <- seq(from = from, to = to, length.out = n)
  smooth.y <- dlnorm(smooth.x, meanlog = mean(x), sdlog = sd(x))
  return(smooth.y)
}

test_that("We fall back to the standard smooth", {
  tukey <- smooth(c(4, 1, 3, 6, 6, 4, 1, 6, 2, 4, 2))
  expect_is(tukey, "tukeysmooth")
  expect_equal(as.numeric(tukey), c(3, 3, 3, 3, 4, 4, 4, 4, 2, 2, 2))
})

test_that("smooth with a density function works", {
  smoothed <- smooth(r.ndka, method = "density", density = unif.density, n = 10)
  expect_is(smoothed, "smooth.roc")
  expect_equal(smoothed$sensitivities, c(1, 1, 1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0, 0))
  expect_equal(smoothed$specificities, c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(as.numeric(smoothed$auc), 0.9375)
})

test_that("smooth with two density functions works", {
  smoothed <- smooth(r.ndka, method = "density", density.controls = norm.density, density.cases = lnorm.density, n = 10)
  expect_is(smoothed, "smooth.roc")
  expect_equal(smoothed$sensitivities, c(
    1, 1, 1, 0.635948942024884, 0.460070154191559, 0.344004532431686,
    0.25735248652959, 0.188201024566009, 0.130658598389315, 0.0813814489619488,
    0.0382893349015216, 0
  ))
  expect_equal(smoothed$specificities, c(0, 0, 0.832138478872629, 0.99999996787709, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(as.numeric(smoothed$auc), 0.9694449)
})


test_that("smooth with fitdistr works", {
  testthat::skip_if_not_installed("MASS")
  smoothed <- smooth(r.ndka, method = "fitdistr", n = 10)
  expect_is(smoothed, "smooth.roc")
  expect_equal(smoothed$sensitivities, c(
    1, 0.818181818181818, 0.683063106156167, 0.636363636363636,
    0.629229378820631, 0.591536044948605, 0.555151655671379, 0.51056251520327,
    0.454545454545455, 0.272727272727273, 0.0909090909090908, 0
  ))
  expect_equal(smoothed$specificities, c(
    0, 0.000240703151138483, 0.090909090909091, 0.24225342148073,
    0.272727272727273, 0.454545454545455, 0.636363636363636, 0.818181818181818,
    0.946312019857331, 0.999975065500302, 0.999999999999993, 1
  ))
  expect_equal(as.numeric(smoothed$auc), 0.580858400167636)
})

test_that("smooth with fitdistr different densities works", {
  testthat::skip_if_not_installed("MASS")
  smoothed <- smooth(r.ndka, method = "fitdistr", density.controls = "normal", density.cases = "lognormal", n = 10)
  expect_is(smoothed, "smooth.roc")
  expect_equal(smoothed$sensitivities, c(
    1, 1, 0.832447986812911, 0.818181818181818, 0.636363636363636,
    0.580641313705266, 0.454545454545455, 0.405633021041487, 0.272727272727273,
    0.267512440925326, 0.0909090909090908, 0
  ))
  expect_equal(smoothed$specificities, c(
    0, 0.090909090909091, 0.272727272727273, 0.281542026457948,
    0.407900467583373, 0.454545454545455, 0.579603793209369, 0.636363636363636,
    0.81102519760261, 0.818181818181818, 0.994940831852285, 1
  ))
  expect_equal(as.numeric(smoothed$auc), 0.567273983384952)
})

test_that("smooth with fitdistr with a density function works", {
  testthat::skip_if_not_installed("MASS")
  smoothed <- smooth(r.ndka,
    method = "fitdistr", n = 10,
    density.controls = dnorm, start.controls = list(mean = 10, sd = 10),
    density.cases = dlnorm, start = list(meanlog = 2.7, sdlog = .822)
  )
  expect_is(smoothed, "smooth.roc")
  expect_equal(smoothed$sensitivities, c(
    1, 1, 0.174065542189585, 0.0241224212514905, 0.00565553823693818,
    0.00176442417351747, 0.000654789746505889, 0.000269910020195159,
    0.000116630962648119, 4.8942161699917e-05, 1.65438472509127e-05,
    0
  ))
  expect_equal(smoothed$specificities, c(
    0, 0, 0.961730914432089, 0.999999997253745, 1, 1, 1, 1, 1,
    1, 1, 1
  ))
  expect_equal(as.numeric(smoothed$auc), 0.568359799581078)
})

test_that("logcondens smoothing respects direction", {
  testthat::skip_if_not_installed("logcondens")
  controls <- c(-1.2, -0.8, -0.5, -0.3, 0, 0.1, 0.4, 0.6, 0.9, 1.3)
  cases <- c(0.2, 0.5, 0.7, 1.0, 1.1, 1.4, 1.8, 2.1, 2.5)
  r.lt <- roc(controls = controls, cases = cases, direction = "<", quiet = TRUE)
  r.gt <- roc(controls = -controls, cases = -cases, direction = ">", quiet = TRUE)
  for (method in c("logcondens", "logcondens.smooth")) {
    s.lt <- smooth(r.lt, method = method, n = 10)
    s.gt <- smooth(r.gt, method = method, n = 10)
    expect_equal(s.gt$sensitivities, s.lt$sensitivities)
    expect_equal(s.gt$specificities, s.lt$specificities)
    expect_equal(as.numeric(s.gt$auc), as.numeric(s.lt$auc))
    expect_gt(as.numeric(s.gt$auc), 0.5)
  }
})

test_that("smooth with reuse.ci recomputes the CI on the smoothed curve", {
  r <- r.s100b
  r$ci <- ci.auc(r, method = "bootstrap", boot.n = 10, conf.level = 0.9, progress = "none")
  s <- smooth(r, reuse.ci = TRUE)
  expect_is(s, "smooth.roc")
  expect_is(s$ci, "ci.auc")
  expect_equal(attr(s$ci, "conf.level"), 0.9)
  expect_equal(attr(s$ci, "boot.n"), 10)
  expect_true(s$ci[1] <= s$ci[3])

  r$ci <- ci.se(r, specificities = c(0.5, 0.9), boot.n = 10, progress = "none")
  s <- smooth(r, reuse.ci = TRUE)
  expect_is(s$ci, "ci.se")
  expect_equal(attr(s$ci, "specificities"), c(0.5, 0.9))
})

test_that("binormal smoothing rejects curves with a single finite specificity", {
  # All finite points share sp = 0.75: lm(sp ~ se) has a slope of 0 up to rounding
  r <- roc(controls = c(1, 2, 3, 10), cases = c(5, 6, 7, 11), direction = "<", quiet = TRUE)
  expect_error(smooth(r, n = 6), "not smoothable")
  r2 <- roc(controls = c(4, 5, 6, 12), cases = c(9, 6, 12, 7), direction = "<", quiet = TRUE)
  expect_error(smooth(r2, n = 20), "not smoothable")
})

test_that("roc(smooth=TRUE, smooth.method='fitdistr') does not pass auc arguments to fitdistr", {
  testthat::skip_if_not_installed("MASS")
  response <- rep(0:1, each = 6)
  predictor <- c(0.5, 1.1, 1.6, 2.0, 2.4, 3.1, 1.9, 2.8, 3.5, 4.2, 5.0, 6.3)
  r <- roc(response, predictor, quiet = TRUE)
  expected <- auc(smooth(r, method = "fitdistr", density = "weibull"), partial.auc = c(1, .8))
  s <- roc(response, predictor,
    quiet = TRUE, smooth = TRUE, smooth.method = "fitdistr",
    density = "weibull", partial.auc = c(1, .8), boot.n = 10
  )
  expect_is(s, "smooth.roc")
  expect_equal(as.numeric(s$auc), as.numeric(expected))
})

test_that("fitdistr smoothing works with the t distribution", {
  skip_if_not_installed("MASS")
  set.seed(42)
  r <- roc(rep(0:1, each = 30), c(rt(30, 5), rt(30, 5) * 2 + 2), quiet = TRUE)
  s <- suppressWarnings(smooth(r, method = "fitdistr", density = "t"))
  expect_s3_class(s, "smooth.roc")
  fit.controls <- suppressWarnings(MASS::fitdistr(r$controls, "t"))$estimate
  x <- seq(min(r$controls), max(r$controls), length.out = 10)
  expect_equal(
    dt((x - fit.controls[["m"]]) / fit.controls[["s"]], fit.controls[["df"]]) / fit.controls[["s"]],
    pROC:::dt_location_scale(x, fit.controls[["m"]], fit.controls[["s"]], fit.controls[["df"]])
  )
})

test_that("fitdistr smoothing gives the ROC curve of the fitted distributions", {
  testthat::skip_if_not_installed("MASS")
  r <- roc(controls = c(-0.5, -0.2, 0, 0.1, 0.3, 0.6), cases = c(-2, 0.5, 1, 1.5, 3, 4.5), direction = "<", quiet = TRUE)
  s <- smooth(r, method = "fitdistr", n = 1000)
  m0 <- s$fit.controls$estimate
  m1 <- s$fit.cases$estimate
  # Binormal AUC of the fitted normal distributions
  expect_equal(as.numeric(s$auc), unname(pnorm((m1["mean"] - m0["mean"]) / sqrt(m0["sd"]^2 + m1["sd"]^2))), tolerance = 1e-5)
  expect_equal(
    coords(s, 0.5, input = "specificity", ret = "sensitivity")[1, 1],
    unname(1 - pnorm(m0["mean"], m1["mean"], m1["sd"])),
    tolerance = 1e-5
  )
  # Same with direction ">" on the negated predictor
  r2 <- roc(controls = -r$controls, cases = -r$cases, direction = ">", quiet = TRUE)
  expect_equal(as.numeric(smooth(r2, method = "fitdistr", n = 1000)$auc), as.numeric(s$auc), tolerance = 1e-6)
})

test_that("binormal smoothing of ordered predictors does not depend on unused levels", {
  resp <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1)
  lv <- c("a", "b", "c", "d", "e")
  pred <- c("a", "a", "b", "b", "c", "c", "d", "e", "b", "c", "c", "d", "d", "e", "e", "e")
  o <- roc(resp, factor(pred, levels = lv, ordered = TRUE), quiet = TRUE)
  oe <- roc(resp, factor(pred, levels = c("a", "b", "bc", "c", "d", "e"), ordered = TRUE), quiet = TRUE)
  rn <- roc(resp, as.integer(factor(pred, levels = lv)), quiet = TRUE)
  expect_equal(coef(smooth(oe)$model), coef(smooth(o)$model))
  expect_equal(coef(smooth(o)$model), coef(smooth(rn)$model))
})

test_that("fitdistr smoothing with a density function accepts other ... arguments", {
  testthat::skip_if_not_installed("MASS")
  response <- rep(0:1, each = 6)
  predictor <- c(0.5, 1.1, 1.6, 2.0, 2.4, 3.1, 1.9, 2.8, 3.5, 4.2, 5.0, 6.3)
  r <- roc(response, predictor, quiet = TRUE)
  # location-scale t density with a fixed df: only m and s are fitted
  dens <- function(x, m, s, df = 30) dt((x - m) / s, df) / s
  expected <- suppressWarnings(smooth(r, method = "fitdistr", density = dens, start = list(m = 2, s = 1), df = 5))
  # with an optim() argument for fitdistr (maxit = 500 is the default)
  s <- suppressWarnings(smooth(r,
    method = "fitdistr", density = dens, start = list(m = 2, s = 1), df = 5,
    control = list(maxit = 500)
  ))
  expect_equal(s$sensitivities, expected$sensitivities)
  expect_equal(s$specificities, expected$specificities)
  # with an auc argument forwarded by roc()
  r2 <- suppressWarnings(roc(response, predictor,
    quiet = TRUE, smooth = TRUE, smooth.method = "fitdistr",
    density = dens, start = list(m = 2, s = 1), df = 5, partial.auc = c(1, .8)
  ))
  expect_equal(as.numeric(r2$auc), as.numeric(auc(expected, partial.auc = c(1, .8))))
})

test_that("density smoothing passes only density arguments from roc() to density", {
  response <- rep(0:1, each = 10)
  predictor <- c(1, 2, 2.5, 3, 3.5, 4, 5, 5.5, 6, 8, 3, 4.5, 5, 6, 6.5, 7, 7.5, 8, 9, 10)
  r <- smooth(roc(response, predictor, quiet = TRUE), method = "density", adjust = 2)
  set.seed(42)
  expected <- ci.coords(r, x = 0.9, input = "specificity", ret = "se", boot.n = 3)
  # ci.coords' x must not reach density() as its data
  set.seed(42)
  obtained <- ci.coords(response, predictor, smooth = TRUE, smooth.method = "density", adjust = 2, x = 0.9, input = "specificity", ret = "se", boot.n = 3, quiet = TRUE)
  expect_equal(obtained$sensitivity, expected$sensitivity)
  set.seed(42)
  obtained <- roc(response, predictor, smooth = TRUE, smooth.method = "density", adjust = 2, ci = TRUE, of = "coords", x = 0.9, input = "specificity", ret = "se", boot.n = 3, quiet = TRUE)$ci
  expect_equal(obtained$sensitivity, expected$sensitivity)
  # same with a density function
  dens <- function(x, ...) stats::density(x, ...)
  s <- smooth(roc(response, predictor, quiet = TRUE), method = "density", density = dens, adjust = 2)
  s2 <- roc(response, predictor, quiet = TRUE, smooth = TRUE, smooth.method = "density", density = dens, adjust = 2, x = 0.9, of = "coords")
  expect_equal(s2$sensitivities, s$sensitivities)
  expect_equal(s$sensitivities, r$sensitivities)
})

test_that("fitdistr smoothing of a class with a single observation fails", {
  testthat::skip_if_not_installed("MASS")
  # fitdistr estimates an NA sd for a single observation
  r <- roc(controls = c(0.5, 0.5, 0.5, 1.5), cases = 2.5, quiet = TRUE)
  expect_error(smooth(r, method = "fitdistr"), "not smoothable")
  expect_error(smooth(r, method = "fitdistr", density = "lognormal"), "not smoothable")
  r <- roc(controls = 0.5, cases = c(1.5, 2.5, 3), quiet = TRUE)
  expect_error(smooth(r, method = "fitdistr"), "not smoothable")
})

test_that("non-stratified bootstraps of fitdistr smoothed curves ignore resamples with a single case or control", {
  testthat::skip_if_not_installed("MASS")
  s <- smooth(roc(c(0, 0, 1, 1, 1), c(0.5, 1.5, 1.0, 2.5, 3.1), quiet = TRUE), method = "fitdistr")
  expect_error(suppressWarnings(ci.sp(s, 0.5, boot.n = 20, boot.stratified = FALSE)), NA)
  expect_error(suppressWarnings(ci.coords(s, 0.5, input = "sensitivity", ret = "specificity", boot.n = 20, boot.stratified = FALSE)), NA)
  expect_error(suppressWarnings(roc.test(s, s, method = "sensitivity", sensitivity = 0.5, boot.n = 20, boot.stratified = FALSE)), NA)
})

test_that("roc(smooth = TRUE) records the name of a fitdistr density function", {
  testthat::skip_if_not_installed("MASS")
  st <- list(meanlog = -2, sdlog = 1)
  s <- suppressWarnings(roc(aSAH$outcome, aSAH$s100b,
    smooth = TRUE, smooth.method = "fitdistr",
    density = dlnorm, start.controls = st, start.cases = st, quiet = TRUE
  ))
  expect_identical(s$fit.controls$densfun, "dlnorm")
  expect_identical(s$fit.cases$densfun, "dlnorm")
  s <- suppressWarnings(roc(outcome ~ s100b, data = aSAH,
    smooth = TRUE, smooth.method = "fitdistr",
    density = dlnorm, start.controls = st, start.cases = st, quiet = TRUE
  ))
  expect_true(any(grepl("fitting dlnorm distribution for controls", capture.output(print(s)))))
})

test_that("roc(smooth = TRUE) keeps only the smoothing arguments for the replays", {
  pdf(NULL)
  on.exit(dev.off())
  title <- quote(italic(S100) * beta)
  s <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, plot = TRUE, main = title, xlab = quote(alpha), col = "red", quiet = TRUE)
  s.char <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, quiet = TRUE)
  expect_identical(s$smoothing.args, s.char$smoothing.args)
  expect_false(any(c("main", "xlab", "col") %in% names(s$smoothing.args)))
  s2 <- roc(aSAH$outcome, aSAH$ndka, smooth = TRUE, quiet = TRUE)
  seed <- sample.int(1e6, 1)
  set.seed(seed)
  expected <- roc.test(s.char, s2, boot.n = 5, progress = "none")$p.value
  set.seed(seed)
  expect_equal(roc.test(s, s2, boot.n = 5, progress = "none")$p.value, expected)
  set.seed(seed)
  expected <- as.numeric(ci.auc(s.char, boot.n = 5, progress = "none"))
  set.seed(seed)
  expect_equal(as.numeric(ci.auc(s, boot.n = 5, progress = "none")), expected)
  set.seed(seed)
  ci.plotmath <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, ci = TRUE, boot.n = 5, progress = "none", plot = TRUE, main = title, quiet = TRUE)$ci
  expect_equal(as.numeric(ci.plotmath), expected)
})

test_that("the graphical adj argument does not change density smoothing", {
  pdf(NULL)
  on.exit(dev.off())
  s <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, smooth.method = "density", quiet = TRUE)
  s.adj <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, smooth.method = "density", plot = TRUE, main = "S100b", adj = 0.5, quiet = TRUE)
  expect_equal(s.adj$sensitivities, s$sensitivities)
  expect_equal(as.numeric(s.adj$auc), as.numeric(s$auc))
  s.adj <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, smooth.method = "density", plot = TRUE, main = "S100b", adj = 0, quiet = TRUE)
  expect_equal(as.numeric(s.adj$auc), as.numeric(s$auc))
  s.adj <- plot.roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, smooth.method = "density", adj = 0.5, quiet = TRUE)
  expect_equal(as.numeric(s.adj$auc), as.numeric(s$auc))
})

test_that("invalid density values give a clear error", {
  x <- seq(-4, 6, length.out = 200)
  dc <- dnorm(x)
  dk <- dnorm(x, 1.5)
  msg <- "'density.controls' and 'density.cases' must be finite, non-negative and not all zero"
  # tiny negative values, like those of FFT-based estimators, also make the AUC wrong
  tiny.negative <- ifelse(dc < 1e-6, -1e-17, dc)
  invalid <- list(c(NA, dc[-1]), c(NaN, dc[-1]), c(Inf, dc[-1]), 0 * dc, dc - 0.01, tiny.negative)
  for (d in invalid) {
    expect_error(roc(density.controls = d, density.cases = dk), msg)
    expect_error(roc(density.controls = dc, density.cases = d), msg)
    expect_error(roc(density.controls = dc, density.cases = d, direction = "<"), msg)
    expect_error(smooth(r.s100b, method = "density", density.controls = dc, density.cases = d), msg)
    expect_error(smooth(r.s100b, method = "density", density.controls = d, density.cases = dk), msg)
  }
})

test_that("a single density value is rejected", {
  msg <- "'density.controls' and 'density.cases' must contain at least two density values."
  expect_error(smooth(r.s100b, method = "density", density = 3), msg, fixed = TRUE)
  expect_error(smooth(r.s100b, method = "density", density.controls = 0.2, density.cases = 0.5), msg, fixed = TRUE)
  expect_error(roc(density.controls = 0.2, density.cases = 0.5), msg, fixed = TRUE)
  # two values are still a (coarse) density
  expect_s3_class(roc(density.controls = c(0.8, 0.2), density.cases = c(0.2, 0.8)), "smooth.roc")
})

test_that("smooth with reuse.ci skips a stored ci.thresholds", {
  r <- r.s100b
  r$ci <- ci.thresholds(r, thresholds = 0.5, boot.n = 2)
  expect_warning(s <- smooth(r, reuse.ci = TRUE), "not available for smoothed ROC curves", fixed = TRUE)
  expect_s3_class(s, "smooth.roc")
  expect_null(s$ci)
  expect_equal(as.numeric(s$auc), as.numeric(smooth(r.s100b)$auc))
})

test_that("smooth with reuse.ci skips a stored ci.coords of thresholds", {
  r <- r.s100b
  # default x = "best": input = "threshold"
  r$ci <- ci.coords(r, "best", best.policy = "random", boot.n = 2, progress = "none")
  expect_warning(s <- smooth(r, reuse.ci = TRUE), "not available for smoothed ROC curves", fixed = TRUE)
  expect_s3_class(s, "smooth.roc")
  expect_null(s$ci)
  r$ci <- ci.coords(r, 0.2, input = "threshold", ret = c("se", "sp"), boot.n = 2, progress = "none")
  expect_warning(s <- smooth(r, reuse.ci = TRUE), "not available for smoothed ROC curves", fixed = TRUE)
  expect_null(s$ci)
  # other inputs are still recomputed
  r$ci <- ci.coords(r, 0.9, input = "specificity", ret = "se", boot.n = 2, progress = "none")
  s <- smooth(r, reuse.ci = TRUE)
  expect_s3_class(s$ci, "ci.coords")
})

test_that("smooth with reuse.ci recomputes a stored ci.coords of the best point", {
  r <- r.s100b
  # x = "best": the default input = "threshold" is irrelevant
  r$ci <- ci.coords(r, "best", ret = c("specificity", "sensitivity"), best.policy = "random", boot.n = 2, progress = "none")
  expect_identical(attr(r$ci, "input"), "threshold")
  warnings <- character(0)
  s <- withCallingHandlers(smooth(r, reuse.ci = TRUE), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("not re-used", warnings, fixed = TRUE)))
  expect_s3_class(s$ci, "ci.coords")
  expect_identical(attr(s$ci, "x"), "best")
  expect_identical(names(s$ci), c("specificity", "sensitivity"))
})

test_that("logcondens smoothing of a class with a single value is not smoothable", {
  skip_if_not_installed("logcondens")
  one <- roc(controls = 1, cases = c(2, 3, 4), quiet = TRUE)
  flat <- roc(controls = c(1, 2, 3), cases = c(4, 4, 4), quiet = TRUE)
  msg <- "ROC curve not smoothable (logcondens needs at least two distinct values in each class)."
  for (r in list(one, flat)) {
    expect_error(smooth(r, method = "logcondens"), msg, fixed = TRUE)
    expect_error(smooth(r, method = "logcondens.smooth"), msg, fixed = TRUE)
  }
})

test_that("fitdistr smoothing of a constant class is not smoothable", {
  skip_if_not_installed("MASS")
  flat <- roc(controls = c(1, 1, 1), cases = c(2, 3, 4), quiet = TRUE)
  msg <- "ROC curve not smoothable (fitdistr could not estimate the distribution parameters)."
  expect_error(suppressWarnings(smooth(flat, method = "fitdistr", density = "logistic")), msg, fixed = TRUE)
  expect_error(smooth(flat, method = "fitdistr", density = "weibull"), msg, fixed = TRUE)
  # Errors unrelated to a constant class are kept
  expect_error(smooth(r.s100b, method = "fitdistr", density = dnorm), "'start' must be a named list")
})

test_that("fitdistr smoothing rejects fixed parameters of closed-form densities", {
  skip_if_not_installed("MASS")
  # fitdistr estimates all the parameters of these densities, ignoring the fixed ones
  expect_error(smooth(r.s100b, method = "fitdistr", density = "normal", mean = 0.2), "'mean' cannot be held fixed with density \"normal\"")
  expect_error(smooth(r.s100b, method = "fitdistr", density.controls = "normal", density.cases = "lognormal", sdlog = 1), "'sdlog' cannot be held fixed with density \"lognormal\"")
  expect_error(smooth(r.s100b, method = "fitdistr", density = "exponential", rate = 2), "'rate' cannot be held fixed with density \"exponential\"")
  # Parameters held fixed by fitdistr are still accepted
  s <- suppressWarnings(smooth(r.s100b, method = "fitdistr", density = dnorm, start = list(sd = 0.1), mean = 0.2))
  expect_named(s$fit.controls$estimate, "sd")
  s <- suppressWarnings(smooth(r.s100b, method = "fitdistr", density = "weibull", shape = 2))
  expect_named(s$fit.controls$estimate, "scale")
})
