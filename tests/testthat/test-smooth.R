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
  # Same values with direction ">": the point reflection of the "<" curve
  r.gt <- roc(controls = controls, cases = cases, direction = ">", quiet = TRUE)
  # Negated values with direction ">": the same curve as r.lt
  r.neg <- roc(controls = -controls, cases = -cases, direction = ">", quiet = TRUE)
  for (method in c("logcondens", "logcondens.smooth")) {
    s.lt <- smooth(r.lt, method = method, n = 10)
    s.gt <- smooth(r.gt, method = method, n = 10)
    expect_identical(s.gt$specificities, s.lt$specificities)
    expect_identical(s.gt$sensitivities, 1 - rev(s.lt$sensitivities))
    expect_equal(as.numeric(s.gt$auc), 1 - as.numeric(s.lt$auc))
    expect_gt(as.numeric(smooth(r.neg, method = method, n = 10)$auc), 0.5)
  }
})

test_that("smooth with reuse.ci recomputes the CI on the smoothed curve", {
  set.seed(42)
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

test_that("fitdistr smoothing works with the t distribution", {
  skip_if_not_installed("MASS")
  r <- roc(rep(0:1, each = 30), c(qt(ppoints(30), 5), qt(ppoints(30), 5) * 2 + 2), quiet = TRUE)
  # lower bounds keep optim() in the parameter space (s > 0, df >= 1): no NaNs produced
  s <- smooth(r, method = "fitdistr", density = "t", lower = c(-Inf, 0.001, 1))
  expect_s3_class(s, "smooth.roc")
  # The density is the location-scale t that MASS::fitdistr fits
  fit <- MASS::fitdistr(r$controls, "t", lower = c(-Inf, 0.001, 1))
  expect_equal(
    sum(log(pROC:::dt_location_scale(r$controls, fit$estimate[["m"]], fit$estimate[["s"]], fit$estimate[["df"]]))),
    fit$loglik
  )
  # The exact AUC of the two fitted distributions: P(case > control)
  e0 <- s$fit.controls$estimate
  e1 <- s$fit.cases$estimate
  exact <- integrate(function(x) {
    pROC:::dt_location_scale(x, e1[["m"]], e1[["s"]], e1[["df"]]) * pROC:::pt_location_scale(x, e0[["m"]], e0[["s"]], e0[["df"]])
  }, -Inf, Inf)$value
  expect_equal(as.numeric(s$auc), exact, tolerance = 1e-4)
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
  # abs(s): optim() may try a negative scale, which would give NaNs
  dens <- function(x, m, s, df = 30) dt((x - m) / abs(s), df) / abs(s)
  expected <- smooth(r, method = "fitdistr", density = dens, start = list(m = 2, s = 1), df = 5)
  # with an optim() argument for fitdistr (maxit = 500 is the default)
  s <- smooth(r,
    method = "fitdistr", density = dens, start = list(m = 2, s = 1), df = 5,
    control = list(maxit = 500)
  )
  expect_equal(s$sensitivities, expected$sensitivities)
  expect_equal(s$specificities, expected$specificities)
})

test_that("fitdistr smoothing of a single case gives an NA curve", {
  testthat::skip_if_not_installed("MASS")
  # fitdistr estimates an NA sd for a single observation, as in 1.19.1
  r <- roc(controls = c(0.5, 0.5, 0.5, 1.5), cases = 2.5, quiet = TRUE)
  for (density in c("normal", "lognormal")) {
    s <- smooth(r, method = "fitdistr", density = density)
    expect_true(anyNA(s$sensitivities))
    expect_true(is.na(s$auc))
  }
  # Also for a single control, and one observation in each class (not an AUC of 0.5)
  for (r in list(roc(controls = 1, cases = 2, quiet = TRUE), roc(controls = 1, cases = c(2, 3, 4), quiet = TRUE))) {
    for (density in c("normal", "lognormal")) {
      s <- smooth(r, method = "fitdistr", density = density)
      expect_true(anyNA(s$sensitivities))
      expect_true(is.na(s$auc))
    }
  }
  # Constant classes have an sd of 0, not NA
  expect_equal(as.numeric(smooth(roc(controls = c(1, 1, 1), cases = c(2, 2, 2), quiet = TRUE), method = "fitdistr")$auc), 1)
  set.seed(1)
  r <- roc(controls = c(1, 1, 1), cases = rnorm(20, 3), quiet = TRUE)
  expect_equal(as.numeric(smooth(r, method = "fitdistr")$auc), 0.994585999818758, tolerance = 1e-12)
})

test_that("missing, infinite or all-zero density values give a NA/NaN AUC", {
  x <- seq(-4, 6, length.out = 200)
  dc <- dnorm(x)
  dk <- dnorm(x, 1.5)
  aucs <- function(d) {
    c(
      roc(density.controls = dc, density.cases = d)$auc,
      roc(density.controls = dc, density.cases = d, direction = "<")$auc,
      smooth(r.s100b, method = "density", density.controls = dc, density.cases = d)$auc,
      smooth(r.s100b, method = "density", density.controls = d, density.cases = dk)$auc
    )
  }
  expect_true(all(is.na(aucs(c(NA, dc[-1])))))
  # An all-negative density becomes all zero
  for (d in list(c(NaN, dc[-1]), c(Inf, dc[-1]), 0 * dc, -dc)) {
    expect_true(all(is.nan(aucs(d))))
  }
  # -Inf is a negative value: it becomes 0
  expect_identical(aucs(c(-Inf, dc[-1])), aucs(c(0, dc[-1])))
  expect_equal(aucs(c(-Inf, dc[-1])), c(0.5, 0.5, 0.5, 0.8555), tolerance = 1e-4)
})

test_that("negative density values are set to 0, as density() does", {
  d1 <- dnorm(seq(-3, 3, length.out = 512))
  d2 <- dnorm(seq(-2, 4, length.out = 512))
  with.value <- function(v) {
    d <- d2
    d[100] <- v
    d
  }
  ref <- roc(density.controls = d1, density.cases = with.value(0))
  # 1.19.1's AUC for the rounding-sized values
  for (v in c(-1e-17, -1e-14, -1e-12, -0.05)) {
    r <- roc(density.controls = d1, density.cases = with.value(v))
    expect_identical(r$direction, ref$direction)
    expect_identical(r$sensitivities, ref$sensitivities)
    expect_identical(r$specificities, ref$specificities)
    expect_identical(as.numeric(r$auc), as.numeric(ref$auc))
  }
  expect_equal(as.numeric(ref$auc), 0.7547, tolerance = 1e-4)
  # Also through smooth()
  s <- smooth(r.s100b, method = "density", density.controls = d1, density.cases = with.value(-1e-14))
  s.ref <- smooth(r.s100b, method = "density", density.controls = d1, density.cases = with.value(0))
  expect_identical(s$sensitivities, s.ref$sensitivities)
  expect_identical(s$specificities, s.ref$specificities)
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

test_that("roc(smooth = TRUE) records the name of a fitdistr density function", {
  testthat::skip_if_not_installed("MASS")
  st <- list(meanlog = -2, sdlog = 1)
  # lower keeps sdlog > 0 in optim(): no NaNs produced
  s <- roc(aSAH$outcome, aSAH$s100b,
    smooth = TRUE, smooth.method = "fitdistr",
    density = dlnorm, start.controls = st, start.cases = st, lower = c(-Inf, 0.001), quiet = TRUE
  )
  expect_identical(s$fit.controls$densfun, "dlnorm")
  expect_identical(s$fit.cases$densfun, "dlnorm")
  s <- roc(outcome ~ s100b, data = aSAH,
    smooth = TRUE, smooth.method = "fitdistr",
    density = dlnorm, start.controls = st, start.cases = st, lower = c(-Inf, 0.001), quiet = TRUE
  )
  expect_true(any(grepl("fitting dlnorm distribution for controls", capture.output(print(s)))))
})

test_that("fitdistr smoothing of a constant class is not smoothable", {
  skip_if_not_installed("MASS")
  flat <- roc(controls = c(1, 1, 1), cases = c(2, 3, 4), quiet = TRUE)
  msg <- "ROC curve not smoothable (fitdistr could not estimate the distribution parameters)."
  # optim() tries a negative scale on the constant class: dlogis() gives NaNs
  expect_warning(expect_error(smooth(flat, method = "fitdistr", density = "logistic"), msg, fixed = TRUE), "NaNs produced")
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
  # (lower: optim() warns that Nelder-Mead is unreliable in one dimension)
  s <- smooth(r.s100b, method = "fitdistr", density = dnorm, start = list(sd = 0.1), mean = 0.2, lower = 0.001)
  expect_named(s$fit.controls$estimate, "sd")
  s <- smooth(r.s100b, method = "fitdistr", density = "weibull", shape = 2, lower = 0.001)
  expect_named(s$fit.controls$estimate, "scale")
})
