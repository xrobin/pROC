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
    1, 1, 0.65584212882921, 0.303849532306639, 0.0922807400203477,
    0.017547821937714, 0.00203415264061833, 0.000141550295211778,
    5.86072275643637e-06, 1.43622216786009e-07, 2.05997195401133e-09,
    0
  ))
  expect_equal(smoothed$specificities, c(
    0, 0, 0.961731211013412, 0.999999997253703, 1, 1, 1, 1, 1,
    1, 1, 1
  ))
  expect_equal(as.numeric(smoothed$auc), 0.814600645965216)
})

test_that("smooth with fitdistr different densities works", {
  testthat::skip_if_not_installed("MASS")
  smoothed <- smooth(r.ndka, method = "fitdistr", density.controls = "normal", density.cases = "lognormal", n = 10)
  expect_is(smoothed, "smooth.roc")
  expect_equal(smoothed$sensitivities, c(
    1, 1, 0.174065394158716, 0.0241224684680268, 0.00565556180305715,
    0.0017644346804079, 0.000654794610631603, 0.000269912354252342,
    0.000116632088037343, 4.89426737202444e-05, 1.6544031070368e-05,
    0
  ))
  expect_equal(smoothed$specificities, c(
    0, 0, 0.961731211013412, 0.999999997253703, 1, 1, 1, 1, 1,
    1, 1, 1
  ))
  expect_equal(as.numeric(smoothed$auc), 0.568359871182632)
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
