library(pROC)
data(aSAH)

context("print")

test_that("print.auc works", {
  expect_output(print(auc(r.wfns)), "^Area under the curve: 0.8237$")
  expect_output(print(auc(r.ndka.percent)), "^Area under the curve: 61.2%$")

  expect_output(print(r.ndka.partial1$auc), "^Partial area under the curve \\(specificity 0\\.99-0\\.9\\): 0\\.01046$")
  expect_output(print(r.s100b.percent.partial1$auc), "^Partial area under the curve \\(specificity 99%-90%\\): 2.983%$")
  expect_output(print(r.s100b.partial2$auc), "^Partial area under the curve \\(sensitivity 0.99-0.9\\): 0.01376$")
})

test_that("print.roc works", {
  expect_known_output(print(r.wfns), "print_output/r.wfns")
  expect_known_output(print(r.ndka), "print_output/r.ndka")
  expect_known_output(print(r.s100b), "print_output/r.s100b")
  expect_known_output(print(r.wfns.percent), "print_output/r.wfns.percent")
  expect_known_output(print(r.ndka.percent), "print_output/r.ndka.percent")
  expect_known_output(print(r.s100b.percent), "print_output/r.s100b.percent")
  expect_known_output(print(r.wfns.partial1), "print_output/r.wfns.partial1")
  expect_known_output(print(r.ndka.partial1), "print_output/r.ndka.partial1")
  expect_known_output(print(r.s100b.partial1), "print_output/r.s100b.partial1")
  expect_known_output(print(r.wfns.percent.partial1), "print_output/r.wfns.percent.partial1")
  expect_known_output(print(r.ndka.percent.partial1), "print_output/r.ndka.percent.partial1")
  expect_known_output(print(r.s100b.percent.partial1), "print_output/r.s100b.percent.partial1")
  expect_known_output(print(r.s100b.partial2), "print_output/r.s100b.partial2")
  expect_known_output(print(roc(outcome ~ ndka, aSAH)), "print_output/ndka_formula")
})

test_that("print.multiclass.roc works", {
  expect_warning(expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka)), "print_output/multiclass"), "2")
  expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, levels = c(3, 4, 5))), "print_output/multiclass_levels")
  expect_warning(expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, percent = TRUE)), "print_output/multiclass_percent"), "2")
})

test_that("print.multiclass.roc works", {
  expect_warning(expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka)), "print_output/multiclass"), "2")
  expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, levels = c(3, 4, 5))), "print_output/multiclass_levels")
  expect_warning(expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, percent = TRUE)), "print_output/multiclass_percent"), "2")
  expect_warning(expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, partial.auc = c(1, .9))), "print_output/multiclass_partial"), "2")
  expect_warning(expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$ndka, partial.auc = c(1, .9), partial.auc.focus = "se")), "print_output/multiclass_partial_se"), "2")
  expect_warning(expect_known_output(print(multiclass.roc(aSAH$gos6, aSAH$wfns, partial.auc = c(1, .9), partial.auc.correct = TRUE)), "print_output/multiclass_partial_correct"), "2")
})

test_that("print.multiclass_roc_multivariate works", {
  n <- c(100, 80, 150)
  responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
  set.seed(42)
  preds <- lapply(n, function(x) runif(x, 0.4, 0.6))
  predictor <- as.matrix(data.frame(
    "X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.7)),
    "X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0.2, 0.8)),
    "X3" = c(runif(n[1] + n[2], 0.3, 0.7), preds[[3]])
  ))
  expect_known_output(print(multiclass.roc(responses, predictor)), "print_output/mv_multiclass")

  expect_warning(expect_known_output(print(multiclass.roc(responses, predictor, levels = c("X2", "X3"))), "print_output/mv_multiclass_levels"), NA)
  expect_known_output(print(multiclass.roc(responses, predictor, percent = TRUE)), "print_output/mv_multiclass_percent")
  expect_known_output(print(multiclass.roc(responses, predictor, partial.auc = c(1, .9))), "print_output/mv_multiclass_partial")
  expect_known_output(print(multiclass.roc(responses, predictor, partial.auc = c(1, .9), partial.auc.focus = "se")), "print_output/mv_multiclass_partial_se")
  expect_known_output(print(multiclass.roc(responses, predictor, partial.auc = c(1, .9), partial.auc.correct = TRUE)), "print_output/mv_multiclass_partial_correct")
})

test_that("print works with a formula", {
  expect_known_output(print(roc(outcome ~ ndka, aSAH)), "print_output/r.ndka.formula")
  expect_warning(expect_known_output(print(multiclass.roc(gos6 ~ ndka, aSAH)), "print_output/mv_multiclass.ndka.formula"), "2")
})

test_that("print works without the auc", {
  expect_known_output(print(roc(outcome ~ ndka, aSAH, auc = FALSE)), "print_output/r.ndka.formula.no_auc")
})

test_that("print works with the CI", {
  skip_slow()
  if (getRversion() > "3.6.0") {
    suppressWarnings(RNGkind(sample.kind = "Rounding"))
  }
  set.seed(42) # For reproducible CI
  expect_known_output(print(roc(outcome ~ ndka, aSAH, ci = TRUE)), "print_output/r.ndka.formula.ci")
})

test_that("print.smooth.roc works", {
  expect_known_output(print(smooth(roc(aSAH$outcome, aSAH$ndka))), "print_output/smooth.ndka")
  expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH))), "print_output/smooth.s100b.formula")
  expect_known_output(print(smooth(roc(aSAH$outcome, aSAH$ndka))), "print_output/smooth.wfns")
  expect_known_output(print(smooth(roc(aSAH$outcome, aSAH$ndka), method = "binormal")), "print_output/smooth.s100b.binormal")
  expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH), method = "density")), "print_output/smooth.s100b.density")
})

test_that("print.smooth.roc with fitdistr works", {
  testthat::skip_if_not_installed("MASS")
  expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH), method = "fitdistr")), "print_output/smooth.s100b.fitdistr")
})

test_that("print.smooth.roc with logcondens works", {
  testthat::skip_if_not_installed("logcondens")
  expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH), method = "logcondens")), "print_output/smooth.s100b.logcondens")
  expect_known_output(print(smooth(roc(outcome ~ s100b, aSAH), method = "logcondens.smooth")), "print_output/smooth.s100b.logcondens.smooth")
})

test_that("print works with ci.auc", {
  skip_slow()
  skip_if(getRversion() < "3.6.0") # added sample.kind
  RNGkind(sample.kind = "Rejection")

  set.seed(42) # For reproducible CI
  expect_known_output(print(ci.auc(r.ndka, method = "bootstrap", boot.n = 3)), "print_output/r.ndka.ci.auc")
})

test_that("print works with ci.coords", {
  skip_slow()
  skip_if(getRversion() < "3.6.0") # added sample.kind
  RNGkind(sample.kind = "Rejection")

  set.seed(42) # For reproducible CI
  expect_known_output(print(ci.coords(r.ndka, x = c(0.5, 0.2), boot.n = 3)), "print_output/r.ndka.ci.coords")
})

test_that("print works with ci.thresholds", {
  skip_slow()
  skip_if(getRversion() < "3.6.0") # added sample.kind
  RNGkind(sample.kind = "Rejection")

  set.seed(42) # For reproducible CI
  expect_known_output(print(ci.thresholds(r.ndka, thresholds = c(0.5, 0.2), boot.n = 3)), "print_output/r.ndka.ci.thresholds")
})

test_that("print works with ci.se", {
  skip_slow()
  skip_if(getRversion() < "3.6.0") # added sample.kind
  RNGkind(sample.kind = "Rejection")

  set.seed(42) # For reproducible CI
  expect_known_output(print(ci.se(r.ndka, boot.n = 3)), "print_output/r.ndka.ci.se")
})

test_that("print works with ci.sp", {
  skip_slow()
  skip_if(getRversion() < "3.6.0") # added sample.kind
  RNGkind(sample.kind = "Rejection")

  set.seed(42) # For reproducible CI
  expect_known_output(print(ci.sp(r.ndka, boot.n = 3)), "print_output/r.ndka.ci.sp")
})

test_that("print works with a formula passed as variable", {
  x <- outcome ~ ndka
  expect_known_output(print(roc(x, aSAH)), "print_output/ndka_formula_var")
})

test_that("print works with a formula with data attached", {
  attach(aSAH)
  x <- outcome ~ ndka
  expect_known_output(print(roc(x)), "print_output/ndka_formula_var_attached")
  expect_known_output(print(roc(outcome ~ ndka)), "print_output/ndka_formula_attached")
  detach(aSAH)
})

test_that("print works with a formula with data attached with 'with'", {
  x <- outcome ~ ndka
  expect_known_output(with(aSAH, print(roc(x))), "print_output/ndka_formula_var_attached")
  expect_known_output(with(aSAH, print(roc(outcome ~ ndka))), "print_output/ndka_formula_attached")
})

test_that("smoothed curves from the formula and data.frame interfaces print the data names", {
  d <- data.frame(y = c(0, 0, 0, 0, 1, 1, 1, 1), x = c(1, 3, 5, 7, 2, 4, 6, 8), z = c(2, 1, 4, 3, 6, 5, 8, 7))
  dataline <- "Data: x in 4 controls (y 0) < 4 cases (y 1)."
  expect_output(print(roc(y ~ x, d, smooth = TRUE, smooth.method = "density", quiet = TRUE)), dataline, fixed = TRUE)
  expect_output(print(roc(y ~ x + z, d, smooth = TRUE, smooth.method = "density", quiet = TRUE)$x), dataline, fixed = TRUE)
  expect_output(print(roc(d, y, x, smooth = TRUE, smooth.method = "density", quiet = TRUE)), dataline, fixed = TRUE)
  expect_output(print(roc_(d, "y", "x", smooth = TRUE, smooth.method = "density", quiet = TRUE)), dataline, fixed = TRUE)
})

test_that("smoothed curves from the plot.roc and lines.roc formula methods print the variable names", {
  pdf(NULL)
  on.exit(dev.off())
  r <- plot.roc(outcome ~ s100b, aSAH, smooth = TRUE, quiet = TRUE)
  expect_output(print(r), "Data: s100b in 72 controls (outcome Good) < 41 cases (outcome Poor).", fixed = TRUE)
  r <- lines.roc(outcome ~ s100b, aSAH, smooth = TRUE, quiet = TRUE)
  expect_output(print(r), "Data: s100b in 72 controls (outcome Good) < 41 cases (outcome Poor).", fixed = TRUE)
})

test_that("print.mv.multiclass.roc separates formula predictor names", {
  d <- data.frame(
    resp = factor(rep(c("x", "y", "z"), 4)),
    x = c(1, 5, 9, 2, 6, 7, 3, 4, 8, 2, 6, 9),
    y = c(9, 5, 1, 8, 4, 2, 7, 6, 3, 8, 5, 1),
    z = c(2, 3, 9, 1, 4, 8, 3, 2, 7, 1, 5, 9)
  )
  expect_output(print(multiclass.roc(resp ~ x + y + z, d)), "Data: multivariate predictor x, y, z with 3 levels of resp: x, y, z.", fixed = TRUE)
  expect_output(print(multiclass.roc(resp ~ ., d)), "Data: multivariate predictor x, y, z with 3 levels of resp: x, y, z.", fixed = TRUE)
})

test_that("print.smooth.roc of a curve built with numeric densities has no data line", {
  x <- seq(-3, 5, length.out = 50)
  s <- roc(density.controls = dnorm(x), density.cases = dnorm(x, 1))
  out <- capture.output(print(s))
  expect_false(any(grepl("^Data:", out)))
  expect_true(any(grepl("^Smoothing: density with controls: dnorm\\(x\\); and cases: dnorm\\(x, 1\\)", out)))
})

test_that("no partial matching of list elements", {
  old <- options(warnPartialMatchDollar = TRUE)
  on.exit(options(old))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_warning(smooth(r.s100b), NA)
  expect_warning(r <- roc(outcome ~ s100b, aSAH, smooth = TRUE, quiet = TRUE), NA)
  expect_warning(ci.t <- ci.thresholds(r.s100b, boot.n = 2, thresholds = 0.5), NA)
  expect_warning(capture.output(print(ci.t)), NA)
  expect_warning(plot(r.s100b), NA)
  expect_warning(plot(ci.t), NA)
  for (f in list(auc, ci, ci.auc)) {
    expect_warning(f(outcome ~ s100b, aSAH, quiet = TRUE), NA)
  }
  expect_warning(ci.se(outcome ~ s100b, aSAH, quiet = TRUE, boot.n = 2), NA)
  expect_warning(ci.sp(outcome ~ s100b, aSAH, quiet = TRUE, boot.n = 2), NA)
  expect_warning(ci.thresholds(outcome ~ s100b, aSAH, quiet = TRUE, boot.n = 2), NA)
  expect_warning(ci.coords(outcome ~ s100b, aSAH, x = 0.5, quiet = TRUE, boot.n = 2), NA)
  expect_warning(plot.roc(outcome ~ s100b, aSAH, quiet = TRUE), NA)
  expect_warning(lines.roc(outcome ~ s100b, aSAH, quiet = TRUE), NA)
  expect_warning(roc.test(outcome ~ s100b + ndka, aSAH, quiet = TRUE), NA)
})
