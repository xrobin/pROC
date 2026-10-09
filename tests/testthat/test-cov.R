library(pROC)
data(aSAH)

test_that("cov with delong works", {
  expect_equal(cov(r.wfns, r.ndka), -0.000532967856762438)
  expect_equal(cov(r.ndka, r.s100b), -0.000756164938056579)
  expect_equal(cov(r.s100b, r.wfns), 0.00119615567376754)
})


test_that("cov with obuchowski works", {
  expect_equal(cov(r.wfns, r.ndka, method = "obuchowski"), -3.917223e-06)
  expect_equal(cov(r.ndka, r.s100b, method = "obuchowski"), 0.0007945308)
  expect_equal(cov(r.s100b, r.wfns, method = "obuchowski"), 0.0008560803)
})


test_that("cov works with auc and mixed roc/auc", {
  expect_equal(cov(auc(r.wfns), auc(r.ndka)), -0.000532967856762438)
  expect_equal(cov(auc(r.ndka), r.s100b), -0.000756164938056579)
  expect_equal(cov(r.s100b, auc(r.wfns)), 0.00119615567376754)
})


test_that("cov with delong and percent works", {
  expect_equal(cov(r.wfns.percent, r.ndka.percent), -5.32967856762438)
  expect_equal(cov(r.ndka.percent, r.s100b.percent), -7.56164938056579)
  expect_equal(cov(r.s100b.percent, r.wfns.percent), 11.9615567376754)
})


test_that("cov with delong, percent and mixed roc/auc works", {
  expect_equal(cov(auc(r.wfns.percent), r.ndka.percent), -5.32967856762438)
  expect_equal(cov(r.ndka.percent, auc(r.s100b.percent)), -7.56164938056579)
  expect_equal(cov(auc(r.s100b.percent), auc(r.wfns.percent)), 11.9615567376754)
})


test_that("cov with obuchowski, percent and mixed roc/auc works", {
  expect_equal(cov(auc(r.wfns.percent), r.ndka.percent, method = "obuchowski"), -0.03917223)
  expect_equal(cov(r.ndka.percent, auc(r.s100b.percent), method = "obuchowski"), 7.9453082)
  expect_equal(cov(auc(r.s100b.percent), auc(r.wfns.percent), method = "obuchowski"), 8.560803)
})


test_that("cov with different auc specifications warns", {
  expect_warning(cov(r.wfns, r.ndka.percent))
  expect_warning(cov(r.wfns.percent, r.ndka))
  # Also mixing auc/roc
  expect_warning(cov(auc(r.wfns), r.ndka.percent))
  expect_warning(cov(r.wfns, auc(r.ndka.percent)))
  expect_warning(cov(r.wfns, auc(r.ndka.percent)))
})


test_that("cov with delong, percent and direction = >", {
  expect_equal(cov(r.ndka.percent, r.s100b.percent), -7.56164938056579)
})


test_that("cov with delong, percent, direction = > and mixed roc/auc", {
  r1 <- roc(aSAH$outcome, -aSAH$ndka, percent = TRUE)
  r2 <- roc(aSAH$outcome, -aSAH$s100b, percent = TRUE)
  expect_equal(cov(r1, r2), -7.56164938056579)
  expect_equal(cov(auc(r1), auc(r2)), -7.56164938056579)
  expect_equal(cov(auc(r1), r2), -7.56164938056579)
  expect_equal(cov(r1, auc(r2)), -7.56164938056579)
})


test_that("cov with bootstrap works", {
  skip_slow()
  skip_if(getRversion() < "3.6.0") # added sample.kind
  RNGkind(sample.kind = "Rejection")
  set.seed(42)
  expect_equal(cov(r.wfns, r.ndka, method = "bootstrap", boot.n = 100), -0.000648524)
  expect_equal(cov(r.ndka.percent, r.s100b.percent, method = "bootstrap", boot.n = 100), -7.17528365)
  expect_equal(cov(r.s100b.partial1, r.wfns.partial1, method = "bootstrap", boot.n = 100), 2.294465e-05)
  expect_equal(cov(r.wfns, r.ndka, method = "bootstrap", boot.n = 100, boot.stratified = FALSE), -0.0007907488)
})

test_that("bootstrap cov works with mixed roc, auc and smooth.roc objects", {
  skip_slow()
  for (roc1 in list(r.s100b, auc(r.s100b), smooth(r.s100b), r.s100b.partial2, r.s100b.partial2$auc)) {
    for (roc2 in list(r.wfns, auc(r.wfns), smooth(r.wfns), r.wfns.partial1, r.wfns.partial1$auc)) {
      n <- round(runif(1, 3, 9)) # keep boot.n small
      stratified <- sample(c(TRUE, FALSE), 1)
      suppressWarnings( # All sorts of warnings are expected
        obtained <- cov(roc1, roc2,
          method = "bootstrap",
          boot.n = n, boot.stratified = stratified
        )
      )
      expect_is(obtained, "numeric")
      expect_false(is.na(obtained))
    }
  }
})

test_that("bootstrap cov works with smooth and !reuse.auc", {
  skip_slow()
  skip_if(getRversion() < "3.6.0") # added sample.kind
  # First calculate cov by giving full curves
  roc1 <- smooth(roc(aSAH$outcome, aSAH$wfns, partial.auc = c(0.9, 1), partial.auc.focus = "sensitivity"))
  roc2 <- smooth(roc(aSAH$outcome, aSAH$s100b, partial.auc = c(0.9, 1), partial.auc.focus = "sensitivity"))

  suppressWarnings(RNGkind(sample.kind = "Rejection"))

  set.seed(42) # For reproducible CI
  expected_cov <- cov(roc1, roc2, boot.n = 100)
  expect_equal(expected_cov, -3.033016e-06)

  # Now with reuse.auc
  set.seed(42) # For reproducible CI
  obtained_cov <- cov(smooth(r.wfns), smooth(r.s100b),
    reuse.auc = FALSE,
    partial.auc = c(0.9, 1), partial.auc.focus = "sensitivity",
    boot.n = 100
  )
  expect_equal(expected_cov, obtained_cov)
})

test_that("cov keeps the partial AUC of auc objects of smoothed curves", {
  a1 <- auc(smooth(r.s100b), partial.auc = c(1, .8))
  a2 <- auc(smooth(r.ndka), partial.auc = c(1, .8))
  s1 <- smooth(roc(aSAH$outcome, aSAH$s100b, quiet = TRUE, partial.auc = c(1, .8)))
  s2 <- smooth(roc(aSAH$outcome, aSAH$ndka, quiet = TRUE, partial.auc = c(1, .8)))
  set.seed(42)
  c.auc <- cov(a1, a2, boot.n = 10)
  set.seed(42)
  c.smooth <- cov(s1, s2, boot.n = 10)
  expect_equal(c.auc, c.smooth)
})

test_that("cov errors on curves smoothed with numeric densities", {
  x <- seq(0, 1, length.out = 64)
  s1 <- smooth(r.s100b, method = "density", density.controls = dnorm(x, .2, .2), density.cases = dnorm(x, .5, .2))
  s2 <- smooth(r.ndka, method = "density", density.controls = dnorm(x, .2, .2), density.cases = dnorm(x, .4, .2))
  expect_error(suppressWarnings(cov(s1, s2, boot.n = 2)), "smoothed with numeric density.controls and density.cases")
})

test_that("cov errors on curves built from numeric densities", {
  x <- seq(-4, 6, length.out = 64)
  d1 <- roc(density.controls = dnorm(x), density.cases = dnorm(x, 1))
  d2 <- roc(density.controls = dnorm(x), density.cases = dnorm(x, 2))
  expect_error(suppressWarnings(cov(d1, d2, boot.n = 2)), "smoothed with numeric density.controls and density.cases")
})

test_that("obuchowski cov of percent curves with partial AUC is the fraction cov times 100^2", {
  expect_equal(
    cov(r.s100b.percent.partial1, r.ndka.percent.partial1, method = "obuchowski"),
    cov(r.s100b.partial1, r.ndka.partial1, method = "obuchowski") * 100^2
  )
})

test_that("cov re-pairs curves with NAs at different positions", {
  x1 <- aSAH$s100b
  x2 <- aSAH$ndka
  x1[c(3, 50)] <- NA
  x2[c(10, 90)] <- NA
  r1 <- roc(aSAH$outcome, x1, quiet = TRUE)
  r2 <- roc(aSAH$outcome, x2, quiet = TRUE)
  ok <- !is.na(x1) & !is.na(x2)
  q1 <- roc(aSAH$outcome[ok], x1[ok], quiet = TRUE)
  q2 <- roc(aSAH$outcome[ok], x2[ok], quiet = TRUE)
  expect_equal(cov(r1, r2), cov(q1, q2))
  expect_equal(cov(r1, r2, method = "obuchowski"), cov(q1, q2, method = "obuchowski"))
  set.seed(42)
  c.na <- cov(r1, r2, method = "bootstrap", boot.n = 10)
  set.seed(42)
  c.ok <- cov(q1, q2, method = "bootstrap", boot.n = 10)
  expect_equal(c.na, c.ok)
})

test_that("cov selects the bootstrap when only roc2 has a partial AUC", {
  set.seed(42)
  res.12 <- suppressWarnings(cov(r.ndka, r.s100b.partial1, boot.n = 10))
  set.seed(42)
  res.21 <- suppressWarnings(cov(r.s100b.partial1, r.ndka, boot.n = 10))
  expect_equal(res.12, res.21)
})
