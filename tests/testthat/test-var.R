library(pROC)
data(aSAH)


test_that("var with delong works", {
  expect_equal(var(r.wfns), 0.00146991470882363)
  expect_equal(var(r.ndka), 0.0031908105493913)
  expect_equal(var(r.s100b), 0.00266868245717244)
})


test_that("var works with auc", {
  expect_equal(var(auc(r.wfns)), 0.00146991470882363)
  expect_equal(var(auc(r.ndka)), 0.0031908105493913)
  expect_equal(var(auc(r.s100b)), 0.00266868245717244)
})


test_that("var with delong and percent works", {
  expect_equal(var(r.wfns.percent), 14.6991470882363)
  expect_equal(var(r.ndka.percent), 31.908105493913)
  expect_equal(var(r.s100b.percent), 26.6868245717244)
})


test_that("var works with auc and percent", {
  expect_equal(var(auc(r.wfns.percent)), 14.6991470882363)
  expect_equal(var(auc(r.ndka.percent)), 31.908105493913)
  expect_equal(var(auc(r.s100b.percent)), 26.6868245717244)
})


test_that("var with delong and percent works", {
  expect_equal(var(roc(aSAH$outcome, -aSAH$ndka, percent = TRUE)), 31.908105493913)
  expect_equal(var(roc(aSAH$outcome, -aSAH$s100b, percent = TRUE)), 26.6868245717244)
})

# Only test whether var runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.
test_that("bootstrap var runs with roc, auc and smooth.roc objects", {
  skip_slow()
  for (roc1 in list(r.s100b, auc(r.s100b), smooth(r.s100b), r.s100b.partial1, r.s100b.partial1$auc)) {
    n <- round(runif(1, 3, 9)) # keep boot.n small
    for (stratified in c(TRUE, FALSE)) {
      stratified <- sample(c(TRUE, FALSE), 1)
      obtained <- var(roc1,
        method = "bootstrap",
        boot.n = n, boot.stratified = stratified
      )
      expect_is(obtained, "numeric")
      expect_false(is.na(obtained))
    }
  }
})

test_that("var errors on curves smoothed with numeric densities", {
  x <- seq(0, 1, length.out = 64)
  s1 <- smooth(r.s100b, method = "density", density.controls = dnorm(x, .2, .2), density.cases = dnorm(x, .5, .2))
  s2 <- smooth(r.ndka, method = "density", density.controls = dnorm(x, .2, .2), density.cases = dnorm(x, .4, .2))
  expect_error(suppressWarnings(var(s1, boot.n = 2)), "smoothed with numeric density.controls and density.cases")
})

test_that("var errors on curves built from numeric densities", {
  x <- seq(-4, 6, length.out = 64)
  d1 <- roc(density.controls = dnorm(x), density.cases = dnorm(x, 1))
  d2 <- roc(density.controls = dnorm(x), density.cases = dnorm(x, 2))
  expect_error(suppressWarnings(var(d1, boot.n = 2)), "smoothed with numeric density.controls and density.cases")
})

test_that("var with unsupported obuchowski cases errors without claiming to use bootstrap", {
  expect_error(var(smooth(r.s100b), method = "obuchowski"), "Use method=\"bootstrap\" instead")
  expect_error(var(r.s100b.partial2, method = "obuchowski"), "Use method=\"bootstrap\" instead")
})

test_that("var works on smoothed percent curves", {
  s <- smooth(r.s100b)
  s.percent <- smooth(r.s100b.percent)
  seed <- sample.int(1e6, 1)
  set.seed(seed)
  v <- var(s, boot.n = 20, progress = "none")
  set.seed(seed)
  v.percent <- var(s.percent, boot.n = 20, progress = "none")
  expect_equal(v.percent, v * 100^2)
  set.seed(seed)
  expect_equal(var(auc(s.percent), boot.n = 20, progress = "none"), v * 100^2)
})
