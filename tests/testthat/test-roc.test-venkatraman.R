library(pROC)
data(aSAH)

test_that("paired venkatraman works as expected", {
  skip_slow()
  ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12)
  expect_venkatraman_htest(ht)
  expect_equal(ht$alternative, "two.sided")
  expect_equal(ht$method, "Venkatraman's test for two paired ROC curves")
  expect_equal(unname(ht$parameter), 12)
  # Test output
  ht$statistic <- c(E = 42)
  ht$p.value <- 0
  expect_known_output(print(ht), "print_output/roc.test-venkatraman.paired")
})

test_that("unpaired venkatraman works as expected", {
  skip_slow()
  expect_warning(ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, paired = FALSE), "paired")
  expect_venkatraman_htest(ht)
  expect_equal(ht$alternative, "two.sided")
  expect_equal(ht$method, "Venkatraman's test for two unpaired ROC curves")
  expect_equal(unname(ht$parameter), 12)
  # Test output
  ht$statistic <- c(E = 41)
  ht$p.value <- 0.548347196932
  expect_known_output(print(ht), "print_output/roc.test-venkatraman.unpaired")
})

test_that("non stratified venkatraman works as expected", {
  skip_slow()
  ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, boot.stratified = FALSE)
  expect_venkatraman_htest(ht)
  expect_equal(ht$alternative, "two.sided")
  expect_equal(ht$method, "Venkatraman's test for two paired ROC curves")
  expect_equal(unname(ht$parameter), 12)
  # Test output
  ht$statistic <- c(E = 43)
  ht$p.value <- 0.05
  expect_known_output(print(ht), "print_output/roc.test-venkatraman.unstratified")
})

test_that("non stratified, unpaired venkatraman works as expected", {
  skip_slow()
  expect_warning(ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, boot.stratified = FALSE, paired = FALSE), "paired")
  expect_venkatraman_htest(ht)
  expect_equal(ht$alternative, "two.sided")
  expect_equal(ht$method, "Venkatraman's test for two unpaired ROC curves")
  expect_equal(unname(ht$parameter), 12)
  # Test output
  ht$statistic <- c(E = 43)
  ht$p.value <- 0.05
  expect_known_output(print(ht), "print_output/roc.test-venkatraman.unpaired.unstratified")
})

test_that("unpaired venkatraman statistic is symmetric with curves of different sizes", {
  set.seed(42)
  r1 <- roc(c(0, 0, 0, 1, 1, 1), c(1, 2, 4, 3, 5, 6), quiet = TRUE)
  r2 <- roc(c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), c(1, 3, 2, 4:12), quiet = TRUE)
  t12 <- roc.test(r1, r2, method = "venkatraman", paired = FALSE, boot.n = 2)
  t21 <- roc.test(r2, r1, method = "venkatraman", paired = FALSE, boot.n = 2)
  expect_equal(t12$statistic, t21$statistic)
})

test_that("venkatraman breaks ties at random by default", {
  set.seed(42)
  t.default <- roc.test(r.wfns, r.s100b, method = "venkatraman", boot.n = 10)
  set.seed(42)
  t.random <- roc.test(r.wfns, r.s100b, method = "venkatraman", boot.n = 10, ties.method = "random")
  expect_equal(t.default$statistic, t.random$statistic)
  expect_equal(t.default$p.value, t.random$p.value)
})

test_that("paired venkatraman of identical tied predictors gives E = 0 and p = 1", {
  set.seed(42)
  for (r in list(r.wfns, r.s100b)) {
    t <- roc.test(r, r, method = "venkatraman", boot.n = 20)
    expect_equal(unname(t$statistic), 0)
    expect_equal(t$p.value, 1)
  }
  # Same predictor values in two separately built curves
  r2 <- roc(aSAH$outcome, as.numeric(aSAH$wfns), levels = c("Good", "Poor"), direction = "<", quiet = TRUE)
  r3 <- roc(aSAH$outcome, as.numeric(aSAH$wfns), levels = c("Good", "Poor"), direction = "<", quiet = TRUE)
  expect_equal(roc.test(r2, r3, method = "venkatraman", boot.n = 20)$p.value, 1)
})

test_that("paired venkatraman draws no random numbers without ties", {
  # Untied data give the result of 1.19.1 (ties.method = "first") for a seed
  set.seed(42)
  d <- rep(0:1, each = 30)
  x <- rnorm(60) + d
  y <- rnorm(60) + 0.8 * d
  r1 <- roc(d, x, quiet = TRUE, direction = "<")
  r2 <- roc(d, y, quiet = TRUE, direction = "<")
  set.seed(1)
  t.random <- roc.test(r1, r2, method = "venkatraman", boot.n = 100)
  set.seed(1)
  t.first <- roc.test(r1, r2, method = "venkatraman", boot.n = 100, ties.method = "first")
  expect_identical(t.random$statistic, t.first$statistic)
  expect_identical(t.random$p.value, t.first$p.value)
})

test_that("paired venkatraman reproduces Example 4.1 of Venkatraman and Begg (1996)", {
  # Melanoma data of Table 4 (72 patients), one row per patient
  ex41 <- data.frame(
    x = c(-1.571, -3.412, 2.398, -1.622, -4.363, -4.788, -2.254, -2.763, -4.576, -2.716, -4.764, -0.055, -1.677,
          1.183, -2.293, 0.183, -3.503, -2.245, -4.06, -0.091, -2.089, -4.952, -0.079, -2.389, -4.023, -1.43, -3.221,
          -4.376, 0.406, 0.786, 0.314, -2.871, -5.881, 3.032, -1.112, 0.298, -1.738, 0.359, -3.221, 1.407, -3.471,
          -1.057, -2.408, -3.06, 1.631, -2.102, -0.052, 1.751, -1.623, -1.714, 1.939, -1.597, -4.717, 1.491, -2.457,
          1.007, -3.776, -1.924, -3.835, -0.441, -1.174, -5.164, 1.195, -1.75, -0.973, -1.726, -0.961, -4.412, 2.762,
          -1.193, -0.421, 1.939),
    y = c(-1.941, -4.933, 0.916, -2.857, -6.524, -6.568, -1.443, 0.032, -4.131, -5.297, -3.476, 2.146, 1.279,
          -0.805, -0.939, -0.314, -1.121, -3.641, -4.399, 0.483, -2.73, -3.808, -0.656, -1.118, -3.493, 1.548, -5.783,
          -3.799, -0.435, 2.111, 0.887, -3.244, -7.103, 3.258, -0.395, -0.323, -3.24, 0.258, -5.418, -0.476, -3.502,
          -2.609, -3.831, -4.621, 1.695, -3.545, -0.119, 1.262, -1.932, -2.22, 2.407, -1.728, -4.315, 0.247, -2.808,
          0.593, -4.108, -0.676, -4.782, -1.847, -0.744, -6.344, -0.77, -4.206, 0.241, -2.085, -1.472, -4.986, 1.563,
          -0.998, 1.873, 0.399),
    d = c(1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1,
          0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
          1, 1))
  r1 <- roc(ex41$d, ex41$x, quiet = TRUE, direction = "<")
  r2 <- roc(ex41$d, ex41$y, quiet = TRUE, direction = "<")
  # The paper reports E = 64 and p = 0.88 with 1000 permutations
  for (ties.method in c("random", "first")) {
    set.seed(42)
    t <- roc.test(r1, r2, method = "venkatraman", boot.n = 1000, ties.method = ties.method)
    expect_equal(unname(t$statistic), 64)
  }
  # The permutation p-value is random: over 40 seeds it ranges from 0.88 to
  # 0.92 here, and its Monte Carlo standard error at 1000 permutations is
  # about 0.01. Check that it is within 0.05 of the paper's value.
  set.seed(1)
  t <- roc.test(r1, r2, method = "venkatraman", boot.n = 1000)
  expect_lt(abs(t$p.value - 0.88), 0.05)
})

test_that("venkatraman rejects a partial AUC on either curve", {
  set.seed(42)
  expect_warning(
    expect_error(roc.test(r.s100b.partial1, r.ndka, method = "venkatraman", boot.n = 2), "Partial AUC is not supported"),
    "Different AUC specifications"
  )
  expect_warning(
    expect_error(roc.test(r.s100b, r.ndka.partial1, method = "venkatraman", boot.n = 2), "Partial AUC is not supported"),
    "Different AUC specifications"
  )
})

test_that("unpaired venkatraman statistic integrates the error curves over [0, 1]", {
  # The error curves e(p) = mp * Fcases + (1 - mp) * (1 - Fcontrols) against the
  # mixture quantile p start at (0, 1 - mp) and end at (1, mp); E is the area
  # between them.
  stat <- function(controls1, cases1, controls2, cases2) {
    r1 <- roc(controls = controls1, cases = cases1, direction = "<", quiet = TRUE)
    r2 <- roc(controls = controls2, cases = cases2, direction = "<", quiet = TRUE)
    unname(suppressWarnings(roc.test(r1, r2, method = "venkatraman", paired = FALSE, boot.n = 1))$statistic)
  }
  # e1: (0, .5) (.25, .25) (.5, 0) (.75, .25) (1, .5)
  # e2: (0, .5) (.25, .75) (.5, .5) (.75, .25) (1, .5)
  # The area between p = 0 and .25 (1/16) used to be left out
  expect_equal(stat(c(1, 2), c(3, 4), c(2, 3), c(1, 4)), 1 / 4)
  # e1: (0, .5) (1/6, 1/3) (2/3, 5/6) (5/6, 2/3) (1, .5)
  # e2: (0, .5) (.5, 0) (2/3, 1/6) (5/6, 1/3) (1, .5)
  expect_equal(stat(c(1, 3, 4), 2, 1, c(2, 3, 4)), 1 / 3)
  # e1: (0, .5) (.25, .25) (.75, .75) (1, .5)
  # e2: (0, .5) (.25, .75) (.75, .25) (1, .5)
  # The curves cross at p = .5: 4 triangles of area 1/16
  expect_equal(stat(c(1, 3), 2, 2, c(1, 3)), 1 / 4)
  # Identical curves
  expect_equal(stat(c(1, 3), 2, c(1, 3), 2), 0)
})

test_that("unpaired venkatraman permutes curves of different designs on a common scale", {
  # Both curves separate cases from controls perfectly. On the mixture scale
  # of the statistic, every control of either curve lies below every case,
  # so every permutation must give two perfect curves again (E = 0). Permuting
  # raw ranks of curves of different sizes or case mixes did not.
  r1 <- roc(controls = 1:8, cases = 9:10, direction = "<", quiet = TRUE)
  r2 <- roc(controls = 1:20, cases = 21:100, direction = "<", quiet = TRUE)
  set.seed(42)
  res <- pROC:::venkatraman.unpaired.test(r1, r2, boot.n = 200)
  expect_lt(res[[1]], 1e-12)
  expect_true(all(res[[2]] < 1e-12))
})
