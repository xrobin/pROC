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

test_that("venkatraman rejects a partial AUC on either curve", {
  expect_error(roc.test(r.s100b.partial1, r.ndka, method = "venkatraman", boot.n = 2), "Partial AUC is not supported")
  expect_error(suppressWarnings(roc.test(r.s100b, r.ndka.partial1, method = "venkatraman", boot.n = 2)), "Partial AUC is not supported")
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
