# How smoothed curves store their smoothing arguments (smoothing.args), how
# roc.test and cov compare them, and how the smoothing is replayed.
# Expected behaviour: 1.19.1's, except where noted (B7, B8, B9).

msg.smoothing <- "Different smoothing parameters"

test_that("B1: roc(smooth = TRUE) stores all its ... in smoothing.args, as in 1.19.1", {
  pdf(NULL)
  on.exit(dev.off())
  s <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, plot = TRUE, col = "red", main = "s100b", quiet = TRUE)
  expect_identical(names(s$smoothing.args), c("col", "main", "method", "n", "bw"))
})

test_that("B2: curves smoothed the same way are compared without warning", {
  set.seed(42)
  expect_no_warning(roc.test(smooth(r.s100b), smooth(r.ndka), boot.n = 5), message = msg.smoothing)
})

test_that("B3: the plot arguments of roc(smooth = TRUE) count as smoothing arguments, as in 1.19.1", {
  pdf(NULL)
  on.exit(dev.off())
  s1 <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, plot = TRUE, col = "red", quiet = TRUE)
  set.seed(42)
  expected <- roc.test(smooth(r.s100b), smooth(r.ndka), boot.n = 5)$p.value
  set.seed(42)
  expect_warning(t <- roc.test(s1, smooth(r.ndka), boot.n = 5), msg.smoothing)
  expect_equal(t$p.value, expected)
})

test_that("B4: an explicit default counts as a difference, as in 1.19.1", {
  set.seed(42)
  expected <- roc.test(smooth(r.s100b, method = "density"), smooth(r.ndka, method = "density"), boot.n = 5)$p.value
  set.seed(42)
  expect_warning(t <- roc.test(smooth(r.s100b, method = "density", adjust = 1), smooth(r.ndka, method = "density"), boot.n = 5), msg.smoothing)
  expect_equal(t$p.value, expected)
})

test_that("B5: n = 512L counts as a difference from the default 512, as in 1.19.1", {
  set.seed(42)
  expect_warning(cov(smooth(r.s100b, n = 512L), smooth(r.ndka), boot.n = 5), msg.smoothing)
})

test_that("B6: an abbreviated kernel counts as a difference, as in 1.19.1", {
  set.seed(42)
  expected <- roc.test(smooth(r.s100b, method = "density", kernel = "epanechnikov"), smooth(r.ndka, method = "density", kernel = "epanechnikov"), boot.n = 5)$p.value
  set.seed(42)
  expect_warning(t <- roc.test(smooth(r.s100b, method = "density", kernel = "e"), smooth(r.ndka, method = "density", kernel = "epanechnikov"), boot.n = 5), msg.smoothing)
  expect_equal(t$p.value, expected)
})

test_that("B7: ci of a curve with an expression in smoothing.args (1.19.1: NA)", {
  pdf(NULL)
  on.exit(dev.off())
  s <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, plot = TRUE, main = quote(italic(s100b)), quiet = TRUE)
  expect_identical(s$smoothing.args$main, quote(italic(s100b)))
  set.seed(42)
  expected <- ci(smooth(r.s100b), boot.n = 5)
  set.seed(42)
  expect_equal(as.numeric(ci(s, boot.n = 5)), as.numeric(expected))
})

test_that("B8: roc.test of a curve with an expression in smoothing.args (1.19.1: error)", {
  pdf(NULL)
  on.exit(dev.off())
  s <- roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, plot = TRUE, main = quote(italic(s100b)), quiet = TRUE)
  set.seed(42)
  expected <- roc.test(smooth(r.s100b), smooth(r.ndka), boot.n = 5)$p.value
  set.seed(42)
  expect_warning(t <- roc.test(s, smooth(r.ndka), boot.n = 5), msg.smoothing)
  expect_equal(t$p.value, expected)
  # The returned curve keeps the stored expression
  expect_identical(t$roc1$smoothing.args$main, quote(italic(s100b)))
})

test_that("B9: roc.test of the partial AUC of smoothed curves (1.19.1: full AUC)", {
  a1 <- auc(smooth(r.s100b), partial.auc = c(1, .8))
  a2 <- auc(smooth(r.ndka), partial.auc = c(1, .8))
  set.seed(42)
  t <- roc.test(a1, a2, boot.n = 5)
  expect_equal(unname(t$estimate), c(as.numeric(a1), as.numeric(a2)))
})
