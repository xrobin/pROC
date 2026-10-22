library(pROC)
data(aSAH)

context("are.paired")

test_that("are.paired works", {
  # most basic example
  expect_true(are.paired(r.wfns, r.ndka))

  # Missing values shouldn't screw up
  aSAH.missing <- aSAH
  aSAH.missing$wfns[1:20] <- NA
  expect_true(are.paired(roc(aSAH.missing$outcome, aSAH.missing$wfns), roc(aSAH.missing$outcome, aSAH.missing$ndka)))
  # Also with different data.frames
  expect_true(are.paired(roc(aSAH.missing$outcome, aSAH.missing$wfns), r.ndka))

  # The following should fail though
  expect_false(are.paired(roc(aSAH$outcome[21:113], aSAH$wfns[21:113]), roc(aSAH$outcome, aSAH$ndka)))

  # Opposite levels should probably fail
  expect_false(are.paired(roc(aSAH$outcome, aSAH$wfns, levels = c("Good", "Poor")), roc(aSAH$outcome, aSAH$ndka, levels = c("Poor", "Good"))))
})

test_that("are.paired works with formula", {
  r.wfns.f <- roc(outcome ~ wfns, aSAH)
  r.ndka.f <- roc(outcome ~ ndka, aSAH)
  # most basic example
  expect_true(are.paired(r.wfns.f, r.ndka.f))

  # Missing values shouldn't screw up
  aSAH.missing <- aSAH
  aSAH.missing$wfns[1:20] <- NA
  expect_true(are.paired(roc(outcome ~ wfns, aSAH.missing), roc(outcome ~ ndka, aSAH.missing)))
  # Also with different data.frames
  expect_true(are.paired(roc(outcome ~ wfns, aSAH.missing), r.ndka.f))

  # The following should fail though
  expect_false(are.paired(roc(outcome ~ wfns, aSAH.missing[21:113, ]), r.ndka))

  # Opposite levels should probably fail
  expect_false(are.paired(roc(outcome ~ wfns, aSAH, levels = c("Good", "Poor")), roc(outcome ~ ndka, aSAH, levels = c("Poor", "Good"))))
})


test_that("are.paired works with auc and mixed roc", {
  expect_true(are.paired(auc(aSAH$outcome, aSAH$wfns), auc(aSAH$outcome, aSAH$ndka)))
  expect_true(are.paired(roc(aSAH$outcome, aSAH$wfns), auc(aSAH$outcome, aSAH$ndka)))
  expect_true(are.paired(auc(aSAH$outcome, aSAH$wfns), roc(aSAH$outcome, aSAH$ndka)))
})

test_that("are.paired return.paired.rocs works", {
  pair <- are.paired(r.wfns, r.ndka, return.paired.rocs = TRUE)
  expect_true(pair)
  expect_identical(attr(pair, "roc1"), r.wfns)
  expect_identical(attr(pair, "roc2"), r.ndka)
})

test_that("are.paired return.paired.rocs works with missing values", {
  aSAH.missing <- aSAH
  aSAH.missing$ndka[1:20] <- NA
  r1 <- roc(aSAH.missing$outcome, aSAH.missing$ndka)
  pair <- are.paired(r1, r.wfns, return.paired.rocs = TRUE)
  expect_true(pair)
  expect_identical(attr(pair, "roc1")$thresholds, roc(aSAH$outcome[21:113], aSAH$ndka[21:113])$thresholds)
  expect_identical(attr(pair, "roc2")$thresholds, roc(aSAH$outcome[21:113], aSAH$wfns[21:113])$thresholds)
})

test_that("are.paired return.paired.rocs doesn't return when unpaired", {
  pair <- are.paired(roc(aSAH$outcome[21:113], aSAH$wfns[21:113]), r.ndka, return.paired.rocs = TRUE)
  expect_null(attributes(pair))
})

test_that("are.paired works with smooth.roc curves", {
  expect_true(are.paired(smooth(r.wfns), smooth(r.ndka)))

  # Missing values shouldn't screw up
  aSAH.missing <- aSAH
  aSAH.missing$wfns[1:20] <- NA
  expect_true(are.paired(smooth(roc(aSAH.missing$outcome, aSAH.missing$wfns)), smooth(roc(aSAH.missing$outcome, aSAH.missing$ndka))))
  # Also with different data.frames
  expect_true(are.paired(smooth(roc(aSAH.missing$outcome, aSAH.missing$wfns)), smooth(r.ndka)))

  # The following should fail though
  expect_false(are.paired(smooth(roc(aSAH$outcome[21:113], aSAH$wfns[21:113])), smooth(roc(aSAH$outcome, aSAH$ndka))))

  # Opposite levels should probably fail
  expect_false(are.paired(smooth(roc(aSAH$outcome, aSAH$wfns, levels = c("Good", "Poor"))), smooth(roc(aSAH$outcome, aSAH$ndka, levels = c("Poor", "Good")))))
})

test_that("are.paired works with auc and mixed roc and smooth", {
  expect_true(are.paired(auc(aSAH$outcome, aSAH$wfns), smooth(roc(aSAH$outcome, aSAH$ndka))))
  expect_true(are.paired(smooth(roc(aSAH$outcome, aSAH$wfns)), auc(aSAH$outcome, aSAH$ndka)))
  expect_true(are.paired(roc(aSAH$outcome, aSAH$wfns), smooth(roc(aSAH$outcome, aSAH$ndka))))
  expect_true(are.paired(smooth(roc(aSAH$outcome, aSAH$wfns)), roc(aSAH$outcome, aSAH$ndka)))
})

test_that("are.paired return.paired.rocs returns smooth curves", {
  aSAH.missing <- aSAH
  aSAH.missing$ndka[1:20] <- NA
  r1 <- roc(aSAH.missing$outcome, aSAH.missing$ndka, smooth = TRUE)
  pair <- are.paired(r1, smooth(r.wfns), return.paired.rocs = TRUE)
  expect_true(pair)
  expect_is(attr(pair, "roc1"), "smooth.roc")
  expect_is(attr(pair, "roc2"), "smooth.roc")
})

test_that("are.paired return.paired.rocs smoothes curves with the right method", {
  skip_if_not_installed("logcondens")
  skip_slow()
  aSAH.missing <- aSAH
  aSAH.missing$ndka[1:20] <- NA
  smooth.methods <- c("binormal", "density", "fitdistr", "logcondens", "logcondens.smooth")

  for (smooth.method in smooth.methods) {
    r1 <- smooth(roc(aSAH.missing$outcome, aSAH.missing$ndka), method = smooth.method)
    pair <- are.paired(r1, smooth(r.s100b, method = smooth.method), return.paired.rocs = TRUE)
    expect_true(pair)
    expect_identical(attr(pair, "roc1")$smoothing.args$method, smooth.method)
    expect_identical(attr(pair, "roc2")$smoothing.args$method, smooth.method)
  }
})

test_that("are.paired return.paired.rocs doesn't return when unpaired and smooth", {
  pair <- are.paired(smooth(roc(aSAH$outcome[21:113], aSAH$wfns[21:113])), r.ndka, return.paired.rocs = TRUE)
  expect_null(attributes(pair))
  pair <- are.paired(roc(aSAH$outcome[21:113], aSAH$wfns[21:113]), smooth(r.ndka), return.paired.rocs = TRUE)
  expect_null(attributes(pair))
  pair <- are.paired(smooth(roc(aSAH$outcome[21:113], aSAH$wfns[21:113])), smooth(r.ndka), return.paired.rocs = TRUE)
  expect_null(attributes(pair))
})

test_that("formula and default interfaces give paired curves", {
  d <- data.frame(
    y = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
    a = c(1, 2, 3, 5, 7, 4, 6, 8, 9, 10),
    b = c(6, 1, 4, 3, 8, 2, 5, 9, 7, 10)
  )
  ra <- roc(y ~ a, d, quiet = TRUE)
  rb <- roc(d$y, d$b, quiet = TRUE)
  # The formula keeps the row names of model.response(), as in 1.19.1
  expect_identical(names(ra$response), rownames(d))
  expect_true(are.paired(ra, rb))
  expect_identical(roc.test(ra, rb)$p.value, roc.test(ra, roc(y ~ b, d, quiet = TRUE))$p.value)
  # Same data with different row names
  d2 <- d
  rownames(d2) <- letters[1:10]
  expect_true(are.paired(ra, roc(y ~ b, d2, quiet = TRUE)))
  # A named response
  expect_true(are.paired(roc(setNames(d$y, letters[1:10]), d$a, quiet = TRUE), rb))
  # A factor and a numeric response are not identical
  expect_false(are.paired(roc(factor(d$y), d$a, quiet = TRUE), rb))
})

test_that("curves that lost different observations to NA are paired on the common ones", {
  # The responses left after removing the NAs are identical, but not the
  # observations: the na.action attribute tells them apart
  y <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1)
  a <- c(NA, 2, 3, 5, 7, 4, 6, 8, 9, 10, 1, 11)
  b <- c(6, NA, 4, 3, 8, 2, 5, 9, 7, 10, 1, 12)
  ra <- roc(y, a, quiet = TRUE)
  rb <- roc(y, b, quiet = TRUE)
  expect_identical(as.vector(ra$response), as.vector(rb$response))
  pair <- are.paired(ra, rb, return.paired.rocs = TRUE)
  expect_true(pair)
  expect_length(attr(pair, "roc1")$response, 10)
  expected <- roc.test(roc(y[-(1:2)], a[-(1:2)], quiet = TRUE), roc(y[-(1:2)], b[-(1:2)], quiet = TRUE))$p.value
  expect_equal(roc.test(ra, rb)$p.value, expected)
})

test_that("are.paired doesn't use objects named oroc1/oroc2 from the user's workspace", {
  oroc1 <- oroc2 <- smooth(r.ndka)
  assign("oroc1", oroc1, envir = globalenv())
  assign("oroc2", oroc2, envir = globalenv())
  on.exit(rm("oroc1", "oroc2", envir = globalenv()))
  r1 <- roc(aSAH$outcome, c(NA, aSAH$s100b[-1]), quiet = TRUE)
  pair <- are.paired(r1, r.wfns, return.paired.rocs = TRUE)
  expect_true(pair)
  expect_s3_class(attr(pair, "roc1"), "roc")
  expect_false(inherits(attr(pair, "roc1"), "smooth.roc"))
  expect_false(inherits(attr(pair, "roc2"), "smooth.roc"))
  expect_equal(
    roc.test(r1, r.wfns)$p.value,
    roc.test(roc(aSAH$outcome[-1], aSAH$s100b[-1], quiet = TRUE), roc(aSAH$outcome[-1], aSAH$wfns[-1], quiet = TRUE))$p.value
  )
})

test_that("are.paired return.paired.rocs returns smoothed curves with their AUC", {
  # Without NA: the smoothed curves are returned as is
  s1 <- smooth(r.s100b)
  s2 <- smooth(r.wfns)
  pair <- are.paired(s1, s2, return.paired.rocs = TRUE)
  expect_true(pair)
  expect_identical(attr(pair, "roc1"), s1)
  expect_identical(attr(pair, "roc2"), s2)
  pair <- are.paired(s1, s2, return.paired.rocs = TRUE, reuse.smooth = FALSE)
  expect_false(inherits(attr(pair, "roc1"), "smooth.roc"))

  # With NA: the AUC specification of the smoothed curve is reused
  s100b.na <- c(NA, aSAH$s100b[-1])
  sa <- roc(aSAH$outcome, s100b.na, smooth = TRUE, partial.auc = c(1, 0.8), quiet = TRUE)
  sb <- roc(aSAH$outcome, aSAH$wfns, smooth = TRUE, partial.auc = c(1, 0.8), quiet = TRUE)
  pair <- are.paired(sa, sb, return.paired.rocs = TRUE)
  expect_true(pair)
  expect_s3_class(attr(pair, "roc1"), "smooth.roc")
  expect_equal(attr(attr(pair, "roc1")$auc, "partial.auc"), c(1, 0.8))
  expect_equal(as.numeric(attr(pair, "roc1")$auc), as.numeric(sa$auc))
  expect_equal(
    as.numeric(attr(pair, "roc2")$auc),
    as.numeric(roc(aSAH$outcome[-1], aSAH$wfns[-1], smooth = TRUE, partial.auc = c(1, 0.8), quiet = TRUE)$auc)
  )
})

test_that("explicit numeric or factor levels pair with the default character levels", {
  y <- as.numeric(aSAH$outcome == "Poor")
  r.default <- roc(y, aSAH$ndka, quiet = TRUE)
  r.numeric <- roc(y, aSAH$s100b, levels = c(0, 1), quiet = TRUE)
  expect_true(are.paired(r.numeric, r.default))
  expect_true(are.paired(r.default, r.numeric))
  expect_equal(roc.test(r.numeric, r.default)$method, "DeLong's test for two correlated ROC curves")
  expect_equal(cov(r.numeric, r.default), cov(r.s100b, r.ndka))
  expect_equal(power.roc.test(r.numeric, r.default)$power, power.roc.test(r.s100b, r.ndka)$power)
  r.factor <- roc(aSAH$outcome, aSAH$s100b, levels = unique(aSAH$outcome), quiet = TRUE)
  expect_true(are.paired(r.factor, r.ndka))
  # different levels are still unpaired
  expect_false(are.paired(roc(y, aSAH$s100b, levels = c(1, 0), quiet = TRUE), r.default))
})

test_that("re-paired curves print the data of the input curves", {
  y <- c(0, 0, 0, 0, 1, 1, 1, 1)
  a <- c(1, NA, 3, 4, 2, 5, 6, 7)
  b <- c(2, 1, 4, NA, 3, 6, 5, 8)
  ra <- roc(y, a, quiet = TRUE)
  rb <- roc(y, b, quiet = TRUE)
  printed <- function(x) paste(capture.output(print(x)), collapse = "\n")
  pair <- are.paired(ra, rb, return.paired.rocs = TRUE)
  expect_false(grepl("idx.exclude", printed(attr(pair, "roc1"))))
  expect_match(printed(attr(pair, "roc1")), "Data: a in 2 controls (y 0) < 4 cases (y 1).", fixed = TRUE)
  expect_match(printed(attr(pair, "roc2")), "Data: b in 2 controls (y 0) < 4 cases (y 1).", fixed = TRUE)
  expect_match(printed(roc.test(ra, rb)$roc1), "roc.default(response = y, predictor = a, quiet = TRUE)", fixed = TRUE)
  # formula and smoothed curves
  d <- data.frame(y = y, a = a, b = b)
  rf <- roc(y ~ a, d, quiet = TRUE)
  expect_match(printed(attr(are.paired(rf, rb, return.paired.rocs = TRUE), "roc1")), "Data: a in 2 controls (y 0) < 4 cases (y 1).", fixed = TRUE)
  s <- smooth(ra, method = "density")
  out <- printed(attr(are.paired(s, rb, return.paired.rocs = TRUE), "roc1"))
  expect_false(grepl("paired|idx.exclude", out))
  expect_match(out, "Data: a in 2 controls (y 0) < 4 cases (y 1).", fixed = TRUE)
})


test_that("re-paired curves keep the best.method, best.weights and best.policy of a stored ci.coords", {
  set.seed(42)
  a <- aSAH$s100b
  a[1] <- NA
  b <- aSAH$ndka
  b[2] <- NA
  ra <- roc(aSAH$outcome, a, quiet = TRUE)
  rb <- roc(aSAH$outcome, b, quiet = TRUE)
  ra$ci <- ci.coords(ra, "best",
    input = "specificity", ret = c("specificity", "sensitivity"),
    best.method = "closest.topleft", best.weights = c(10, 0.1), best.policy = "random",
    boot.n = 2, progress = "none"
  )
  pair <- are.paired(ra, rb, return.paired.rocs = TRUE, reuse.ci = TRUE)
  ci <- attr(pair, "roc1")$ci
  expect_s3_class(ci, "ci.coords")
  expect_identical(attr(ci, "best.method"), "closest.topleft")
  expect_identical(attr(ci, "best.weights"), c(10, 0.1))
  expect_identical(attr(ci, "best.policy"), "random")
})
