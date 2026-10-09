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
  # Save the user's objects of these names and restore them on exit
  names <- c("oroc1", "oroc2")
  existing <- names[vapply(names, exists, logical(1), envir = globalenv(), inherits = FALSE)]
  saved <- mget(existing, envir = globalenv())
  on.exit({
    rm(list = names, envir = globalenv())
    list2env(saved, envir = globalenv())
  })
  assign("oroc1", oroc1, envir = globalenv())
  assign("oroc2", oroc2, envir = globalenv())
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
