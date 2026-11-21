context("case weights")

# Whole-number weights, with zeros: the curve of the replicated observations
k <- rep_len(c(2, 0, 1, 3, 1), nrow(aSAH))
# Fractional weights
w <- rep_len(c(0.3, 1.7, 0.9, 2.2, 0.1, 1), nrow(aSAH))

test_that("unit weights give the unweighted curve", {
  for (marker in c("s100b", "wfns")) {
    r <- roc(aSAH$outcome, aSAH[[marker]], quiet = TRUE)
    r.w <- roc(aSAH$outcome, aSAH[[marker]], weights = rep(1, nrow(aSAH)), quiet = TRUE)
    expect_identical(r.w$sensitivities, r$sensitivities)
    expect_identical(r.w$specificities, r$specificities)
    expect_identical(r.w$thresholds, r$thresholds)
    expect_identical(as.numeric(r.w$auc), as.numeric(r$auc))
  }
})

test_that("whole-number weights give the curve of the replicated observations", {
  for (marker in c("s100b", "wfns")) {
    for (percent in c(FALSE, TRUE)) {
      r.w <- roc(aSAH$outcome, aSAH[[marker]], weights = k, percent = percent, quiet = TRUE)
      r.rep <- roc(rep(aSAH$outcome, k), rep(aSAH[[marker]], k), percent = percent, quiet = TRUE)
      expect_identical(r.w$sensitivities, r.rep$sensitivities)
      expect_identical(r.w$specificities, r.rep$specificities)
      expect_identical(r.w$thresholds, r.rep$thresholds)
      expect_identical(r.w$direction, r.rep$direction)
      expect_identical(as.numeric(r.w$auc), as.numeric(r.rep$auc))
      for (x in c("all", "local maximas", "best")) {
        expect_identical(coords(r.w, x, ret = "all"), coords(r.rep, x, ret = "all"))
      }
    }
  }
  # the observations with a weight of 0 are not in the curve
  r.w <- roc(aSAH$outcome, aSAH$s100b, weights = k, quiet = TRUE)
  expect_identical(r.w$predictor, aSAH$s100b[k > 0])
  expect_identical(r.w$weights, k[k > 0])
  expect_identical(r.w$original.weights, k)
})

test_that("the weighted AUC is the weighted proportion of concordant pairs", {
  r.w <- roc(aSAH$outcome, aSAH$s100b, weights = w, quiet = TRUE)
  cases <- which(aSAH$outcome == "Poor")
  controls <- which(aSAH$outcome == "Good")
  concordant <- outer(aSAH$s100b[cases], aSAH$s100b[controls], ">") +
    outer(aSAH$s100b[cases], aSAH$s100b[controls], "==") / 2
  expected <- sum(concordant * outer(w[cases], w[controls])) / (sum(w[cases]) * sum(w[controls]))
  expect_equal(as.numeric(r.w$auc), expected)
  # weighted counts
  co <- coords(r.w, 0.205, ret = c("tp", "fn", "tn", "fp"))
  expect_equal(co$tp, sum(w[cases][aSAH$s100b[cases] > 0.205]))
  expect_equal(co$tp + co$fn, sum(w[cases]))
  expect_equal(co$tn, sum(w[controls][aSAH$s100b[controls] < 0.205]))
  expect_equal(co$tn + co$fp, sum(w[controls]))
})

test_that("the weighted curve starts and ends exactly", {
  r.w <- roc(aSAH$outcome, aSAH$s100b, weights = w, percent = TRUE, quiet = TRUE)
  expect_identical(range(r.w$sensitivities), c(0, 100))
  expect_identical(range(r.w$specificities), c(0, 100))
})

test_that("direction = 'auto' compares the weighted medians", {
  # Unweighted medians: controls 2 < cases 3. Weighted: controls 2 > cases 1
  response <- c(0, 0, 0, 1, 1, 1)
  predictor <- c(1, 2, 3, 1, 3, 4)
  expect_identical(roc(response, predictor, quiet = TRUE)$direction, "<")
  expect_identical(roc(response, predictor, weights = c(1, 1, 1, 10, 1, 1), quiet = TRUE)$direction, ">")
})

test_that("weights work with a formula, missing values and subset", {
  aSAH.w <- aSAH
  aSAH.w$w <- w
  aSAH.w$w[3] <- NA
  expected <- roc(aSAH.w$outcome[-3], aSAH.w$s100b[-3], weights = w[-3], quiet = TRUE)
  r.f <- roc(outcome ~ s100b, aSAH.w, weights = w, quiet = TRUE)
  expect_identical(r.f$sensitivities, expected$sensitivities)
  expect_identical(as.numeric(r.f$auc), as.numeric(expected$auc))
  expect_identical(roc(aSAH.w$outcome, aSAH.w$s100b, weights = aSAH.w$w, na.rm = FALSE, quiet = TRUE), NA)

  female <- aSAH$gender == "Female"
  r.sub <- roc(outcome ~ s100b + ndka, aSAH, weights = w, subset = gender == "Female", quiet = TRUE)
  expect_identical(
    as.numeric(r.sub$ndka$auc),
    as.numeric(roc(aSAH$outcome[female], aSAH$ndka[female], weights = w[female], quiet = TRUE)$auc)
  )
  # the other formula methods refuse them
  expect_error(auc(outcome ~ s100b, aSAH, weights = w), "weights are not supported")
})

test_that("invalid weights are refused", {
  expect_error(roc(aSAH$outcome, aSAH$s100b, weights = w[-1], quiet = TRUE), "same length")
  expect_error(roc(aSAH$outcome, aSAH$s100b, weights = -w, quiet = TRUE), "non-negative")
  expect_error(roc(aSAH$outcome, aSAH$s100b, weights = as.character(w), quiet = TRUE), "numeric")
  expect_error(roc(controls = 1:3, cases = 2:4, weights = rep(1, 6), quiet = TRUE), "response")
  expect_error(roc(aSAH$outcome, aSAH$s100b, weights = as.numeric(aSAH$outcome == "Good"), quiet = TRUE), "No case observation")
})

test_that("print shows the case weights", {
  expect_output(print(roc(aSAH$outcome, aSAH$s100b, weights = k, quiet = TRUE)), "Case weights: 91 controls, 66 cases.", fixed = TRUE)
})

test_that("the bootstrap works with case weights", {
  r.w <- roc(aSAH$outcome, aSAH$s100b, weights = w, quiet = TRUE)
  r.w2 <- roc(aSAH$outcome, aSAH$ndka, weights = w, quiet = TRUE)
  for (stratified in c(TRUE, FALSE)) {
    ci <- ci.auc(r.w, boot.n = 10, boot.stratified = stratified, progress = FALSE)
    expect_identical(attr(ci, "method"), "bootstrap")
    expect_false(is.unsorted(ci))
    expect_true(is.finite(var(r.w, boot.n = 10, boot.stratified = stratified, progress = FALSE)))
    expect_true(is.finite(cov(r.w, r.w2, boot.n = 10, boot.stratified = stratified, progress = FALSE)))
  }
  expect_identical(dim(ci.thresholds(r.w, thresholds = c(0.2, 0.5), boot.n = 10, progress = FALSE)$sensitivity), c(2L, 3L))
  expect_identical(dim(ci.se(r.w, specificities = c(0.5, 0.9), boot.n = 10, progress = FALSE)), c(2L, 3L))
  expect_identical(dim(ci.sp(r.w, sensitivities = c(0.5, 0.9), boot.n = 10, progress = FALSE)), c(2L, 3L))
  expect_identical(dim(ci.coords(r.w, x = 0.2, ret = "tp", boot.n = 10, progress = FALSE)$tp), c(1L, 3L))

  # two curves
  expect_true(are.paired(r.w, r.w2))
  expect_message(test <- roc.test(r.w, r.w2, boot.n = 10, progress = FALSE), "paired: TRUE")
  expect_identical(unname(test$estimate), c(as.numeric(r.w$auc), as.numeric(r.w2$auc)))
  expect_match(test$method, "Bootstrap")
  # the same resamples for both curves
  expect_identical(unname(roc.test(r.w, r.w, boot.n = 10, progress = FALSE, paired = TRUE)$statistic), 0)
  expect_match(roc.test(r.w, r.w2, paired = FALSE, boot.n = 10, progress = FALSE)$method, "Bootstrap")
  expect_match(roc.test(r.w, r.w2, method = "specificity", specificity = 0.9, paired = TRUE, boot.n = 10, progress = FALSE)$method, "Specificity")
})

test_that("curves with different missing or null weights are paired on their common observations", {
  k2 <- k
  k2[c(1, 3)] <- c(NA, 0)
  r.w <- roc(aSAH$outcome, aSAH$s100b, weights = k, quiet = TRUE)
  r.w2 <- roc(aSAH$outcome, aSAH$ndka, weights = k2, quiet = TRUE)
  paired <- are.paired(r.w, r.w2, return.paired.rocs = TRUE)
  expect_true(paired)
  keep <- k > 0 & !is.na(k2) & k2 > 0
  expect_identical(attr(paired, "roc1")$predictor, aSAH$s100b[keep])
  expect_identical(attr(paired, "roc2")$predictor, aSAH$ndka[keep])
  expect_identical(attr(paired, "roc1")$weights, k[keep])
  expect_identical(attr(paired, "roc2")$weights, k2[keep])
})

test_that("what has no weighted version is refused", {
  r.w <- roc(aSAH$outcome, aSAH$s100b, weights = w, quiet = TRUE)
  refused <- "not supported for ROC curves with case weights"
  expect_error(ci.auc(r.w, method = "delong"), refused)
  expect_error(var(r.w, method = "delong"), refused)
  expect_error(var(r.w, method = "obuchowski"), refused)
  expect_error(smooth(r.w), refused)
  expect_error(roc(aSAH$outcome, aSAH$s100b, weights = w, smooth = TRUE, quiet = TRUE), refused)
  expect_error(roc.test(r.w, r.ndka, method = "delong", paired = TRUE), refused)
  expect_error(roc.test(r.ndka, r.w, method = "venkatraman", paired = TRUE), refused)
  expect_error(cov(r.w, r.ndka, method = "delong", paired = TRUE), refused)
  expect_error(cov(r.ndka, r.w, method = "obuchowski", paired = TRUE), refused)
  expect_error(power.roc.test(r.w), refused)
  expect_error(power.roc.test(r.ndka, r.w), refused)
  expect_error(multiclass.roc(aSAH$gos6, aSAH$s100b, weights = w), "not implemented")
})
