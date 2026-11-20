library(pROC)
data(aSAH)

context("roc.utils.percent")

test_that("roc_utils_unpercent works on full AUC", {
  expect_equal_ignore_call(pROC:::roc_utils_unpercent.roc(r.wfns.percent), r.wfns)
})

test_that("roc_utils_unpercent works on partial AUC", {
  expect_equal_ignore_call(pROC:::roc_utils_unpercent.roc(r.wfns.percent.partial1), r.wfns.partial1)
})

test_that("roc_utils_unpercent works with CI", {
  r <- roc(aSAH$outcome, aSAH$wfns, ci = TRUE)
  r.percent <- roc(aSAH$outcome, aSAH$wfns, ci = TRUE, percent = TRUE)
  expect_equal_ignore_call(pROC:::roc_utils_unpercent.roc(r.percent), r)
})

test_that("roc_utils_unpercent works without AUC", {
  r <- roc(aSAH$outcome, aSAH$wfns, auc = FALSE)
  r.percent <- roc(aSAH$outcome, aSAH$wfns, auc = FALSE, percent = TRUE)
  expect_equal_ignore_call(pROC:::roc_utils_unpercent.roc(r.percent), r)
})

test_that("roc_utils_unpercent works with ci .thresholds, .sp, .se", {
  for (of in c("thresholds", "sp", "se")) {
    set.seed(42)
    r <- roc(aSAH$outcome, aSAH$wfns, ci = TRUE, of = of, boot.n = 10)
    set.seed(42)
    r.percent <- roc(aSAH$outcome, aSAH$wfns, ci = TRUE, percent = TRUE, of = of, boot.n = 10)
    expect_equal_ignore_call(pROC:::roc_utils_unpercent.roc(r.percent), r)
  }
})

test_that("a percent curve with a ci.coords can be used in ci.auc, var, cov and power.roc.test", {
  set.seed(42)
  r <- roc(aSAH$outcome, aSAH$s100b,
    percent = TRUE, ci = TRUE, of = "coords",
    x = 90, input = "specificity", ret = "sensitivity", boot.n = 10, quiet = TRUE
  )
  expect_s3_class(r$ci, "ci.coords")
  expect_equal(as.numeric(ci.auc(r)), as.numeric(ci.auc(r.s100b.percent)))
  expect_equal(var(r), var(r.s100b.percent))
  expect_equal(var(r, method = "obuchowski"), var(r.s100b.percent, method = "obuchowski"))
  expect_equal(
    cov(r, r.ndka.percent, method = "obuchowski"),
    cov(r.s100b.percent, r.ndka.percent, method = "obuchowski")
  )
  expect_equal(power.roc.test(r)$power, power.roc.test(r.s100b.percent)$power)
  expect_equal(
    power.roc.test(r, r.ndka.percent)$power,
    power.roc.test(r.s100b.percent, r.ndka.percent)$power
  )
})

test_that("a smoothed percent curve with a ci.coords can be used in var", {
  s <- roc(aSAH$outcome, aSAH$s100b,
    percent = TRUE, smooth = TRUE, ci = TRUE, of = "coords",
    x = 90, input = "specificity", ret = "sensitivity", boot.n = 10, quiet = TRUE
  )
  expect_s3_class(s$ci, "ci.coords")
  s.noci <- smooth(r.s100b.percent)
  set.seed(42)
  v <- var(s, boot.n = 10)
  set.seed(42)
  expect_equal(v, var(s.noci, boot.n = 10))
})
