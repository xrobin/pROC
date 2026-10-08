library(pROC)
data(aSAH)

# 'progress' and 'parallel' were deprecated in 1.19 when the bootstrap became
# sequential. Both are still accepted so old scripts keep running, and both are
# ignored. The documentation promises a warning in each case.

# Keep boot.n tiny: these tests check the warnings, not the bootstrap.
B <- 3

test_that("'parallel' warns and is ignored", {
  expect_warning(ci.auc(r.wfns, method = "bootstrap", boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
  expect_warning(ci.se(r.wfns, boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
  expect_warning(ci.sp(r.wfns, boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
  expect_warning(ci.thresholds(r.wfns, boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
  expect_warning(var(r.wfns, method = "bootstrap", boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
  expect_warning(cov(r.wfns, r.ndka, method = "bootstrap", boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
  expect_warning(roc.test(r.wfns, r.ndka, method = "bootstrap", boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
})


test_that("'parallel' on a smoothed curve warns too", {
  skip_if_not_installed("MASS")
  s.wfns <- smooth(r.ndka)
  expect_warning(ci.auc(s.wfns, boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
  expect_warning(ci.se(s.wfns, boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
  expect_warning(ci.sp(s.wfns, boot.n = B, parallel = TRUE),
                 "Parallel processing is deprecated")
})


test_that("'parallel = FALSE' is silent", {
  expect_no_warning(ci.auc(r.wfns, method = "bootstrap", boot.n = B, parallel = FALSE))
  expect_no_warning(ci.se(r.wfns, boot.n = B, parallel = FALSE))
  # and the default is FALSE, so the plain call is silent as well
  expect_no_warning(ci.auc(r.wfns, method = "bootstrap", boot.n = B))
})


test_that("any non-FALSE 'parallel' warns, not just TRUE", {
  # Before 1.19 'parallel' was a logical, but accept anything truthy here: the
  # point is that a request to go parallel is never silently dropped.
  expect_warning(ci.auc(r.wfns, method = "bootstrap", boot.n = B, parallel = 2),
                 "Parallel processing is deprecated")
  expect_warning(ci.auc(r.wfns, method = "bootstrap", boot.n = B, parallel = "snow"),
                 "Parallel processing is deprecated")
  expect_warning(ci.auc(r.wfns, method = "bootstrap", boot.n = B, parallel = NULL),
                 "Parallel processing is deprecated")
})


test_that("'progress' warns and is ignored", {
  expect_warning(ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = "text"),
                 "Progress bars are deprecated")
  expect_warning(ci.se(r.wfns, boot.n = B, progress = "text"),
                 "Progress bars are deprecated")
  expect_warning(ci.sp(r.wfns, boot.n = B, progress = "text"),
                 "Progress bars are deprecated")
  expect_warning(ci.thresholds(r.wfns, boot.n = B, progress = "text"),
                 "Progress bars are deprecated")
  expect_warning(ci.coords(r.wfns, x = 0.5, input = "specificity", boot.n = B,
                           progress = "text"),
                 "Progress bars are deprecated")
  expect_warning(var(r.wfns, method = "bootstrap", boot.n = B, progress = "text"),
                 "Progress bars are deprecated")
  expect_warning(cov(r.wfns, r.ndka, method = "bootstrap", boot.n = B, progress = "text"),
                 "Progress bars are deprecated")
  expect_warning(roc.test(r.wfns, r.ndka, method = "bootstrap", boot.n = B,
                          progress = "text"),
                 "Progress bars are deprecated")
})


test_that("deprecated arguments do not change the result", {
  set.seed(42)
  with.arg <- suppressWarnings(
    ci.auc(r.wfns, method = "bootstrap", boot.n = 20, parallel = TRUE, progress = "text"))
  set.seed(42)
  without.arg <- ci.auc(r.wfns, method = "bootstrap", boot.n = 20)
  expect_equal(as.numeric(with.arg), as.numeric(without.arg))
})
