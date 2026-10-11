library(pROC)
data(aSAH)

context("ci.se")

# Only test whether ci.se runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.


for (stratified in c(TRUE, FALSE)) {
  for (test.roc in list(r.s100b, smooth(r.s100b))) {
    test_that("ci.se with default specificities", {
      n <- round(runif(1, 3, 9)) # keep boot.n small
      obtained <- ci.se(test.roc,
        boot.n = n,
        boot.stratified = stratified, conf.level = .91
      )
      expect_is(obtained, "ci.se")
      expect_is(obtained, "ci")
      expect_equal(dim(obtained), c(11, 3))
      expect_equal(attr(obtained, "conf.level"), .91)
      expect_equal(attr(obtained, "boot.n"), n)
      expect_equal(colnames(obtained), c("4.5%", "50%", "95.5%"))
      expect_equal(attr(obtained, "boot.stratified"), stratified)
    })

    test_that("ci.se accepts one specificity", {
      n <- round(runif(1, 3, 9)) # keep boot.n small
      obtained <- ci.se(test.roc,
        specificities = 0.9, boot.n = n,
        boot.stratified = stratified, conf.level = .91
      )
      expect_is(obtained, "ci.se")
      expect_is(obtained, "ci")
      expect_equal(dim(obtained), c(1, 3))
      expect_equal(attr(obtained, "conf.level"), .91)
      expect_equal(attr(obtained, "boot.n"), n)
      expect_equal(colnames(obtained), c("4.5%", "50%", "95.5%"))
      expect_equal(attr(obtained, "boot.stratified"), stratified)
    })
  }
}

test_that("ci.se names itself in the multiclass error", {
  mc <- multiclass.roc(aSAH$gos6, aSAH$s100b, quiet = TRUE)
  expect_error(ci.se(mc), "'ci.se' not available for multiclass ROC curves.", fixed = TRUE)
})

test_that("ci.se of a smoothed curve gives an NA interval when no resample can be smoothed", {
  s <- smooth(r.s100b)
  # Every resample fails to smooth
  s$smoothing.args$method <- "density"
  s$smoothing.args$density.controls <- function(x, ...) stop("not smoothable")
  expect_warning(obtained <- ci.se(s, specificities = c(0.5, 0.9), boot.n = 3, progress = "none"), "NA value")
  expect_equal(dim(obtained), c(2, 3))
  expect_true(all(is.na(obtained)))
})

test_that("ci.se of a smoothed curve with a corrected partial AUC does not warn for each resample", {
  # smooth() recomputed the corrected pAUC of every resample, warning when below the diagonal
  resp <- rep(0:1, each = 8)
  pred <- c(0.1, 0.5, 0.9, 1.3, 1.6, 2.0, 2.4, 2.9, 0.3, 0.8, 1.2, 1.7, 2.2, 2.6, 3.1, 3.4)
  s <- smooth(roc(resp, pred, partial.auc = c(1, 0.9), partial.auc.correct = TRUE, quiet = TRUE))
  warnings <- character(0)
  withCallingHandlers(ci.se(s, boot.n = 50), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("Partial AUC correction not defined", warnings)))
})
