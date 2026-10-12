library(pROC)
data(aSAH)

context("ci.thresholds")

# Only test whether ci.thresholds runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.

for (stratified in c(TRUE, FALSE)) {
  test_that("ci.threshold accepts thresholds=best", {
    n <- round(runif(1, 3, 9)) # keep boot.n small
    obtained <- ci.thresholds(r.wfns,
      thresholds = "best", boot.n = n,
      boot.stratified = stratified, conf.level = .91
    )
    expect_is(obtained, "ci.thresholds")
    expect_is(obtained, "ci")
    expect_equal(names(obtained), c("specificity", "sensitivity"))
    expect_equal(dim(obtained$specificity), c(1, 3))
    expect_equal(dim(obtained$sensitivity), c(1, 3))
    expect_equal(attr(obtained, "conf.level"), .91)
    expect_equal(attr(obtained, "boot.n"), n)
    expect_equal(colnames(obtained$specificity), c("4.5%", "50%", "95.5%"))
    expect_equal(colnames(obtained$sensitivity), c("4.5%", "50%", "95.5%"))
    expect_equal(attr(obtained, "boot.stratified"), stratified)
  })

  test_that("ci.threshold accepts thresholds=best", {
    n <- round(runif(1, 3, 9)) # keep boot.n small
    obtained <- ci.thresholds(r.ndka,
      thresholds = "local maximas", boot.n = n,
      boot.stratified = stratified, conf.level = .91
    )
    expected.thresholds <- coords(r.ndka, x = "l", ret = "t", transpose = FALSE)$threshold
    expect_is(obtained, "ci.thresholds")
    expect_is(obtained, "ci")
    expect_equal(names(obtained), c("specificity", "sensitivity"))
    expect_equal(dim(obtained$specificity), c(length(expected.thresholds), 3))
    expect_equal(dim(obtained$sensitivity), c(length(expected.thresholds), 3))
    expect_equal(attr(obtained, "conf.level"), .91)
    expect_equal(attr(obtained, "boot.n"), n)
    expect_equal(colnames(obtained$specificity), c("4.5%", "50%", "95.5%"))
    expect_equal(colnames(obtained$sensitivity), c("4.5%", "50%", "95.5%"))
    expect_equal(attr(obtained, "boot.stratified"), stratified)
  })

  test_that("ci.threshold accepts numeric thresholds", {
    n <- round(runif(1, 3, 9)) # keep boot.n small
    obtained <- ci.thresholds(r.ndka,
      thresholds = c(0.5, 0.2), boot.n = n,
      boot.stratified = stratified, conf.level = .91
    )
    expected.thresholds <- coords(r.ndka, x = "l", ret = "t", transpose = FALSE)$threshold
    expect_is(obtained, "ci.thresholds")
    expect_is(obtained, "ci")
    expect_equal(names(obtained), c("specificity", "sensitivity"))
    expect_equal(dim(obtained$specificity), c(2, 3))
    expect_equal(dim(obtained$sensitivity), c(2, 3))
    expect_equal(attr(obtained, "conf.level"), .91)
    expect_equal(attr(obtained, "boot.n"), n)
    expect_equal(colnames(obtained$specificity), c("4.5%", "50%", "95.5%"))
    expect_equal(colnames(obtained$sensitivity), c("4.5%", "50%", "95.5%"))
    expect_equal(attr(obtained, "boot.stratified"), stratified)
  })
}

test_that("ci.thresholds rejects character thresholds on a numeric ROC", {
  expect_error(
    ci.thresholds(r.ndka, thresholds = "10", boot.n = 3),
    "arg"
  )
  expect_error(
    ci.thresholds(r.ndka, thresholds = c("9", "10"), boot.n = 3),
    "'thresholds' of class character must be of length 1"
  )
})

test_that("ci.thresholds with NA thresholds gives a clear error", {
  msg <- "Missing values are not allowed in 'thresholds'."
  expect_error(ci.thresholds(r.s100b, thresholds = c(0.5, NA), boot.n = 2), msg, fixed = TRUE)
  expect_error(ci.thresholds(r.s100b, thresholds = NA, boot.n = 2), msg, fixed = TRUE)
  expect_error(ci.thresholds(r.s100b, thresholds = c(TRUE, NA, rep(FALSE, length(r.s100b$thresholds) - 2)), boot.n = 2), msg, fixed = TRUE)
})

test_that("ci.thresholds 'all' on a partial AUC curve skips the window limits", {
  # coords(r, "all") adds the interpolated window end with an NA threshold
  r <- roc(controls = c(1, 2, 3, 4, 5), cases = c(3, 4, 5, 6, 7), partial.auc = c(1, 0.7), quiet = TRUE)
  ci <- ci.thresholds(r, boot.n = 10, thresholds = "all")
  expect_identical(rownames(ci$sensitivity), c("4.5", "5.5", "6.5", "Inf"))
  expect_identical(rownames(ci$specificity), c("4.5", "5.5", "6.5", "Inf"))
  expect_false(anyNA(ci$sensitivity))
  expect_false(anyNA(ci$specificity))
  expect_output(print(ci), "6.5")
})

test_that("ci.thresholds gives a clear error with no threshold in the partial AUC window", {
  # specificities 0, 0.5, 0.5, 1, 1: none in [0.6, 0.8]
  r <- roc(controls = c(1, 3), cases = c(2, 4), partial.auc = c(0.8, 0.6), quiet = TRUE)
  for (thresholds in c("local maximas", "all", "best")) {
    expect_error(
      suppressWarnings(ci.thresholds(r, boot.n = 2, thresholds = thresholds)),
      "No threshold of the ROC curve in the partial AUC range"
    )
  }
})
