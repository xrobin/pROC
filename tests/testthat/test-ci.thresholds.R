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

test_that("ci.thresholds with NA thresholds gives NA rows", {
  set.seed(42)
  ci <- ci.thresholds(r.s100b, thresholds = c(0.5, NA), boot.n = 10)
  set.seed(42)
  ci1 <- ci.thresholds(r.s100b, thresholds = 0.5, boot.n = 10)
  expect_identical(ci$specificity[1, ], ci1$specificity[1, ])
  expect_identical(ci$sensitivity[1, ], ci1$sensitivity[1, ])
  expect_identical(unname(ci$specificity[2, ]), rep(NA_real_, 3))
  expect_identical(unname(ci$sensitivity[2, ]), rep(NA_real_, 3))
  expect_output(print(ci), "NA")
  # A logical NA selects a NA threshold, as in base R indexing
  ci <- ci.thresholds(r.s100b, thresholds = c(TRUE, NA, rep(FALSE, length(r.s100b$thresholds) - 2)), boot.n = 2)
  expect_identical(as.vector(attr(ci, "thresholds")), c(r.s100b$thresholds[1], NA))
  expect_identical(unname(ci$sensitivity[2, ]), rep(NA_real_, 3))
})
