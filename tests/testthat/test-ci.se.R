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

test_that("ci.se gives the same result on percent and fraction curves", {
  # sp = 23/40 is a vertical segment: (23/40) * 100 != 57.5 in floating point
  controls <- c(1:23, 30:46)
  cases <- c(23.6, 23.7, 23.8, 50:55)
  r <- roc(controls = controls, cases = cases, quiet = TRUE)
  rp <- roc(controls = controls, cases = cases, percent = TRUE, quiet = TRUE)
  seed <- sample.int(1e6, 1)
  set.seed(seed)
  obtained <- ci.se(r, specificities = c(0.575, 0.6), boot.n = 50)
  set.seed(seed)
  obtained.percent <- ci.se(rp, specificities = c(57.5, 60), boot.n = 50)
  expect_equal(as.numeric(obtained.percent), as.numeric(obtained) * 100)
})
