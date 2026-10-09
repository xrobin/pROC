library(pROC)
data(aSAH)

context("ci.coords")

test_that("ci.coords accepts threshold output with x=best", {
  expect_warning(
    expect_error(ci.coords(r.wfns, x = "best", input = "specificity", ret = c("threshold", "specificity", "sensitivity"), boot.n = 1), NA),
    "not available for ordered"
  )
})

test_that("ci.coords rejects threshold output except with x=best", {
  expect_error(ci.coords(r.s100b, x = 0.9, input = "specificity", ret = c("threshold", "specificity", "sensitivity"), boot.n = 1))
})

test_that("ci.coords accepts threshold output with x=best or if input was threshold", {
  expect_warning(
    expect_s3_class(ci.coords(r.wfns, x = "2", input = "threshold", ret = c("threshold", "specificity", "sensitivity"), boot.n = 1), "ci.coords"),
    "not available for ordered"
  )
  expect_warning(
    expect_s3_class(ci.coords(r.wfns, x = "best", ret = c("threshold", "specificity", "sensitivity"), boot.n = 1), "ci.coords"),
    "not available for ordered"
  )
})

test_that("ci.coords merges identical best points from empty ordered levels", {
  # Level 2 is unused: thresholds 2 and 3 are the same point in every replicate,
  # and clearly the best one, so that "stop" fails deterministically
  predictor <- ordered(c(rep("1", 20), rep("4", 2), rep("3", 20), rep("1", 2)), levels = c("1", "2", "3", "4"))
  response <- rep(c(0, 1), each = 22)
  r <- roc(response, predictor, quiet = TRUE)
  expect_s3_class(ci.coords(r, x = "best", boot.n = 3), "ci.coords")
  expect_warning(
    expect_s3_class(ci.coords(r, x = "best", ret = c("threshold", "specificity", "sensitivity"), boot.n = 3), "ci.coords"),
    "not available for ordered"
  )
  expect_error(ci.coords(r, x = "best", boot.n = 3, best.policy = "stop"), "More than one")
})

test_that("best.policy unique.stop still stops on distinct best points", {
  res <- data.frame(threshold = c(1.5, 2.5), specificity = c(0.5, 0.5), sensitivity = c(0.8, 0.8))
  expect_error(enforce.best.policy(res, "unique.stop"), "More than one distinct")
  res <- data.frame(specificity = c(0.5, 0.6), sensitivity = c(0.8, 0.7))
  expect_error(enforce.best.policy(res, "unique.stop"), "More than one distinct")
  res <- data.frame(specificity = c(0.5, 0.5), sensitivity = c(0.8, 0.8))
  expect_equal(enforce.best.policy(res, "unique.stop"), res[1, , drop = FALSE])
})

# Only test whether ci.coords runs and returns without error.
# Uses a very small number of iterations for speed
# Doesn't test whether the results are correct.
valid_coords_input <- coord.is.monotone <- c(
  "threshold", "sensitivity", "specificity", "tn", "tp", "fn", "fp", "tpr",
  "tnr", "fpr", "fnr", "1-specificity", "1-sensitivity", "recall"
)
for (input in valid_coords_input) {
  for (stratified in c(TRUE, FALSE)) {
    for (test.roc in list(r.s100b, smooth(r.s100b))) {
      # A smoothed curve has no thresholds, so "threshold" is not a valid
      # input for it: coords() rejects that combination and ci.coords() now
      # agrees instead of silently working from absent thresholds.
      if (input == "threshold" && methods::is(test.roc, "smooth.roc")) next
      context(sprintf("input: %s, stratified: %s, class: %s", input, stratified, class(test.roc)))
      test_that("ci.coords accepts one x and one ret", {
        skip_slow()
        obtained <- ci.coords(test.roc,
          x = 0.8, input = input, ret = "sp",
          boot.n = 3, conf.level = .91, boot.stratified = stratified
        )
        expect_equal(attr(obtained, "ret"), "specificity")
        expect_equal(names(obtained), attr(obtained, "ret"))
        for (ci.mat in obtained) {
          expect_equal(dim(ci.mat), c(1, 3))
          expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
        }
      })

      test_that("ci.coords accepts one x and multiple ret", {
        skip_slow()
        obtained <- ci.coords(test.roc,
          x = 0.8, input = input, ret = c("sp", "ppv", "tp", "1-sensitivity"),
          boot.n = 3, conf.level = .91, boot.stratified = stratified
        )
        expect_equal(attr(obtained, "ret"), c("specificity", "ppv", "tp", "1-sensitivity"))
        expect_equal(names(obtained), attr(obtained, "ret"))
        for (ci.mat in obtained) {
          expect_equal(dim(ci.mat), c(1, 3))
          expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
        }
      })

      test_that("ci.coords accepts multiple x and one ret", {
        skip_slow()
        obtained <- ci.coords(test.roc,
          x = c(0.8, 0.9), input = input, ret = "sp",
          boot.n = 3, conf.level = .91, boot.stratified = stratified
        )
        expect_equal(attr(obtained, "ret"), "specificity")
        expect_equal(names(obtained), attr(obtained, "ret"))
        for (ci.mat in obtained) {
          expect_equal(dim(ci.mat), c(2, 3))
          expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
        }
      })

      test_that("ci.coords accepts multiple x and ret", {
        skip_slow()
        obtained <- ci.coords(test.roc,
          x = c(0.9, 0.8), input = input, ret = c("sp", "ppv", "tp", "1-se"),
          boot.n = 3, conf.level = .91, boot.stratified = stratified
        )
        expect_equal(attr(obtained, "ret"), c("specificity", "ppv", "tp", "1-sensitivity"))
        expect_equal(names(obtained), attr(obtained, "ret"))
        for (ci.mat in obtained) {
          expect_equal(dim(ci.mat), c(2, 3))
          expect_equal(colnames(ci.mat), c("4.5%", "50%", "95.5%"))
        }
      })
    }
  }
}


test_that("ci.coords works on a smoothed curve with x = 'best'", {
  # ci.coords() resolved 'input' from a two-element default through
  # match.arg(several.ok = FALSE), so it failed for every x on a smoothed
  # curve unless 'input' was given explicitly:
  #   Error in match.arg(x, valid.args, several.ok = FALSE):
  #     'arg' must be of length 1
  # and once past that, x = "best" hit a second problem: the bootstrap worker
  # called coords.roc() on a smooth.roc, bypassing the smooth method that
  # fills the absent thresholds with NA.
  s <- smooth(r.s100b)

  set.seed(1)
  best <- ci.coords(s, "best", ret = "sensitivity", boot.n = 20)
  expect_s3_class(best, "ci.coords")
  # The bootstrap median should sit near the point estimate.
  point <- as.numeric(coords(s, "best", ret = "sensitivity"))
  expect_equal(as.numeric(best$sensitivity)[2], point, tolerance = 0.05)

  # A numeric x works without naming 'input', which defaults to specificity.
  set.seed(1)
  by.default <- ci.coords(s, 0.5, ret = "sensitivity", boot.n = 20)
  set.seed(1)
  explicit <- ci.coords(s, 0.5, input = "specificity", ret = "sensitivity", boot.n = 20)
  expect_equal(as.numeric(unlist(by.default)), as.numeric(unlist(explicit)))
})


test_that("ci.coords rejects input = 'threshold' on a smoothed curve", {
  # As coords() does: a smoothed curve has no thresholds.
  s <- smooth(r.s100b)
  expect_error(ci.coords(s, 0.5, input = "threshold", ret = "sp", boot.n = 3))
  expect_error(coords(s, 0.5, input = "threshold", ret = "sp"))
})
