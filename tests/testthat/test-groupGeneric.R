library(pROC)
data(aSAH)

test_that("Ops on auc and ci objects drop the attributes, including unary operators", {
  set.seed(42)
  a <- auc(r.s100b)
  ci1 <- ci.auc(r.s100b)
  expect_identical(-a, -as.numeric(a))
  expect_identical(+a, as.numeric(a))
  expect_identical(-ci1, -as.numeric(ci1))
  expect_identical(1 - a, 1 - as.numeric(a))
})

test_that("Ops between auc and ci objects return plain numbers without warning", {
  # R < 4.1.0 only dispatches Ops when both operands find the same method object;
  # Ops.auc and Ops.ci.auc are identical() but distinct once lazy-loaded.
  skip_if(getRversion() < "4.1")
  set.seed(42)
  a <- auc(r.s100b)
  ci1 <- ci.auc(r.s100b)
  expect_silent(d <- ci1 - a)
  expect_identical(d, as.numeric(ci1) - as.numeric(a))
  expect_silent(d <- a - ci1)
  expect_identical(d, as.numeric(a) - as.numeric(ci1))
})
