library(pROC)
data(aSAH)

test_that("paired venkatraman works as expected", {
  skip_slow()
  ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12)
  expect_venkatraman_htest(ht)
  expect_equal(ht$alternative, "two.sided")
  expect_equal(ht$method, "Venkatraman's test for two paired ROC curves")
  expect_equal(unname(ht$parameter), 12)
  # Test output
  ht$statistic <- c(E = 42)
  ht$p.value <- 0
  expect_known_output(print(ht), "print_output/roc.test-venkatraman.paired")
})

test_that("unpaired venkatraman works as expected", {
  skip_slow()
  expect_warning(ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, paired = FALSE), "paired")
  expect_venkatraman_htest(ht)
  expect_equal(ht$alternative, "two.sided")
  expect_equal(ht$method, "Venkatraman's test for two unpaired ROC curves")
  expect_equal(unname(ht$parameter), 12)
  # Test output
  ht$statistic <- c(E = 41)
  ht$p.value <- 0.548347196932
  expect_known_output(print(ht), "print_output/roc.test-venkatraman.unpaired")
})

test_that("non stratified venkatraman works as expected", {
  skip_slow()
  ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, boot.stratified = FALSE)
  expect_venkatraman_htest(ht)
  expect_equal(ht$alternative, "two.sided")
  expect_equal(ht$method, "Venkatraman's test for two paired ROC curves")
  expect_equal(unname(ht$parameter), 12)
  # Test output
  ht$statistic <- c(E = 43)
  ht$p.value <- 0.05
  expect_known_output(print(ht), "print_output/roc.test-venkatraman.unstratified")
})

test_that("non stratified, unpaired venkatraman works as expected", {
  skip_slow()
  expect_warning(ht <- roc.test(r.s100b, r.wfns, method = "venkatraman", boot.n = 12, boot.stratified = FALSE, paired = FALSE), "paired")
  expect_venkatraman_htest(ht)
  expect_equal(ht$alternative, "two.sided")
  expect_equal(ht$method, "Venkatraman's test for two unpaired ROC curves")
  expect_equal(unname(ht$parameter), 12)
  # Test output
  ht$statistic <- c(E = 43)
  ht$p.value <- 0.05
  expect_known_output(print(ht), "print_output/roc.test-venkatraman.unpaired.unstratified")
})

test_that("unpaired venkatraman statistic is symmetric with curves of different sizes", {
  r1 <- roc(c(0, 0, 0, 1, 1, 1), c(1, 2, 4, 3, 5, 6), quiet = TRUE)
  r2 <- roc(c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), c(1, 3, 2, 4:12), quiet = TRUE)
  t12 <- roc.test(r1, r2, method = "venkatraman", paired = FALSE, boot.n = 2)
  t21 <- roc.test(r2, r1, method = "venkatraman", paired = FALSE, boot.n = 2)
  expect_equal(t12$statistic, t21$statistic)
})
