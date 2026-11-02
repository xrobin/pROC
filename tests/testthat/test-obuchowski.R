library(pROC)

# Calibration of the Obuchowski variance against Obuchowski and McClish
# (1997), Stat Med 16:1529-1542. The paper's tables print 3 decimals: the
# expectations are checked within 0.001.

# Exactly binormal data with qnorm(TPR) = A + B * qnorm(FPR): controls
# N(0, 1) and cases N(A / B, 1 / B^2), on the same 1000 quantiles
binormal.roc <- function(A, B, partial.auc = FALSE) {
  controls <- qnorm(ppoints(1000))
  roc(controls = controls, cases = A / B + controls / B, direction = "<", partial.auc = partial.auc, quiet = TRUE)
}

# Full area, and the FPR ranges 0-.05, 0-.2 and .05-.1 as specificity ranges
obuchowski.ranges <- list(FALSE, c(1, .95), c(1, .8), c(.95, .9))

# SE of the (partial) AUC with 100 cases and 100 controls
obuchowski.se100 <- function(A, B) {
  vapply(obuchowski.ranges, function(pa) {
    r <- binormal.roc(A, B, pa)
    sqrt(var(r, method = "obuchowski") * length(r$cases) / 100)
  }, numeric(1))
}

test_that("binormal A and B are recovered from exactly binormal data", {
  for (AB in list(c(1.19, 1), c(0.89, 1 / 3))) {
    params <- pROC:::binormal_params_obuchowski(smooth(binormal.roc(AB[1], AB[2]), method = "binormal")$model)
    expect_equal(unname(params), AB, tolerance = 0.01, scale = 1)
  }
})

test_that("var obuchowski matches the Equation (4) columns of Tables II and III", {
  # Table II, A = 1.19, B = 1; full area, FPR 0-.05, 0-.2, .05-.1
  expect_equal(obuchowski.se100(1.19, 1), c(0.035, 0.003, 0.015, 0.004), tolerance = 0.001, scale = 1)
  # Table III, A = 0.89, B = 1/3
  expect_equal(obuchowski.se100(0.89, 1 / 3), c(0.033, 0.003, 0.009, 0.002), tolerance = 0.001, scale = 1)
})

test_that("var_params_obuchowski matches Table IV", {
  # theta = 0.8, N = 200 patients, N_A = 200 / (1 + R) cases. The partial
  # areas take the (upper, lower) FPR bounds
  table.iv <- function(R, B) {
    A <- qnorm(0.8) * sqrt(1 + B^2)
    n.cases <- 200 / (1 + R)
    v <- c(
      pROC:::var_params_obuchowski(A, B, R),
      pROC:::var_params_obuchowski(A, B, R, 0.05, 0),
      pROC:::var_params_obuchowski(A, B, R, 0.2, 0),
      pROC:::var_params_obuchowski(A, B, R, 0.1, 0.05)
    )
    sqrt(v / n.cases)
  }
  expect_equal(table.iv(0.5, 0.5), c(0.031, 0.003, 0.011, 0.003), tolerance = 0.001, scale = 1)
  expect_equal(table.iv(5, 2), c(0.048, 0.002, 0.024, 0.006), tolerance = 0.001, scale = 1)
})
