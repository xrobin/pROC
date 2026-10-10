library(pROC)
data(aSAH)

context("power.roc.test")

# define variables shared among multiple tests here

test_that("power.roc.test basic function", {
  res <- power.roc.test(r.s100b)
  expect_equal(as.numeric(res$auc), as.numeric(r.s100b$auc))
  expect_equal(res$ncases, length(r.s100b$cases))
  expect_equal(res$ncontrols, length(r.s100b$controls))
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.9904833, tolerance = 0.000001)
})

test_that("power.roc.test with percent works", {
  res <- power.roc.test(r.s100b.percent)
  expect_equal(as.numeric(res$auc), as.numeric(r.s100b$auc))
  expect_equal(res$ncases, length(r.s100b$cases))
  expect_equal(res$ncontrols, length(r.s100b$controls))
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.9904833, tolerance = 0.000001)
})

test_that("power.roc.test with given auc function", {
  res <- power.roc.test(ncases = 41, ncontrols = 72, auc = 0.73, sig.level = 0.05)
  expect_equal(as.numeric(res$auc), 0.73)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.9897453, tolerance = 0.000001)
})


test_that("power.roc.test sig.level can be omitted", {
  res <- power.roc.test(ncases = 41, ncontrols = 72, auc = 0.73)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.9897453, tolerance = 0.000001)
})

test_that("power.roc.test can determine ncases & ncontrols", {
  res <- power.roc.test(auc = r.s100b$auc, sig.level = 0.05, power = 0.95, kappa = 1.7)
  expect_equal(as.numeric(res$auc), as.numeric(r.s100b$auc))
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$power, 0.95)
  expect_equal(res$ncases, 29.29764, tolerance = 0.000001)
  expect_equal(res$ncontrols, 49.806, tolerance = 0.000001)
})

test_that("power.roc.test can determine sig.level", {
  res <- power.roc.test(ncases = 41, ncontrols = 72, auc = 0.73, power = 0.95, sig.level = NULL)
  expect_equal(as.numeric(res$auc), 0.73)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(res$power, 0.95)
  expect_equal(res$sig.level, 0.009238584, tolerance = 0.000001)
})

test_that("power.roc.test can determine AUC", {
  res <- power.roc.test(ncases = 41, ncontrols = 72, sig.level = 0.05, power = 0.95)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(res$power, 0.95)
  expect_equal(res$sig.level, 0.05)
  expect_equal(as.numeric(res$auc), 0.6961054, tolerance = 0.000001)
})

test_that("power.roc.test can take 2 ROC curves with DeLong variance", {
  res <- power.roc.test(r.ndka, r.wfns)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.6872839, tolerance = 0.000001)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})


test_that("power.roc.test can take 2 percent ROC curves with DeLong variance", {
  res <- power.roc.test(r.ndka.percent, r.wfns.percent)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.6872839, tolerance = 0.000001)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test can take 2 ROC curves with Obuchowski variance", {
  res <- power.roc.test(r.ndka, r.wfns, method = "obuchowski")
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.7869306, tolerance = 0.000001)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test ncases/ncontrols can take 2 ROC curves with DeLong variance", {
  res <- power.roc.test(r.ndka, r.wfns, power = 0.9)
  expect_equal(res$ncases, 67.55038, tolerance = 0.000001)
  expect_equal(res$ncontrols, 118.6251, tolerance = 0.000001)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test ncases/ncontrols can take 2 ROC curves with Obuchowski variance", {
  res <- power.roc.test(r.ndka, r.wfns, power = 0.9, method = "obuchowski")
  expect_equal(res$ncases, 55.80446, tolerance = 0.000001)
  expect_equal(res$ncontrols, 97.99807, tolerance = 0.000001)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test sig.level can take 2 ROC curves with DeLong variance", {
  res <- power.roc.test(r.ndka, r.wfns, power = 0.9, sig.level = NULL)
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.1982034, tolerance = 0.000001)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test sig.level can take 2 ROC curves with Obuchowski variance", {
  res <- power.roc.test(r.ndka, r.wfns, power = 0.9, sig.level = NULL, method = "obuchowski")
  expect_equal(res$ncases, 41)
  expect_equal(res$ncontrols, 72)
  expect_equal(as.numeric(res$auc1), as.numeric(r.ndka$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.wfns$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.1308799, tolerance = 0.000001)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test defaults to bootstrap with partial AUC", {
  set.seed(42)
  res.default <- power.roc.test(r.wfns.partial, r.ndka.partial, power = 0.9, boot.n = 100)
  set.seed(42)
  res.bootstrap <- power.roc.test(r.wfns.partial, r.ndka.partial, power = 0.9, boot.n = 100, method = "bootstrap")
  expect_equal(res.default, res.bootstrap)
  expect_true(is.finite(res.default$ncases))
  expect_equal(as.numeric(res.default$auc1), as.numeric(r.wfns.partial$auc))

  # All three modes
  set.seed(42)
  expect_true(is.finite(power.roc.test(r.wfns.partial, r.ndka.partial, boot.n = 100)$power))
  set.seed(42)
  expect_true(is.finite(power.roc.test(r.wfns.partial, r.ndka.partial, power = 0.9, sig.level = NULL, boot.n = 100)$sig.level))
})

test_that("power.roc.test refuses explicit DeLong with partial AUC", {
  expect_error(
    power.roc.test(r.wfns.partial, r.ndka.partial, power = 0.9, method = "delong"),
    "DeLong method is not supported for partial AUC"
  )
})

test_that("power.roc.test still defaults to DeLong with full AUC", {
  expect_equal(
    power.roc.test(r.ndka, r.wfns, power = 0.9),
    power.roc.test(r.ndka, r.wfns, power = 0.9, method = "delong")
  )
})


test_that("power.roc.test works with partial AUC", {
  r.wfns.partial <- roc(aSAH$outcome, aSAH$wfns, quiet = TRUE, partial.auc = c(1, 0.9))
  r.ndka.partial <- roc(aSAH$outcome, aSAH$ndka, quiet = TRUE, partial.auc = c(1, 0.9))
  res <- power.roc.test(r.wfns.partial, r.ndka.partial, power = 0.9, method = "obuchowski")

  expect_equal(res$ncases, 227.7338, tolerance = 0.000001)
  expect_equal(res$ncontrols, 399.9227, tolerance = 0.000001)
  expect_equal(as.numeric(res$auc1), as.numeric(r.wfns.partial$auc))
  expect_equal(as.numeric(res$auc2), as.numeric(r.ndka.partial$auc))
  expect_equal(res$power, 0.9)
  expect_equal(res$sig.level, 0.05)
  expect_equal(res$alternative, "two.sided")
})

test_that("power.roc.test works with binormal parameters", {
  ob.params <- list(
    A1 = 2.6, B1 = 1, A2 = 1.9, B2 = 1, rn = 0.6, ra = 0.6, FPR11 = 0,
    FPR12 = 0.2, FPR21 = 0, FPR22 = 0.2, delta = 0.037
  )

  res1 <- power.roc.test(ob.params, power = 0.8, sig.level = 0.05)
  expect_equal(res1$ncases, 119.7869, tolerance = 0.000001)
  expect_equal(res1$ncontrols, 119.7869, tolerance = 0.000001)
  expect_equal(res1$power, 0.8)
  expect_equal(res1$sig.level, 0.05)

  res2 <- power.roc.test(ob.params, power = 0.8, sig.level = NULL, ncases = 107)
  expect_equal(res2$ncases, 107)
  expect_equal(res2$ncontrols, 107)
  expect_equal(res2$power, 0.8)
  expect_equal(res2$sig.level, 0.07258085, tolerance = 0.000001)

  res3 <- power.roc.test(ob.params, power = NULL, sig.level = 0.05, ncases = 107)
  expect_equal(res3$ncases, 107)
  expect_equal(res3$ncontrols, 107)
  expect_equal(res3$sig.level, 0.05)
  expect_equal(res3$power, 0.7605865, tolerance = 0.000001)
})

## With only binormal parameters given
# From example 2 of Obuchowski and McClish, 1997.


test_that("power.roc.test returns correct results from litterature", {
  context("Check results in Obuchowski 2004 Table 4")
  # Note: the table reports at least 10 in each cell, and adapts
  # the complement value to match kappa. So in 0.25/0.95 we have 10/40
  # although both values are < 10.
  # Note2: some values don't match exactly, specifically
  # expected.ncases[0.5, 0.6] and expected.ncontrols[4, 0.6]
  # are off by 1 (< 1%).
  kappas <- c(0.25, 0.5, 1, 2, 4)
  thetas <- c(0.6, 0.7, 0.8, 0.9, 0.95)
  expected.ncontrols <- matrix(
    c(
      84,     20,     10,     10,     10,
      101,    25,     10,     10,     10,
      135,    33,     14,     10,     10,
      203,    50,     21,     20,     20,
      # 339,	84,	40,	40,	40
      340,    84,     40,     40,     40 # Fixed
    ),
    nrow = 5, byrow = TRUE,
    dimnames = list(kappas, thetas)
  )
  expected.ncases <- matrix(
    c(
      334,    80,     40,     40,     40,
      # 201,	49,	20,	20,	20,
      202,    49,     20,     20,     20, # Fixed
      135,    33,     14,     10,     10,
      102,    25,     11,     10,     10,
      85,     21,     10,     10,     10
    ),
    nrow = 5, byrow = TRUE,
    dimnames = list(kappas, thetas)
  )

  for (kappa in kappas) {
    for (theta in thetas) {
      context(sprintf("kappa: %s, theta: %s", kappa, theta))
      pr <- power.roc.test(auc = theta, sig.level = 0.05, power = 0.9, kappa = kappa, alternative = "one.sided")
      expect_equal(max(10, ifelse(ceiling(pr$ncases) < 10, 10, 0) * kappa, ceiling(pr$ncontrols)), expected.ncontrols[as.character(kappa), as.character(theta)])
      expect_equal(max(10, ifelse(ceiling(pr$ncontrols) < 10, 10, 0) / kappa, ceiling(pr$ncases)), expected.ncases[as.character(kappa), as.character(theta)])
    }
  }
})


test_that("kappa works with a single ROC curve", {
  # kappa from data
  res <- power.roc.test(r.s100b, sig.level = 0.05, power = 0.9)
  expect_equal(res$ncases, 23.5598674)
  expect_equal(res$ncontrols, 41.3734257)
  expect_equal(res$ncases / res$ncontrols, length(r.s100b$cases) / length(r.s100b$controls))
  # set kappa
  res <- power.roc.test(r.s100b, sig.level = 0.05, power = 0.9, kappa = 1)
  expect_equal(res$ncases, 29.5697422)
  expect_equal(res$ncontrols, 29.5697422)
})


test_that("kappa works with two ROC curves", {
  # kappa from data
  res <- power.roc.test(r.s100b, r.ndka, sig.level = 0.05, power = 0.9)
  expect_equal(res$ncases, 210.7168158)
  expect_equal(res$ncontrols, 370.0392862)
  expect_equal(res$ncases / res$ncontrols, length(r.s100b$cases) / length(r.s100b$controls))
  # set kappa
  res <- power.roc.test(r.s100b, r.ndka, sig.level = 0.05, power = 0.9, kappa = 1)
  expect_equal(res$ncases, 210.7168158)
  expect_equal(res$ncases, 210.7168158)
  # ...
})

test_that("power.roc.test with bootstrap uses the variance of each AUC", {
  # The variances must be computed over the replicates (rows of
  # resampled.values), not over the two curves of a single replicate.
  zalpha <- qnorm(1 - 0.05 / 2)
  zbeta <- qnorm(0.9)
  delta <- as.numeric(r.s100b$auc - r.ndka$auc)
  n <- length(r.s100b$cases)
  set.seed(42)
  res <- power.roc.test(r.s100b, r.ndka, power = 0.9, method = "bootstrap", boot.n = 100)
  set.seed(42)
  cv <- cov(r.s100b, r.ndka, method = "bootstrap", boot.n = 100, boot.return = TRUE)
  rv <- attr(cv, "resampled.values")
  var1 <- var(rv[1, ]) * n
  var2 <- var(rv[2, ]) * n
  cov12 <- as.numeric(cv) * n
  cov0 <- cov12 * sqrt(var1 / var2)
  v0 <- 2 * var1 - 2 * cov0
  va <- var1 + var2 - 2 * cov12
  expect_equal(res$ncases, (zalpha * sqrt(v0) + zbeta * sqrt(va))^2 / delta^2)
})

test_that("power.roc.test re-pairs curves with NAs at different positions", {
  x1 <- aSAH$s100b
  x2 <- aSAH$ndka
  x1[c(3, 50)] <- NA
  x2[c(10, 90)] <- NA
  r1 <- roc(aSAH$outcome, x1, quiet = TRUE)
  r2 <- roc(aSAH$outcome, x2, quiet = TRUE)
  ok <- !is.na(x1) & !is.na(x2)
  q1 <- roc(aSAH$outcome[ok], x1[ok], quiet = TRUE)
  q2 <- roc(aSAH$outcome[ok], x2[ok], quiet = TRUE)
  numeric_fields <- function(res) {
    unlist(res[c("ncases", "ncontrols", "auc1", "auc2", "sig.level", "power")])
  }
  for (method in c("delong", "obuchowski")) {
    expect_equal(
      numeric_fields(power.roc.test(r1, r2, method = method)),
      numeric_fields(power.roc.test(q1, q2, method = method))
    )
    expect_equal(
      numeric_fields(power.roc.test(r1, r2, method = method, power = 0.9)),
      numeric_fields(power.roc.test(q1, q2, method = method, power = 0.9))
    )
    expect_equal(
      numeric_fields(power.roc.test(r1, r2, method = method, power = 0.9, sig.level = NULL)),
      numeric_fields(power.roc.test(q1, q2, method = method, power = 0.9, sig.level = NULL))
    )
  }
  set.seed(42)
  res.r <- power.roc.test(r1, r2, method = "bootstrap", boot.n = 20)
  set.seed(42)
  res.q <- power.roc.test(q1, q2, method = "bootstrap", boot.n = 20)
  expect_equal(numeric_fields(res.r), numeric_fields(res.q))
})

test_that("power.roc.test accepts the auc of a percent ROC curve", {
  res <- power.roc.test(auc = r.s100b.percent$auc, ncases = 41, ncontrols = 72)
  expect_equal(as.numeric(res$auc), as.numeric(r.s100b$auc))
  expect_equal(res$power, 0.9904833, tolerance = 0.000001)
  res <- power.roc.test(auc = r.s100b.percent$auc, sig.level = 0.05, power = 0.95, kappa = 1.7)
  expect_equal(res$ncases, 29.29764, tolerance = 0.000001)
  expect_equal(res$ncontrols, 49.806, tolerance = 0.000001)
})

test_that("power.roc.test checks the range of auc", {
  expect_error(power.roc.test(auc = 73, ncases = 41, ncontrols = 72), "'auc' must range from 0 to 1")
  expect_error(power.roc.test(auc = -0.1, ncases = 41, ncontrols = 72), "'auc' must range from 0 to 1")
})

test_that("power.roc.test accepts a vector of AUCs", {
  res <- power.roc.test(auc = c(0.7, 0.8, 0.9), power = 0.9)
  expect_length(res$ncases, 3)
  expect_equal(res$ncases[2], power.roc.test(auc = 0.8, power = 0.9)$ncases)
  expect_error(power.roc.test(auc = c(0.7, 1.2), power = 0.9), "must range from 0 to 1")
})

test_that("power.roc.test with binormal parameters accepts FPR bounds in any order", {
  base <- list(A1 = 2.6, B1 = 1, A2 = 1.9, B2 = 1, rn = 0.6, ra = 0.6, delta = 0.037)
  lower.upper <- c(base, list(FPR11 = 0, FPR12 = 0.2, FPR21 = 0, FPR22 = 0.2))
  upper.lower <- c(base, list(FPR11 = 0.2, FPR12 = 0, FPR21 = 0.2, FPR22 = 0))
  mixed <- c(base, list(FPR11 = 0.2, FPR12 = 0, FPR21 = 0, FPR22 = 0.2))
  expected <- power.roc.test(lower.upper, power = 0.8)$ncases
  expect_equal(power.roc.test(upper.lower, power = 0.8)$ncases, expected)
  expect_equal(power.roc.test(mixed, power = 0.8)$ncases, expected)
  expect_equal(power.roc.test(mixed, ncases = 107)$power, power.roc.test(lower.upper, ncases = 107)$power)
  expect_equal(
    power.roc.test(mixed, ncases = 107, power = 0.8, sig.level = NULL)$sig.level,
    power.roc.test(lower.upper, ncases = 107, power = 0.8, sig.level = NULL)$sig.level
  )
})

test_that("power.roc.test with binormal parameters requires all four FPR bounds", {
  base <- list(A1 = 2.6, B1 = 1, A2 = 1.9, B2 = 1, rn = 0.6, ra = 0.6, delta = 0.037)
  expect_error(power.roc.test(c(base, list(FPR11 = 0.2, FPR12 = 0)), power = 0.8), "FPR21, FPR22")
  expect_error(power.roc.test(c(base, list(FPR11 = 0.2, FPR12 = 0, FPR21 = 0)), power = 0.8), "FPR22")
})

test_that("power.roc.test refuses a partial AUC with one ROC curve", {
  expect_error(power.roc.test(r.s100b.partial), "only available for the full AUC")
  expect_error(power.roc.test(r.s100b.partial, power = 0.9), "only available for the full AUC")
  expect_error(power.roc.test(r.s100b.percent.partial1), "only available for the full AUC")
  expect_error(
    power.roc.test(r.s100b, reuse.auc = FALSE, partial.auc = c(1, 0.8)),
    "only available for the full AUC"
  )
  # Full AUC still works, also when recomputed
  expect_equal(power.roc.test(r.s100b, reuse.auc = FALSE)$power, power.roc.test(r.s100b)$power)
})

test_that("Obuchowski variance and covariance use the binormal A and B parameters", {
  # Binormal data: controls ~ N(0, 1), cases ~ N(2, 2), so that
  # qnorm(TPR) = A + B * qnorm(FPR) with A = 2 / 2 = 1 and B = 1 / 2.
  controls <- qnorm(ppoints(60))
  cases <- 2 + 2 * qnorm(ppoints(60))
  r <- roc(controls = controls, cases = cases, direction = "<", quiet = TRUE)
  rp <- roc(controls = controls, cases = cases, direction = "<", quiet = TRUE, partial.auc = c(1, 0.8))
  # Ratios to the formulas with the true parameters (relative tolerance)
  expect_equal(
    var(r, method = "obuchowski") / (pROC:::var_params_obuchowski(1, 0.5, 1) / 60),
    1,
    tolerance = 0.05
  )
  expect_equal(
    var(rp, method = "obuchowski") / (pROC:::var_params_obuchowski(1, 0.5, 1, 0.2, 0) / 60),
    1,
    tolerance = 0.05
  )
  expect_equal(
    cov(rp, rp, method = "obuchowski") / (pROC:::cov_params_obuchowski(1, 0.5, 1, 0.5, 1, 1, 1, 0.2, 0, 0.2, 0) / 60),
    1,
    tolerance = 0.05
  )
})

test_that("power.roc.test uses the covariance under the null hypothesis with Obuchowski", {
  # Obuchowski & McClish (1997) / DESIGNROC: V0 = 2 * V(A1, B1) - 2 * C(A1, B1, A1, B1)
  p <- list(A1 = 1.5, B1 = 1, A2 = 1, B2 = 1, rn = 0.5, ra = 0.5, delta = pnorm(1.5 / sqrt(2)) - pnorm(1 / sqrt(2)))
  v1 <- pROC:::var_params_obuchowski(p$A1, p$B1, 1)
  v2 <- pROC:::var_params_obuchowski(p$A2, p$B2, 1)
  c0 <- pROC:::cov_params_obuchowski(p$A1, p$B1, p$A1, p$B1, p$rn, p$ra, 1)
  ca <- pROC:::cov_params_obuchowski(p$A1, p$B1, p$A2, p$B2, p$rn, p$ra, 1)
  v0 <- 2 * v1 - 2 * c0
  va <- v1 + v2 - 2 * ca
  expected <- (qnorm(0.975) * sqrt(v0) + qnorm(0.8) * sqrt(va))^2 / p$delta^2
  expect_equal(power.roc.test(p, power = 0.8)$ncases, expected)
  expect_equal(expected, 118.5139, tolerance = 1e-6) # DESIGNROC: 118.535 (single precision)
  # Same with ROC curves: null covariance from the binormal parameters of roc1
  res <- power.roc.test(r.ndka, r.wfns, power = 0.9, method = "obuchowski")
  n <- length(r.ndka$cases)
  v1 <- var(r.ndka, method = "obuchowski") * n
  v2 <- var(r.wfns, method = "obuchowski") * n
  ca <- cov(r.ndka, r.wfns, method = "obuchowski") * n
  c0 <- pROC:::cov0.roc.obuchowski(r.ndka, r.wfns)
  v0 <- 2 * v1 - 2 * c0
  va <- v1 + v2 - 2 * ca
  delta <- as.numeric(r.ndka$auc - r.wfns$auc)
  expect_equal(res$ncases, (qnorm(0.975) * sqrt(v0) + qnorm(0.9) * sqrt(va))^2 / delta^2)
})

test_that("power.roc.test with DeLong uses the AUC correlation under the null hypothesis", {
  # V0 = 2 * var1 * (1 - cor12): both AUCs with the variance of roc1
  res <- power.roc.test(r.s100b, r.ndka, power = 0.9)
  n <- length(r.s100b$cases)
  var1 <- var(r.s100b) * n
  var2 <- var(r.ndka) * n
  cov12 <- cov(r.s100b, r.ndka) * n
  v0 <- 2 * var1 * (1 - cov12 / sqrt(var1 * var2))
  va <- var1 + var2 - 2 * cov12
  delta <- as.numeric(r.s100b$auc - r.ndka$auc)
  expect_equal(res$ncases, (qnorm(0.975) * sqrt(v0) + qnorm(0.9) * sqrt(va))^2 / delta^2)
})
