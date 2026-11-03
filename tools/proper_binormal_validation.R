## Validation of smooth(method = "binormal", fit = "proper") against the
## published results of the "proper" binormal model (issue #152).
##
## Maintainer-only. Run from the package root with the version of pROC to
## check already installed (and mvtnorm):
##
##     Rscript tools/proper_binormal_validation.R [B] [cores]
##
## B: datasets per simulated population (default 2000); cores: for
## parallel::mclapply (default 4; use 1 on Windows). With B = 2000 it takes
## about 45 minutes on 4 cores. The results and their discussion are in
## tools/proper_binormal_validation.md.
##
## Sources:
## - Metz CE, Pan X (1999) "Proper" binormal ROC curves: theory and
##   maximum-likelihood estimation. J Math Psychol 43:1-33 (M&P).
## - Pesce LL, Metz CE (2007) Reliable and computationally efficient
##   maximum-likelihood estimation of "proper" binormal ROC curves. Acad
##   Radiol 14:814-829 (P&M).
## - Owen DB (1956) Tables for computing bivariate normal probabilities.
##   Ann Math Stat 27:1075-1090.
library(pROC)
library(parallel)
options(width = 200)

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 2000L
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
RNGkind("L'Ecuyer-CMRG")
set.seed(20261103)

cat("pROC", as.character(packageVersion("pROC")), "| mvtnorm", as.character(packageVersion("mvtnorm")),
  "|", R.version.string, "| B =", B, "\n\n")

## 1. AUC (M&P p. 12, Eq. 36) ---------------------------------------------
cat("== 1. AUC of Eq. 36\n")
# P&M p. 826, Table 1: population AUC, da and c (6 significant digits)
table1 <- data.frame(
  auc = c(0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.95, 0.95, 0.95),
  shape = rep(c("high", "medium", "symmetric"), 3),
  da = c(0.1, 0.358, 0.358287, 1.24365, 1.46571, 1.46574, 0.0999998, 2.324, 2.32618),
  c = c(-0.154357, -0.046907, 0, -0.428571, -0.14, 0, -0.853738, -0.292369, 0)
)
table1$eq36 <- mapply(pROC:::proper_binormal_auc, table1$da, table1$c)
print(table1, digits = 8)
cat(sprintf("Largest difference to Table 1: %.1e\n", max(abs(table1$eq36 - table1$auc))))
# da = 0: F(0, 0; rho) = 1/4 + arcsin(rho) / (2 pi) (Owen p. 1078), so the AUC is
# 1 + arcsin(rho) / pi with rho = -(1 - c^2) / (1 + c^2)
cs <- c(-0.99, -0.5, -1e-3, -1e-6, 1e-6, 1e-3, 0.5, 0.99)
closed <- 1 + asin(-(1 - cs^2) / (1 + cs^2)) / pi
cat(sprintf("Largest difference to the closed form at da = 0: %.1e\n\n",
  max(abs(vapply(cs, pROC:::proper_binormal_auc, numeric(1), da = 0) - closed))))

## 2. Continuous data: P&M Tables 2-4 (pp. 826-827) -------------------------
# P&M p. 822: x from the conventional binormal model with a and b (negatives
# N(0, 1), positives N(a / b, 1 / b)), and y = -[(bx - a)^2 - x^2] / 2, the
# decision variable of the PBM (M&P Eq. 11); a and b from da and c with
# M&P p. 23. At c = 0 (b = 1), y is linear in x: x is used directly.
sim.pbm <- function(n, da, c) {
  b <- (1 + c) / (1 - c)
  a <- da * sqrt(1 + c^2) / (1 - c)
  y <- if (c == 0) identity else function(x) -((b * x - a)^2 - x^2) / 2
  list(controls = y(rnorm(n)), cases = y(rnorm(n, a / b, 1 / b)))
}
fit.auc <- function(d) {
  r <- roc(controls = d$controls, cases = d$cases, direction = "<", quiet = TRUE)
  c(
    proper = tryCatch(as.numeric(auc(smooth(r, fit = "proper"))), error = function(e) NA_real_),
    lm = tryCatch(as.numeric(auc(smooth(r))), error = function(e) NA_real_),
    empirical = as.numeric(auc(r))
  )
}
# Published means and empirical SDs at 50 + 50 cases
published <- data.frame(
  auc = c(0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.95, 0.95, 0.95),
  shape = rep(c("high", "medium", "symmetric"), 3),
  mean = c(0.617, 0.610, 0.608, 0.860, 0.852, 0.852, 0.957, 0.951, 0.951),
  sd = c(0.049, 0.050, 0.051, 0.034, 0.038, 0.038, 0.018, 0.021, 0.019)
)
cat("== 2. P&M Tables 2-4, 50 + 50 cases\n")
rows <- list()
for (i in seq_len(nrow(table1))) {
  res <- do.call(rbind, mclapply(seq_len(B), function(j) fit.auc(sim.pbm(50, table1$da[i], table1$c[i])), mc.cores = cores))
  p <- res[, "proper"]
  se <- sd(p, na.rm = TRUE) / sqrt(sum(!is.na(p)))
  rows[[i]] <- data.frame(
    auc = table1$auc[i], shape = table1$shape[i],
    mean = mean(p, na.rm = TRUE), se.mean = se, published = published$mean[i],
    z = (mean(p, na.rm = TRUE) - published$mean[i]) / se,
    sd = sd(p, na.rm = TRUE), published.sd = published$sd[i],
    bias = mean(p, na.rm = TRUE) - table1$auc[i],
    failures = sum(is.na(p)),
    lm.bias = mean(res[, "lm"], na.rm = TRUE) - table1$auc[i], lm.sd = sd(res[, "lm"], na.rm = TRUE),
    lm.failures = sum(is.na(res[, "lm"])),
    empirical.bias = mean(res[, "empirical"]) - table1$auc[i]
  )
  print(rows[[i]], digits = 3, row.names = FALSE)
}
rows <- do.call(rbind, rows)
cat(sprintf("Combined z of the means against the published ones: %.1f\n\n", sum(rows$z) / sqrt(nrow(rows))))

## 3. Five-category data: M&P Tables 1-3 (pp. 22-23) ------------------------
# M&P p. 20: conventional binormal a = 1.5, b = 0.8, 30 negative and 40
# positive cases, 5 categories with the population operating points
# (FPF, TPF) = (0.00025, 0.10), (0.000328, 0.25), (0.03040, 0.50),
# (0.28114, 0.85). FPF 0.000328 does not lie on the curve (TPF = 0.25 needs
# FPF = 0.00328): the cut-points are taken from the TPFs.
cat("== 3. M&P Tables 1-3, 5 categories, 30 + 40 cases\n")
a <- 1.5
b <- 0.8
cuts <- sort((a - qnorm(c(0.10, 0.25, 0.50, 0.85))) / b)
fpf.at <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
true.tpf <- pnorm(a + b * qnorm(fpf.at))
true.area <- pnorm(a / sqrt(1 + b^2))
res <- do.call(rbind, mclapply(seq_len(B), function(j) {
  r <- roc(controls = findInterval(rnorm(30), cuts), cases = findInterval(rnorm(40, a / b, 1 / b), cuts), direction = "<", quiet = TRUE)
  lm.area <- tryCatch(as.numeric(auc(smooth(r))), error = function(e) NA_real_)
  s <- tryCatch(smooth(r, fit = "proper"), error = function(e) NULL)
  if (is.null(s)) {
    return(c(rep(NA_real_, 8), lm.area))
  }
  c(approx(1 - s$specificities, s$sensitivities, xout = fpf.at, ties = max)$y, as.numeric(auc(s)), lm.area)
}, mc.cores = cores))
ok <- !is.na(res[, 8])
print(data.frame(
  FPF = fpf.at, true.TPF = true.tpf, mean.TPF = colMeans(res[ok, 1:7]),
  RMSE = sqrt(colMeans(sweep(res[ok, 1:7], 2, true.tpf)^2)),
  published.mean = c(0.615, 0.705, 0.802, 0.857, 0.921, 0.960, 0.988),
  published.RMSE = c(0.1557, 0.1150, 0.0795, 0.0639, 0.0459, 0.0297, 0.0111)
), digits = 3, row.names = FALSE)
cat(sprintf("Area: mean %.4f, RMSE %.4f (published 0.879, 0.0437; true %.4f); failures %d of %d\n",
  mean(res[ok, 8]), sqrt(mean((res[ok, 8] - true.area)^2)), true.area, sum(!ok), B))
cat(sprintf("Least-squares fit: area %.4f, RMSE %.4f; failures %d of %d\n\n",
  mean(res[, 9], na.rm = TRUE), sqrt(mean((res[, 9] - true.area)^2, na.rm = TRUE)), sum(is.na(res[, 9])), B))

## 4. Failures: how often do the data imply perfect performance? ------------
# M&P p. 27: no curve-parameter estimates are possible then (an error in
# pROC). Probability that all the controls are below all the cases, by
# simulation (no fitting).
cat("== 4. Probability of perfect separation (2e5 samples)\n")
N <- 2e5
d <- table1[table1$auc == 0.95 & table1$shape == "high", ]
b2 <- (1 + d$c) / (1 - d$c)
a2 <- d$da * sqrt(1 + d$c^2) / (1 - d$c)
y <- function(x) -((b2 * x - a2)^2 - x^2) / 2
p.high <- mean(replicate(N, max(y(rnorm(50))) < min(y(rnorm(50, a2 / b2, 1 / b2)))))
p.cat <- mean(replicate(N, max(findInterval(rnorm(30), cuts)) < min(findInterval(rnorm(40, a / b, 1 / b), cuts))))
cat(sprintf("AUC 0.95 high skew, 50 + 50: %.5f (expected per %d datasets: %.1f)\n", p.high, B, B * p.high))
cat(sprintf("5 categories, 30 + 40: %.5f (expected per %d datasets: %.1f)\n", p.cat, B, B * p.cat))
