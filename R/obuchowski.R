# Formula 3 from Obuchowski 2004, p. 1123
# Variance of an AUC given kappa
var.theta.obuchowski <- function(theta, kappa) {
    A <- qnorm(theta) * 1.414
    (0.0099 * exp(-A^2/2)) * ((5 * A^2 + 8) + (A^2 + 8)/kappa)
}

# Formulas from Obuchowski 1997, p. 1530-1531
expr1 <- function(A, B) {
  exp(-A^2/(2 * (1 + B^2)))
}
expr2 <- function(B) {
  1 + B^2
}
expr3 <- function(cdagger1, cdagger2) {
  pnorm(cdagger1) - pnorm(cdagger2)
}
expr4 <- function(cddagger1, cddagger2) {
 exp(-cddagger1) - exp(- cddagger2)
}

cdagger <- function (A, B, FPRi) {
  (qnorm(FPRi) + A * B * (1 + B^2)^(-1) )* (1 + B^2)^(1/2)
}

cddagger <- function(cdagger) {
  cdagger^2 / 2
}

f.full <- function(A, B) {
  expr1 <- expr1(A, B)
  expr2 <- expr2(B)
  expr1 * (2 * pi * expr2) ^ (-1/2)
}

f.partial <- function(A, B, FPR1, FPR2) {
  cdagger1 <- cdagger(A, B, FPR1)
  cdagger2 <- cdagger(A, B, FPR2)
  expr1 <- expr1(A, B)
  expr2 <- expr2(B)
  expr3 <- expr3(cdagger1, cdagger2)
  expr1 * (2 * pi * expr2) ^ (-1/2) * expr3
}

g.full <- function(A, B) {
  expr1 <- expr1(A, B)
  expr2 <- expr2(B)
  - expr1 * A * B * (2 * pi * expr2^3) ^ (-1/2)
}

g.partial <- function(A, B, FPR1, FPR2) {
  cdagger1 <- cdagger(A, B, FPR1)
  cdagger2 <- cdagger(A, B, FPR2)
  cddagger1 <- cddagger(cdagger1)
  cddagger2 <- cddagger(cdagger2)
  expr1 <- expr1(A, B)
  expr2 <- expr2(B)
  expr3 <- expr3(cdagger1, cdagger2)
  expr4 <- expr4(cddagger1, cddagger2)
  # WARNING: we have set (-expr4), in contradiction with Obuchowski paper
  expr1 * (2 * pi * expr2) ^ (-1) * (-expr4) - A * B * expr1 * (2 * pi * expr2^3) ^ (-1/2) * expr3
}

# Variance of a ROC curve given a 'roc' object
var.roc.obuchowski <- function(roc) {
  A <- (mean(roc$cases) - mean(roc$controls)) / sd(roc$cases)
  B <- sd(roc$controls) / sd(roc$cases)
  kappa <- length(roc$controls) / length(roc$cases)

  if (!identical(attr(roc$auc, "partial.auc"), FALSE)) {
    FPR1 <- attr(roc$auc, "partial.auc")[2]
    FPR2 <- attr(roc$auc, "partial.auc")[1]
    var.params.obuchowski(A, B, kappa, FPR1, FPR2)
  }
  else {
    var.params.obuchowski(A, B, kappa)
  }
}

var.params.obuchowski <- function(A, B, kappa, FPR1, FPR2) {
  if (!missing(FPR1) && !is.null(FPR1) && !missing(FPR1) && !is.null(FPR2)) {
    f.partial(A, B, FPR1, FPR2)^2 * (1 + B^2 / kappa + A^2/2) + g.partial(A, B, FPR1, FPR2)^2 * B^2 * (1 + kappa) / (2*kappa)
  }
  else {
    f.full(A, B)^2 * (1 + B^2 / kappa + A^2/2) + g.full(A, B)^2 * B^2 * (1 + kappa) / (2*kappa)
  }
}

# Covariance under the alternative hypothesis
cov.roc.obuchowski <- function(roc1, roc2) {
  A1 <- (mean(roc1$cases) - mean(roc1$controls)) / sd(roc1$cases)
  B1 <- sd(roc1$controls) / sd(roc1$cases)
  A2 <- (mean(roc2$cases) - mean(roc2$controls)) / sd(roc2$cases)
  B2 <- sd(roc2$controls) / sd(roc2$cases)
  R <- length(roc1$controls) / length(roc1$cases)
  if (!identical(attr(roc1$auc, "partial.auc"), FALSE)) {
    FPR11 <- attr(roc1$auc, "partial.auc")[2]
    FPR21 <- attr(roc1$auc, "partial.auc")[1]
    FPR12 <- attr(roc2$auc, "partial.auc")[2]
    FPR22 <- attr(roc2$auc, "partial.auc")[1]
    f1 <- f.partial(A1, B1, FPR11, FPR21)
    f2 <- f.partial(A2, B2, FPR12, FPR22)
    g1 <- g.partial(A1, B1, FPR11, FPR21)
    g2 <- g.partial(A2, B2, FPR12, FPR22)
  }
  else {
    f1 <- f.full(A1, B1)
    f2 <- f.full(A2, B2)
    g1 <- g.full(A1, B1)
    g2 <- g.full(A2, B2)
  }
  ra <- cor(roc1$cases, roc2$cases)
  rn <- cor(roc1$controls, roc2$controls)
  co <- f1 * f2 * (ra + rn * B1 * B2 / R + ra^2 * A1 * A2  / 2) +
        g1 * g2 * (B1 * B2 * (rn^2 + R * ra^2) / (2 * R)) + 
        f1 * g2 * (ra^2 * A1 * B2 / 2) + f2 * g1 * (ra^2 * A2 * B1 / 2)
  return(co)
}

# Covariance under the null hypothesis
cov0.roc.obuchowski <- function(roc1, roc2) {
  A <- (mean(roc1$cases) - mean(roc1$controls)) / sd(roc1$cases)
  B <- sd(roc1$controls) / sd(roc1$cases)
  R <- length(roc1$controls) / length(roc1$cases)
  if (!identical(attr(roc1$auc, "partial.auc"), FALSE)) {
    FPR1 <- attr(roc1$auc, "partial.auc")[2]
    FPR2 <- attr(roc1$auc, "partial.auc")[1]
    f <- f.partial(A, B, FPR1, FPR2)
    g <- g.partial(A, B, FPR1, FPR2)
  }
  else {
    f <- f.full(A, B)
    g <- g.full(A, B)
  }
  ra <- cor(roc1$cases, roc2$cases)
  rn <- cor(roc1$controls, roc2$controls)
  co <- f * f * (ra + rn * B * B / R + ra^2 * A * A  / 2) +
        g * g * (B * B * (rn^2 + R * ra^2) / (2 * R)) + 
        f * g * (ra^2 * A * B / 2) + f * g * (ra^2 * A * B / 2)
  return(co)
}

# Variance of a delta under the alternative hypothesis
var.delta.obuchowski <- function(roc1, roc2) {
  var.roc.obuchowski(roc1) + var.roc.obuchowski(roc2) - 2 * cov.roc.obuchowski(roc1, roc2)
}

# Variance of a delta under the null hypothesis
var0.delta.obuchowski <- function(roc1, roc2) {
  if (roc1$auc < roc2$auc) {
    roc.min <- roc1
    roc.max <- roc2
  }
  else {
    roc.min <- roc2
    roc.max <- roc1
  }
  2 * var.roc.obuchowski(roc.min) - 2 * cov0.roc.obuchowski(roc.min, roc.max)
}

