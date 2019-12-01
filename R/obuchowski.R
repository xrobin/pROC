# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2011-2014 Xavier Robin, Alexandre Hainard, Natacha Turck,
# Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez
# and Markus Müller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Formula 3 from Obuchowski 2004, p. 1123
# Variance of an AUC given kappa
var.theta.obuchowski <- function(theta, kappa) {
    A <- qnorm(theta) * 1.414
    (0.0099 * exp(-A^2/2)) * ((5 * A^2 + 8) + (A^2 + 8)/kappa)
}

# Formulas from Obuchowski 1997, table 1 p. 1531
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
  (qnorm(FPRi) + A * B * (1 + B^2)^(-1) ) * (1 + B^2)^(1/2)
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
  binormal <- smooth(roc, method="binormal")$model
  A <- unname(coefficients(binormal)[1])
  B <- unname(coefficients(binormal)[2])
  kappa <- length(roc$controls) / length(roc$cases)

  if (!identical(attr(roc$auc, "partial.auc"), FALSE)) {
    FPR1 <- 1 - attr(roc$auc, "partial.auc")[2]
    FPR2 <- 1 - attr(roc$auc, "partial.auc")[1]
    va <- var.params.obuchowski(A, B, kappa, FPR1, FPR2)
  }
  else {
    va <- var.params.obuchowski(A, B, kappa)
  }
  return(va)
}

# Variance of a ROC curve given the parameters
# Obuchowski 1997, formula 4 p. 1530
# A and B: params of the binormal ROC curve
# kappa: proportion controls / cases
# FPR1, FPR2: the bottom (1) or top (2) bounds of the pAUC interval
var.params.obuchowski <- function(A, B, kappa, FPR1, FPR2) {
  if (!missing(FPR1) && !is.null(FPR1) && !missing(FPR1) && !is.null(FPR2)) {
    f.partial(A, B, FPR1, FPR2)^2 * (1 + B^2 / kappa + A^2/2) + g.partial(A, B, FPR1, FPR2)^2 * B^2 * (1 + kappa) / (2*kappa)
  }
  else {
    f.full(A, B)^2 * (1 + B^2 / kappa + A^2/2) + g.full(A, B)^2 * B^2 * (1 + kappa) / (2*kappa)
  }
}

# Covariance of 2 given 'roc' objects (under the alternative hypothesis)
cov.roc.obuchowski <- function(roc1, roc2) {
  binormal1 <- smooth(roc1, method="binormal")$model
  A1 <- unname(coefficients(binormal1)[1])
  B1 <- unname(coefficients(binormal1)[2])
  binormal2 <- smooth(roc2, method="binormal")$model
  A2 <-unname(coefficients(binormal2)[1])
  B2 <- unname(coefficients(binormal2)[2])
  kappa <- length(roc1$controls) / length(roc1$cases)
  ra <- cor(roc1$cases, roc2$cases)
  rn <- cor(roc1$controls, roc2$controls)
  if (!identical(attr(roc1$auc, "partial.auc"), FALSE)) {
    FPR11 <- 1 - attr(roc1$auc, "partial.auc")[2]
    FPR12 <- 1 - attr(roc1$auc, "partial.auc")[1]
    FPR21 <- 1 - attr(roc2$auc, "partial.auc")[2]
    FPR22 <- 1 - attr(roc2$auc, "partial.auc")[1]
    co <- cov.params.obuchowski(A1, B1, A2, B2, rn, ra, kappa, FPR11, FPR12, FPR21, FPR22)
  }
  else {
    co <- cov.params.obuchowski(A1, B1, A2, B2, rn, ra, kappa)
  }
  return(co)
}

# Covariance under the null hypothesis
# roc1 is taken as null
cov0.roc.obuchowski <- function(roc1, roc2) {
  binormal <- smooth(roc, method="binormal")$model
  A <- unname(coefficients(binormal)[1])
  B <- unname(coefficients(binormal)[2])
  R <- length(roc1$controls) / length(roc1$cases)
  ra <- cor(roc1$cases, roc2$cases)
  rn <- cor(roc1$controls, roc2$controls)
  if (!identical(attr(roc1$auc, "partial.auc"), FALSE)) {
    FPR1 <- attr(roc1$auc, "partial.auc")[2]
    FPR2 <- attr(roc1$auc, "partial.auc")[1]
    co <- cov.params.obuchowski(A, B, A, B, rn, ra, kappa, FPR1, FPR2, FPR1, FPR2)
  }
  else {
    co <- cov.params.obuchowski(A, B, A, B, rn, ra, kappa)
  }
  return(co)
}


# Covariance of a ROC curve given the parameters
# Obuchowski 1997, formula 5 p. 1531
# (A|B)(1|2): A and B params of the binormal ROC curve
# rn, ra: correlation of the results in ROC curves 1 and 2 in controls (n) and cases (a) patients
# kappa: proportion controls / cases
# FPR(1|2)(1|2): the bounds of the pAUC interval:
#    ***** ROC curve 1 or 2
#         ***** bottom (1) or top (2) of the interval
cov.params.obuchowski <- function(A1, B1, A2, B2, rn, ra, kappa, FPR11, FPR12, FPR21, FPR22) {
  if (!missing(FPR11) && !is.null(FPR11) &&
      !missing(FPR12) && !is.null(FPR12) &&
      !missing(FPR21) && !is.null(FPR21) &&
      !missing(FPR22) && !is.null(FPR22)) {
    f1 <- f.partial(A1, B1, FPR11, FPR12)
    f2 <- f.partial(A2, B2, FPR21, FPR22)
    g1 <- g.partial(A1, B1, FPR11, FPR12)
    g2 <- g.partial(A2, B2, FPR21, FPR22)
  }
  else {
    f1 <- f.full(A1, B1)
    f2 <- f.full(A2, B2)
    g1 <- g.full(A1, B1)
    g2 <- g.full(A2, B2)
  }
  f1 * f2 * (ra + rn * B1 * B2 / kappa + ra^2 * A1 * A2  / 2) +
    g1 * g2 * (B1 * B2 * (rn^2 + kappa * ra^2) / (2 *  kappa)) + 
    f1 * g2 * (ra^2 * A1 * B2 / 2) + f2 * g1 * (ra^2 * A2 * B1 / 2)
}

# Variance of a difference between two ROC curves given the parameters
# Obuchowski 1997, formula 4 and 5 p. 1530--1531
# (A|B)(1|2): A and B params of the binormal ROC curve
# rn, ra: correlation of the results in ROC curves 1 and 2 in controls (n) and cases (a) patients
# kappa: proportion controls / cases
# FPR(1|2)(1|2): the bounds of the pAUC interval:
#    ***** ROC curve 1 or 2
#         ***** bottom (1) or top (2) of the interval
vardiff.params.obuchowski <- function(A1, B1, A2, B2, rn, ra, kappa, FPR11, FPR12, FPR21, FPR22) {
  var(A1, B1, kappa, FPR11, FPR12) + var(A2, B2, kappa, FPR21, FPR22) +
    2 * cov(A1, B1, A2, B2, rn, ra, kappa, FPR11, FPR12, FPR21, FPR22)
}

# Variance of a difference between two ROC curves given the parameters
# under the null hypothesis. ROC curve 1 is taken as null
# Obuchowski 1997, formula 4 and 5 p. 1530--1531
# (A|B)(1|2): A and B params of the binormal ROC curve
# rn, ra: correlation of the results in ROC curves 1 and 2 in controls (n) and cases (a) patients
# kappa: proportion controls / cases
# FPR(1|2)(1|2): the bounds of the pAUC interval:
#    ***** ROC curve 1 or 2
#         ***** bottom (1) or top (2) of the interval
vardiff0.params.obuchowski <- function(A1, B1, A2, B2, rn, ra, kappa, FPR11, FPR12, FPR21, FPR22) {
  2 * var(A1, B1, kappa, FPR11, FPR12) + 
    2 * cov(A1, B1, A2, B2, rn, ra, kappa, FPR11, FPR12, FPR21, FPR22)
}
