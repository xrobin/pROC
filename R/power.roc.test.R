# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2011 Xavier Robin, Alexandre Hainard, Natacha Turck,
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

power.roc.test <- function(...)
  UseMethod("power.roc.test")

power.roc.test.roc <- function(roc1, roc2, sig.level = 0.05, power = NULL, alternative = c("two.sided", "one.sided"), reuse.auc=TRUE, ...) {
  # Basic sanity checks
  if (!is.null(power) && (power < 0 || power > 1))
    stop("'power' must range from 0 to 1")
  if (!is.null(sig.level) && (sig.level < 0 || sig.level > 1))
    stop("'sig.level' must range from 0 to 1")
  
  # check that the AUC of roc1 was computed, or do it now
  if (is.null(roc1$auc) | !reuse.auc) {
    roc1$auc <- auc(roc1, ...)
  }
  if (!is.null(attr(roc1$auc, "partial.auc.correct")) && attr(roc1$auc, "partial.auc.correct")) {
    stop("Cannot compute power with corrected partial AUCs")
  }
  if (attr(roc1$auc, "percent")) {
    roc1$auc <- roc1$auc / 100
    attr(roc1$auc, "partial.auc") <- attr(roc1$auc, "partial.auc") / 100
    attr(roc1$auc, "percent") <- FALSE
  }
    

  if (!missing(roc2) && !is.null(roc2)) {
    alternative <- match.arg(alternative)
    if (!is.null(sig.level) && alternative == "two.sided") {
      sig.level <- sig.level / 2
    }
    
    if ("roc" %in% class(roc2)) {
      # check that the AUC of roc2 was computed, or do it now
      if (is.null(roc2$auc) | !reuse.auc) {
        roc2$auc <- auc(roc2, ...)
      }
      if (!is.null(attr(roc2$auc, "partial.auc.correct")) && attr(roc2$auc, "partial.auc.correct")) {
        stop("Cannot compute power with corrected partial AUCs")
      }
      if (attr(roc2$auc, "percent")) {
        roc2$auc <- roc2$auc / 100
        attr(roc2$auc, "partial.auc") <- attr(roc2$auc, "partial.auc") / 100
        attr(roc2$auc, "percent") <- FALSE
      }

      # Make sure the ROC curves are paired
      rocs.are.paired <- are.paired(roc1, roc2)
      if (!rocs.are.paired) {
        stop("The sample size for a difference in AUC cannot be applied to unpaired ROC curves yet.")
      }
      # Make sure the AUC specifications are identical
      attr1 <- attributes(roc1$auc); attr1$roc <- NULL
      attr2 <- attributes(roc2$auc); attr2$roc <- NULL
      if (!identical(attr1, attr2)) {
        stop("Different AUC specifications in the ROC curves.")
      }

      if (!is.null(attr(roc1$auc, "partial.auc.focus")) && attr(roc1$auc, "partial.auc.focus") != "specificity") {
        stop("Power tests for pAUC in sensitivity are not implemented yet!")
      }
    
      # check that the same region was requested in auc. Otherwise, issue a warning
      if (!identical(attributes(roc1$auc)[names(attributes(roc1$auc))!="roc"], attributes(roc2$auc)[names(attributes(roc2$auc))!="roc"]))
        warning("Different AUC specifications in the ROC curves. Enforcing the inconsistency, but unexpected results may be produced.")
 
      ncontrols <- length(roc1$controls)
      ncases <- length(roc1$cases)
      kappa <- ncontrols / ncases

      # Power test
      if (is.null(power)) {
        if (is.null(sig.level))
          stop("'sig.level' or 'power' must be provided.")
        zalpha <- qnorm(sig.level)
        zbeta <- zbeta.obuchowski(roc1, roc2, zalpha)
        power <- 1 - pnorm(zbeta)
      }
      # sig.level
      else if (is.null(sig.level)) {
        zbeta <- qnorm(1 - power)
        zalpha <- zalpha.obuchowski(roc1, roc2, zbeta)
        sig.level <- pnorm(zalpha)
      }
      # Sample size
      else {
        zalpha <- qnorm(sig.level)
        zbeta <- qnorm(1 - power)

        ncases <- ncases.obuchowski(roc1, roc2, zalpha, zbeta)
        ncontrols <- kappa * ncases
      }

      # Restore sig.level if two.sided
      if (alternative == "two.sided") {
        sig.level <- sig.level * 2
      }
      return(structure(list(ncases=ncases, ncontrols=ncontrols, auc1=roc1$auc, auc2=roc2$auc, sig.level=sig.level, power=power, alternative=alternative, method="Two ROC curves power calculation"), class="power.htest"))
    }
    else {
      stop("'roc2' must be an object of class 'roc'.")
    }
  }
  else {
    if (is.null(sig.level) || is.null(power)) {
      ncontrols <- length(roc1$controls)
      ncases <- length(roc1$cases)
    }
    auc <- roc1$auc
    return(power.roc.test(ncontrols = ncontrols, ncases = ncases, auc = auc, sig.level = sig.level, power = power, alternative = alternative, ...))
  }
}

power.roc.test.numeric <- function(ncontrols = NULL, ncases = NULL, auc = NULL, sig.level = 0.05, power = NULL,  kappa = 1, alternative = c("two.sided", "one.sided"), ...) {
  # basic sanity checks
  if (!is.null(ncases) && ncases < 0)
    stop("'ncases' must be positive")
  if (!is.null(ncontrols) && ncontrols < 0)
    stop("'ncontrols' must be positive")
  if (!is.null(kappa) && kappa < 0)
    stop("'kappa' must be positive")
  if (!is.null(power) && (power < 0 || power > 1))
    stop("'power' must range from 0 to 1")
  if (!is.null(sig.level) && (sig.level < 0 || sig.level > 1))
    stop("'sig.level' must range from 0 to 1")
  
  alternative <- match.arg(alternative)
  if (alternative == "two.sided" && !is.null(sig.level)) {
    sig.level <- sig.level / 2
  }

  # determine AUC
  if (is.null(auc)) {
    if (is.null(ncontrols) || is.null(ncases))
      stop("'ncontrols' and 'ncases' or 'auc' must be provided.")
    else if (is.null(power))
      stop("'power' or 'auc' must be provided.")
    else if (is.null(sig.level))
      stop("'sig.level' or 'auc' must be provided.")
    kappa <- ncontrols / ncases
    zalpha <- qnorm(sig.level)
    zbeta <- qnorm(1 - power)

    tryCatch(
             root <- uniroot(power.roc.test.optimize.auc.function, interval=c(0.5, 1-1e-16), ncontrols=ncontrols, ncases=ncases, zalpha=zalpha, zbeta=zbeta),
             error=function(e) {stop(sprintf("AUC could not be solved:\n%s", e))}
             )
    auc <- root$root
  }

  # Determine number of patients (sample size)
  else if (is.null(ncases) && is.null(ncontrols)) {
    if (is.null(power))
      stop("'power' or 'ncases' and 'ncontrols' must be provided.")
    else if (is.null(kappa))
      stop("'kappa' must be provided.")
    else if (is.null(auc))
      stop("'auc' or 'ncases' and 'ncontrols' must be provided.")
    else if (is.null(sig.level))
      stop("'sig.level' or 'ncases' and 'ncontrols' must be provided.")

    theta <- as.numeric(auc)
    Vtheta <- var.theta.obuchowski(theta, kappa)
    zalpha <- qnorm(sig.level)
    zbeta <- qnorm(1 - power)
    ncases <- (zalpha * sqrt(0.0792 * (1 + 1/kappa)) + zbeta * sqrt(Vtheta))^2 / (theta - 0.5)^2
    ncontrols <- kappa * ncases
  }
  
  # Determine power
  else if (is.null(power)) { 
    if (is.null(ncontrols) || is.null(ncases))
      stop("'ncontrols' and 'ncases' or 'power' must be provided.")
    else if (is.null(auc))
      stop("'auc' or 'power' must be provided.")
    else if (is.null(sig.level))
      stop("'sig.level' or 'power' must be provided.")
    kappa <- ncontrols / ncases

    theta <- as.numeric(auc)
    Vtheta <- var.theta.obuchowski(theta, kappa)
    zalpha <- qnorm(sig.level)

    rs.beta <- Vtheta
    rs.alpha <- 0.0792 * (1 + 1 / kappa)
    zar <- zalpha * sqrt(rs.alpha)
    a <- rs.beta
    b <- 2 * zar * sqrt(rs.beta)
    c <- (zar^2 - ncases * (theta^2 - theta + 0.25))
    zbeta <- solve.2deg.eqn(a, b, c)
    power <- 1 - pnorm(zbeta)
  }

  # Determine sig.level
  else  if (is.null(sig.level)) { 
    if (is.null(ncontrols) || is.null(ncases))
      stop("'ncontrols' and 'ncases' or 'sig.level' must be provided.")
    else if (is.null(auc))
      stop("'auc' or 'sig.level' must be provided.")
    else if (is.null(power))
      stop("'power' or 'sig.level' must be provided.")
    kappa <- ncontrols / ncases

    theta <- as.numeric(auc)
    Vtheta <- var.theta.obuchowski(theta, kappa)
    
    zbeta <- qnorm(1 - power)

    zbr <- qnorm(1 - power) * sqrt(Vtheta)
    rs.alpha <- 0.0792 * (1 + 1 / kappa)
    a <- rs.alpha
    b <- 2 * zbr * sqrt(rs.alpha)
    c <- (zbr^2 - ncases * (theta^2 - theta + 0.25))
    zalpha <- solve.2deg.eqn(a, b, c)
    sig.level <- pnorm(zalpha)
  }
  else {
    stop("One of 'power', 'sig.level', 'auc', or both 'ncases' and 'ncontrols' must be NULL.")
  }
  # Restore sig.level if two.sided
  if (alternative == "two.sided") {
    sig.level <- sig.level * 2
  }
  return(structure(list(ncases=ncases, ncontrols=ncontrols, auc=auc, sig.level=sig.level, power=power, method="One ROC curve power calculation"), class="power.htest"))
}

# Formula 3 from Obuchowski 2004, p. 1123
var.theta.obuchowski <- function(theta, kappa) {
    A <- qnorm(theta) * 1.414
    (0.0099 * exp(-A^2/2)) * ((5 * A^2 + 8) + (A^2 + 8)/kappa)
}

power.roc.test.optimize.auc.function <- function(x, ncontrols, ncases, zalpha, zbeta) {
  kappa <- ncontrols / ncases
  Vtheta <- var.theta.obuchowski(x, kappa)
  (zalpha * sqrt(0.0792 * (1 + 1/kappa)) + zbeta * sqrt(Vtheta))^2 / (x - 0.5)^2 - ncases
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
  expr1 * (2 * pi * expr2) ^ (-1) * expr4 - A * B * expr1 * (2 * pi * expr2^3) ^ (-1/2) * expr3
}

var.roc.obuchowski <- function(roc) {
  A <- (mean(roc$cases) - mean(roc$controls)) / sd(roc$cases)
  B <- sd(roc$controls) / sd(roc$cases)
  R <- length(roc$controls) / length(roc$cases)

  if (!identical(attr(roc$auc, "partial.auc"), FALSE)) {
    FPR1 <- attr(roc$auc, "partial.auc")[2]
    FPR2 <- attr(roc$auc, "partial.auc")[1]
    f.partial(A, B, FPR1, FPR2)^2 * (1 + B^2 / R + A^2/2) + g.partial(A, B, FPR1, FPR2)^2 * B^2 * (1 + R) / (2*R)
  }
  else {
    f.full(A, B)^2 * (1 + B^2 / R + A^2/2) + g.full(A, B)^2 * B^2 * (1 + R) / (2*R)    
  }
}

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

var.delta.obuchowski <- function(roc1, roc2) {
  var.roc.obuchowski(roc1) + var.roc.obuchowski(roc2) - 2 * cov.roc.obuchowski(roc1, roc2)
}

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

ncases.obuchowski <- function(roc1, roc2, zalpha, zbeta) {
  delta <- roc1$auc - roc2$auc
  na <- (zalpha * sqrt(var0.delta.obuchowski(roc1, roc2)) +
       zbeta * sqrt(var.delta.obuchowski(roc1, roc2))) ^2 /
       delta^2
  return(na)
}

zalpha.obuchowski <- function(roc1, roc2, zbeta) {
  delta <- roc1$auc - roc2$auc
  ncases <- length(roc1$cases)
  v0 <- var0.delta.obuchowski(roc1, roc2)
  va <- var.delta.obuchowski(roc1, roc2)
  a <- v0
  b <- 2 * zbeta * sqrt(v0) * sqrt(va)
  c <- zbeta^2 * va - ncases * delta ^ 2
  return(solve.2deg.eqn(a, b, c))
}

zbeta.obuchowski <- function(roc1, roc2, zalpha) {
  delta <- roc1$auc - roc2$auc
  ncases <- length(roc1$cases)
  v0 <- var0.delta.obuchowski(roc1, roc2)
  va <- var.delta.obuchowski(roc1, roc2)
  a <- va
  b <- 2 * zalpha * sqrt(va) * sqrt(v0)
  c <- zalpha^2 * v0 - ncases * delta ^ 2
  return(solve.2deg.eqn(a, b, c))
}

solve.2deg.eqn <- function(a, b, c) {
  return((- b - sqrt(b^2 - 4*a*c)) / (2*a))
}
