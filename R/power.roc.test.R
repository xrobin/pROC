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

power.roc.test.roc <- function(roc1, roc2, sig.level = 0.05, power = NULL, alternative = c("two.sided", "one.sided"), ...) {
  # Basic sanity checks
  if (!is.null(power) && (power < 0 || power > 1))
    stop("'power' must range from 0 to 1")
  if (!is.null(sig.level) && (sig.level < 0 || sig.level > 1))
    stop("'sig.level' must range from 0 to 1")
  
  if (!is.null(roc2)) {
    alternative <- match.arg(alternative)
    if (!is.null(sig.level) && alternative == "two.sided") {
      sig.level <- sig.level / 2
    }
    
    if (class(roc2) == "roc") {
      rocs.are.paired <- are.paired(roc1, roc2)
      if (!rocs.are.paired) {
        stop("The sample size for a difference in AUC cannot be applied to unpaired ROC curves yet.")
      }

      auc1 <- as.numeric(auc(roc1))
      auc2 <- as.numeric(auc(roc2))
      ncontrols <- length(roc1$controls)
      ncases <- length(roc1$cases)
      kappa <- ncontrols / ncases
      Vtheta1 <- var.theta.obuchowski(auc1, kappa)
      Vtheta2 <- var.theta.obuchowski(auc1, kappa)
      cov12 <- cov(roc1, roc2)
      vardiff <- Vtheta1 + Vtheta2 + 2 * cov12
      # Variance under null hypothesis
      vardiff0 <- 2 * var.theta.obuchowski(.5, kappa) + 2 * cov12
      Dauc.sq <- (auc1 - auc2)^2
      

      # Power test
      if (is.null(power)) {
        if (is.null(sig.level))
          stop("'sig.level' or 'power' must be provided.")
        zalpha <- qnorm(sig.level)

        a <- vardiff
        b <- 2 * zalpha * sqrt(vardiff) * sqrt(vardiff0)
        c <- zalpha^2 * sqrt(vardiff0) - ncases * Dauc.sq
        zbeta <- (- b - sqrt(b^2 - 4*a*c)) / (2*a)
        power <- 1 - pnorm(zbeta)
      }
      # sig.level
      else if (is.null(sig.level)) {
        zbeta <- qnorm(1 - power)

        a <- vardiff0
        b <- 2 * zbeta * sqrt(vardiff) * sqrt(vardiff0)
        c <- zbeta^2 * sqrt(vardiff0) - ncases * Dauc.sq
        zalpha <- (- b - sqrt(b^2 - 4*a*c)) / (2*a)
        sig.level <- pnorm(zalpha)
      }
      # Sample size
      else {
        zalpha <- qnorm(sig.level)
        zbeta <- qnorm(1 - power)

        ncases <- ((zalpha * sqrt(vardiff0) + zbeta * sqrt(vardiff)) ^ 2) / ((auc1 - auc2) ^ 2)
        ncontrols <- kappa * ncases
      }

      # Restore sig.level if two.sided
      if (alternative == "two.sided") {
        sig.level <- sig.level * 2
      }
      return(structure(list(ncases=ncases, ncontrols=ncontrols, auc1=auc1, auc2=auc2, sig.level=sig.level, power=power, alternative=alternative, method="Two ROC curves power calculation"), class="power.htest"))
    }
    else {
      stop("'roc2' must be an object of class 'roc'.")
    }
  }
  else {
    ncontrols <- length(roc$controls)
    ncases <- length(roc$cases)
    auc <- auc.roc(roc)
    power.roc.test(ncontrols = ncontrols, ncases = ncases, auc = auc, sig.level = sig.level, power = power, alternative = alternative, ...)
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
  if (alternative == "two.sided") {
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

    power.roc.test.optimize.auc.function <- function(x, ncontrols, ncases, sig.level, power) {
      kappa <- ncontrols / ncases
      Vtheta <- var.theta.obuchowski(theta, kappa)
      zalpha <- qnorm(sig.level)
      zbeta <- qnorm(1 - power)
      (zalpha * sqrt(0.0792 * (1 + 1/kappa)) + zbeta * sqrt(Vtheta))^2 / (x - 0.5)^2 - ncases
    }
    tryCatch(root <- uniroot(power.roc.test.optimize.auc.function, interval=c(0.5, 1-1e-16), ncontrols=ncontrols, ncases=ncases, sig.level=sig.level, power=power), error=function(e) stop(sprintf("AUC could not be solved:\n%s", e)))
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
    zbeta <- (- b - sqrt(b^2 - 4*a*c)) / (2*a)
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
    zalpha <- (- b - sqrt(b^2 - 4*a*c)) / (2*a)
    sig.level <- pnorm(zalpha)
  }
  else {
    stop("One of 'power', 'sig.level', 'auc', or both 'ncases' and 'ncontrols' must be NULL.")
  }
  # Restore sig.level if two.sided
  if (alternative == "two.sided") {
    sig.level <- sig.level * 2
  }
  structure(list(ncases=ncases, ncontrols=ncontrols, auc=auc, sig.level=sig.level, power=power, method="One ROC curve power calculation"), class="power.htest")
}

# Formula 3 in Obuchowski, p. 1123
var.theta.obuchowski <- function(theta, kappa) {
    A <- qnorm(theta) * 1.414
    (0.0099 * exp(-A^2/2)) * ((5 * A^2 + 8) + (A^2 + 8)/kappa)
}


# Formulas from table 1, Obuchowski 1997, p. 1531
compute.expr1 <- function(A, B) {
  exp(-A^2/(2 * (1 + B^2)))
}
compute.expr2 <- function(B) {
  1 + B^2
}
compute.expr3 <- function(cdagger1, cdagger2) {
  pnorm(cdagger1) - pnorm(cdagger2)
}
compute.expr4 <- function(cddagger1, cddagger2) {
 exp(-cddagger1) - exp(- cddagger2)
}

cdagger <- function (A, B, FPRi) {
  (qnorm(FPRi) + A * B * (1 + B^2)^(-1) )* (1 + B^2)^(1/2)
}

cddagger <- function(cdagger) {
  cdagger^2 / 2
}

f.full <- function(A, B) {
  expr1 <- compute.expr1(A, B)
  expr2 <- compute.expr2(B)
  expr1 * (2 * pi * expr2) ^ (-1/2)
}

f.partial <- function(A, B, FPR1, FPR2) {
  cdagger1 <- cdagger(A, B, FPR1)
  cdagger2 <- cdagger(A, B, FPR2)
  expr1 <- compute.expr1(A, B)
  expr2 <- compute.expr2(B)
  expr3 <- compute.expr3(cdagger1, cdagger2)
  expr1 * (2 * pi * expr2) ^ (-1/2) * expr3
}
g.full <- function(A, B) {
  expr1 <- compute.expr1(A, B)
  expr2 <- compute.expr2(B)
  - expr1 * A * B * (2 * pi * expr2^3) ^ (-1/2)
}
g.partial <- function(A, B, FPR1, FPR2) {
  cdagger1 <- cdagger(A, B, FPR1)
  cdagger2 <- cdagger(A, B, FPR2)
  cddagger1 <- cddagger(cdagger1)
  cddagger2 <- cddagger(cdagger2)
  expr1 <- compute.expr1(A, B)
  expr2 <- compute.expr2(B)
  expr3 <- compute.expr3(cdagger1, cdagger2)
  expr4 <- compute.expr4(cddagger1, cddagger2)
  expr1 * (2 * pi * expr2) ^ (-1) * expr4 - A * B * expr1 * (2 * pi * expr2^3) ^ (-1/2) * expr3
}
