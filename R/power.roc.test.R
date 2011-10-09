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
    if (is.numeric(attr(roc1$auc, "partial.auc"))) {
      attr(roc1$auc, "partial.auc") <- attr(roc1$auc, "partial.auc") / 100
    }
    attr(roc1$auc, "percent") <- FALSE
    roc1$percent <- FALSE
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
        if (is.numeric(attr(roc2$auc, "partial.auc"))) {
          attr(roc2$auc, "partial.auc") <- attr(roc2$auc, "partial.auc") / 100
        }
        attr(roc2$auc, "percent") <- FALSE
        roc2$percent <- FALSE
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
        zbeta <- zbeta.bootdelong(roc1, roc2, zalpha)
        power <- 1 - pnorm(zbeta)
      }
      # sig.level
      else if (is.null(sig.level)) {
        zbeta <- qnorm(1 - power)
        zalpha <- zalpha.bootdelong(roc1, roc2, zbeta)
        sig.level <- pnorm(zalpha)
      }
      # Sample size
      else {
        zalpha <- qnorm(sig.level)
        zbeta <- qnorm(1 - power)

        ncases <- ncases.bootdelong(roc1, roc2, zalpha, zbeta)
        ncontrols <- kappa * ncases
      }

      # Restore sig.level if two.sided
      if (alternative == "two.sided") {
        sig.level <- sig.level * 2
      }
      return(structure(list(ncases=ceiling(ncases), ncontrols=ceiling(ncontrols), auc1=roc1$auc, auc2=roc2$auc, sig.level=sig.level, power=power, alternative=alternative, method="Two ROC curves power calculation"), class="power.htest"))
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
    else {
      ncontrols <- ncases <- NULL
    }
    auc <- roc1$auc
    return(power.roc.test.numeric(ncontrols = ncontrols, ncases = ncases, auc = auc, sig.level = sig.level, power = power, alternative = alternative, ...))
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
  return(structure(list(ncases=ceiling(ncases), ncontrols=ceiling(ncontrols), auc=auc, sig.level=sig.level, power=power, method="One ROC curve power calculation"), class="power.htest"))
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

var.delta.bootdelong <- function(covvar) {
  covvar$var1 + covvar$var2 - 2 * covvar$cov12
}

var0.delta.bootdelong <- function(covvar) {
  if (covvar$var1 < covvar$var2) {
    varroc <- covvar$var2
  }
  else {
    varroc <- covvar$var1
  }
  2 * varroc - 2 * covvar$cov12
}

ncases.bootdelong <- function(roc1, roc2, zalpha, zbeta) {
  delta <- roc1$auc - roc2$auc
  covvar <- covvar(roc1, roc2)
  na <- (zalpha * sqrt(var0.delta.bootdelong(covvar)) +
       zbeta * sqrt(var.delta.bootdelong(covvar))) ^2 /
       delta^2
  return(as.vector(na))
}

zalpha.bootdelong <- function(roc1, roc2, zbeta) {
  delta <- roc1$auc - roc2$auc
  ncases <- length(roc1$cases)
  covvar <- covvar(roc1, roc2)
  v0 <- var0.delta.bootdelong(covvar)
  va <- var.delta.bootdelong(covvar)
  a <- v0
  b <- 2 * zbeta * sqrt(v0) * sqrt(va)
  c <- zbeta^2 * va - ncases * delta ^ 2
  return(as.vector(solve.2deg.eqn(a, b, c)))
}

zbeta.bootdelong <- function(roc1, roc2, zalpha) {
  delta <- roc1$auc - roc2$auc
  ncases <- length(roc1$cases)
  covvar <- covvar(roc1, roc2)
  v0 <- var0.delta.bootdelong(covvar)
  va <- var.delta.bootdelong(covvar)
  a <- va
  b <- 2 * zalpha * sqrt(va) * sqrt(v0)
  c <- zalpha^2 * v0 - ncases * delta ^ 2
  return(as.vector(solve.2deg.eqn(a, b, c)))
}

solve.2deg.eqn <- function(a, b, c) {
  return((- b - sqrt(b^2 - 4*a*c)) / (2*a))
}

covvar <- function(roc1, roc2) {
  cov12 <- cov(roc1, roc2, boot.return=TRUE)
  if (!is.null(attr(cov12, "resampled.values"))) {
    var1 <- var(attr(cov12, "resampled.values")[,1])
    var2 <- var(attr(cov12, "resampled.values")[,2])
    attr(cov12, "resampled.values") <- NULL
  }
  else {
    var1 <- var(roc1)
    var2 <- var(roc2)
  }
  return(list(var1 = var1, var2 = var2, cov12 = cov12))
}
