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

power.roc.test <- function(...)
  UseMethod("power.roc.test")

power.roc.test.roc <- function(roc1, roc2, sig.level = 0.05, power = NULL, kappa = NULL,
							   alternative = c("two.sided", "one.sided"), 
							   reuse.auc=TRUE, method=c("delong", "bootstrap", "obuchowski"),
							   ...) {
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
  roc1 <- roc.utils.unpercent(roc1)

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
      roc2 <- roc.utils.unpercent(roc2)

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
    
      # check that the same region was requested in auc. Otherwise, issue a warning
      if (!identical(attributes(roc1$auc)[names(attributes(roc1$auc))!="roc"], attributes(roc2$auc)[names(attributes(roc2$auc))!="roc"]))
        warning("Different AUC specifications in the ROC curves. Enforcing the inconsistency, but unexpected results may be produced.")
 
      ncontrols <- length(roc1$controls)
      ncases <- length(roc1$cases)
      if (is.null(kappa)) {
        kappa <- ncontrols / ncases
      }

      # Power test
      if (is.null(power)) {
        if (is.null(sig.level))
          stop("'sig.level' or 'power' must be provided.")
        zalpha <- qnorm(1 - sig.level)
        zbeta <- zbeta.obuchowski(roc1, roc2, zalpha, method=method, ...)
        power <- pnorm(zbeta)
      }
      # sig.level
      else if (is.null(sig.level)) {
        zbeta <- qnorm(power)
        zalpha <- zalpha.obuchowski(roc1, roc2, zbeta, method=method, ...)
        sig.level <- 1 - pnorm(zalpha)
      }
      # Sample size
      else {
        zalpha <- qnorm(1 - sig.level)
        zbeta <- qnorm(power)
        ncases <- ncases.obuchowski(roc1, roc2, zalpha, zbeta, method=method, ...)
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
  	ncontrols <- length(roc1$controls)
  	ncases <- length(roc1$cases)
    if (! is.null(sig.level) && ! is.null(power)) {
      if (is.null(kappa)) {
        kappa <- ncontrols / ncases
      }
      ncontrols <- ncases <- NULL
    }
    auc <- auc(roc1)
    # TODO: implement this with var() and cov() for the given ROC curve
    return(power.roc.test.numeric(ncontrols = ncontrols, ncases = ncases, auc = auc, sig.level = sig.level, power = power, alternative = alternative, kappa = kappa, ...))
  }
}

power.roc.test.numeric <- function(auc = NULL, ncontrols = NULL, ncases = NULL, sig.level = 0.05, power = NULL,  kappa = 1, alternative = c("two.sided", "one.sided"), ...) {
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
  
  # Complete ncontrols and ncases with kappa
  if (is.null(ncontrols) && ! is.null(ncases) && !is.null(kappa))
    ncontrols <- kappa * ncases
  else if (is.null(ncases) && ! is.null(ncontrols) && !is.null(kappa))
    ncases <- ncontrols / kappa
  
  alternative <- match.arg(alternative)
  if (alternative == "two.sided" && !is.null(sig.level)) {
    sig.level <- sig.level / 2
  }

  # determine AUC
  if (is.null(auc)) {
    if (is.null(ncontrols) || is.null(ncases))
      stop("'ncontrols' and 'ncases' (or one of these with 'kappa') or 'auc' must be provided.")
    else if (is.null(power))
      stop("'power' or 'auc' must be provided.")
    else if (is.null(sig.level))
      stop("'sig.level' or 'auc' must be provided.")
    kappa <- ncontrols / ncases
    zalpha <- qnorm(1 - sig.level)
    zbeta <- qnorm(power)

    tryCatch(
             root <- uniroot(power.roc.test.optimize.auc.function, interval=c(0.5, 1-1e-16), ncontrols=ncontrols, ncases=ncases, zalpha=zalpha, zbeta=zbeta),
             error=function(e) {stop(sprintf("AUC could not be solved:\n%s", e))}
             )
    auc <- root$root
  }

  # Determine number of patients (sample size)
  else if (is.null(ncases) && is.null(ncontrols)) {
    if (is.null(power))
      stop("'ncontrols' and 'ncases' (or one of these with 'kappa') or 'auc' must be provided.")
    else if (is.null(kappa))
      stop("'kappa' must be provided.")
    else if (is.null(auc))
      stop("'auc' or 'ncases' and 'ncontrols' must be provided.")
    else if (is.null(sig.level))
      stop("'sig.level' or 'ncases' and 'ncontrols' must be provided.")

    theta <- as.numeric(auc)
    Vtheta <- var.theta.obuchowski(theta, kappa)
    
    ncases <- solve.nd(zalpha = qnorm(1 - sig.level),
    				   zbeta = qnorm(power),
    				   v0 = 0.0792 * (1 + 1 / kappa),
    				   va = Vtheta,
    				   delta = theta - 0.5)
    ncontrols <- kappa * ncases
  }
  
  # Determine power
  else if (is.null(power)) { 
    if (is.null(ncontrols) || is.null(ncases))
      stop("'ncontrols' and 'ncases' (or one of these with 'kappa') or 'auc' must be provided.")
    else if (is.null(auc))
      stop("'auc' or 'power' must be provided.")
    else if (is.null(sig.level))
      stop("'sig.level' or 'power' must be provided.")
  	
    kappa <- ncontrols / ncases
    theta <- as.numeric(auc)
    Vtheta <- var.theta.obuchowski(theta, kappa)
    
    zbeta <- solve.zbeta(nd = ncases,
    					 zalpha = qnorm(1 - sig.level),
    					 v0 = 0.0792 * (1 + 1 / kappa),
    					 va = Vtheta,
    					 delta = theta - 0.5)
    power <- pnorm(zbeta)
  }

  # Determine sig.level
  else  if (is.null(sig.level)) { 
    if (is.null(ncontrols) || is.null(ncases))
      stop("'ncontrols' and 'ncases' (or one of these with 'kappa') or 'auc' must be provided.")
    else if (is.null(auc))
      stop("'auc' or 'sig.level' must be provided.")
    else if (is.null(power))
      stop("'power' or 'sig.level' must be provided.")
  	
    kappa <- ncontrols / ncases
    theta <- as.numeric(auc)
    Vtheta <- var.theta.obuchowski(theta, kappa)
    
    zalpha <- solve.zalpha(nd = ncases,
    					  zbeta = qnorm(power),
    					  v0 = 0.0792 * (1 + 1 / kappa),
    					  va = Vtheta,
    					  delta = theta - 0.5)
    sig.level <- 1 - pnorm(zalpha)
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

power.roc.test.list <- function(parslist, ncontrols = NULL, ncases = NULL, sig.level = 0.05, power = NULL,  kappa = 1, alternative = c("two.sided", "one.sided"), ...) {
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

  
  # Complete ncontrols and ncases with kappa
  if (is.null(ncontrols) && ! is.null(ncases) && !is.null(kappa))
    ncontrols <- kappa * ncases
  else if (is.null(ncases) && ! is.null(ncontrols) && !is.null(kappa))
    ncases <- ncontrols / kappa

  # Warn if anything is passed with ...
  if (length(list(...)) > 0) {
    warning(paste("The following arguments were ignored:", paste(names(list(...)), collapse=", ")))
  }
  
  alternative <- match.arg(alternative)
  if (alternative == "two.sided" && !is.null(sig.level)) {
    sig.level <- sig.level / 2
  }

  # Check required elements of parslist
  required <- c("A1", "B1", "A2", "B2", "rn", "ra", "delta")
  if (any(! required %in% names(parslist))) {
    stop(paste("Missing parameter(s):", paste(required[! required %in% names(parslist) ], collapse=", ")))
  }

  # Determine number of patients (sample size)
  if (is.null(ncases) && is.null(ncontrols)) {
    if (is.null(power))
      stop("'ncontrols' and 'ncases' (or one of these with 'kappa') or 'auc' must be provided.")
    else if (is.null(kappa))
      stop("'kappa' must be provided.")
    else if (is.null(sig.level))
      stop("'sig.level' or 'ncases' and 'ncontrols' must be provided.")

    zalpha <- qnorm(1 - sig.level)
    zbeta <- qnorm(power)
    ncases <- ncases.obuchowski.params(parslist, zalpha, zbeta, kappa)
    ncontrols <- kappa * ncases
  }
  
  # Determine power
  else if (is.null(power)) {
    if (is.null(ncontrols) || is.null(ncases))
      stop("'ncontrols' and 'ncases' (or one of these with 'kappa') or 'auc' must be provided.")
    else if (is.null(sig.level))
      stop("'sig.level' or 'power' must be provided.")
    kappa <- ncontrols / ncases

    zalpha <- qnorm(1 - sig.level)
    zbeta <- zbeta.obuchowski.params(parslist, zalpha, ncases, kappa)
    power <- pnorm(zbeta)
  }

  # Determine sig.level
  else  if (is.null(sig.level)) { 
    if (is.null(ncontrols) || is.null(ncases))
      stop("'ncontrols' and 'ncases' (or one of these with 'kappa') or 'auc' must be provided.")
    else if (is.null(power))
      stop("'power' or 'sig.level' must be provided.")
    kappa <- ncontrols / ncases

    zbeta <- qnorm(power)
    zalpha <- zalpha.obuchowski.params(parslist, zbeta, ncases, kappa)
    sig.level <- 1 - pnorm(zalpha)
  }
  else {
    stop("One of 'power', 'sig.level', 'auc', or both 'ncases' and 'ncontrols' must be NULL.")
  }
  # Restore sig.level if two.sided
  if (alternative == "two.sided") {
    sig.level <- sig.level * 2
  }
  return(structure(list(ncases=ncases, ncontrols=ncontrols, sig.level=sig.level, power=power, method="Two ROC curves power calculation"), class="power.htest"))

}


#### HIDDEN FUNCTIONS ####

# A function to 'optimize' auc
power.roc.test.optimize.auc.function <- function(x, ncontrols, ncases, zalpha, zbeta) {
  kappa <- ncontrols / ncases
  Vtheta <- var.theta.obuchowski(x, kappa)
  (zalpha * sqrt(0.0792 * (1 + 1/kappa)) + zbeta * sqrt(Vtheta))^2 / (x - 0.5)^2 - ncases
}

# Compute variance of a delta from a 'covvar' list (see 'covvar' below)
var.delta.covvar <- function(covvar) {
  covvar$var1 + covvar$var2 - 2 * covvar$cov12
}

# Compute variance of a delta from a 'covvar' list (see 'covvar' below)
# under the null hypothesis
# roc1 taken as reference.
var0.delta.covvar <- function(covvar) {
  2 * covvar$var1 - 2 * covvar$cov12
}

# Compute the number of cases with Obuchowski formula and var(... method=method)
ncases.obuchowski <- function(roc1, roc2, zalpha, zbeta, method, ...) {
  delta <- roc1$auc - roc2$auc
  covvar <- covvar(roc1, roc2, method, ...)
  v0 <- var0.delta.covvar(covvar)
  va <- var.delta.covvar(covvar)
  nd <- solve.nd(zalpha = zalpha,
  			   zbeta = zbeta,
  			   v0 = v0, va = va,
  			   delta = delta)
  return(nd)
}

# Compute the number of cases with Obuchowski formula from params
ncases.obuchowski.params <- function(parslist, zalpha, zbeta, kappa) {
  covvar <- list(
                 var1 = var.params.obuchowski(parslist$A1, parslist$B1, kappa, parslist$FPR11, parslist$FPR12),
                 var2 = var.params.obuchowski(parslist$A2, parslist$B2, kappa, parslist$FPR21, parslist$FPR22),
                 cov12 = cov.params.obuchowski(parslist$A1, parslist$B1, parslist$A2, parslist$B2, parslist$rn, parslist$ra, kappa, parslist$FPR11, parslist$FPR12, parslist$FPR21, parslist$FPR22)
                 )
  v0 <- var0.delta.covvar(covvar)
  va <- var.delta.covvar(covvar)
  nd <- solve.nd(zalpha = zalpha,
  			   zbeta = zbeta,
  			   v0 = v0, va = va,
  			   delta = parslist$delta)
  return(nd)
}

# Compute the z alpha with Obuchowski formula and var(... method=method)
zalpha.obuchowski <- function(roc1, roc2, zbeta, method, ...) {
  delta <- roc1$auc - roc2$auc
  ncases <- length(roc1$cases)
  covvar <- covvar(roc1, roc2, method, ...)
  v0 <- var0.delta.covvar(covvar)
  va <- var.delta.covvar(covvar)
  zalpha <- solve.zalpha(nd=ncases,
  					   zbeta = zbeta,
  					   v0 = v0, va = va,
  					   delta = delta)
  return(zalpha)
}

# Compute the z alpha with Obuchowski formula from params
zalpha.obuchowski.params <- function(parslist, zbeta, ncases, kappa) {
  covvar <- list(
                 var1 = var.params.obuchowski(parslist$A1, parslist$B1, kappa, parslist$FPR11, parslist$FPR12),
                 var2 = var.params.obuchowski(parslist$A2, parslist$B2, kappa, parslist$FPR21, parslist$FPR22),
                 cov12 = cov.params.obuchowski(parslist$A1, parslist$B1, parslist$A2, parslist$B2, parslist$rn, parslist$ra, kappa, parslist$FPR11, parslist$FPR12, parslist$FPR21, parslist$FPR22)
                 )
  v0 <- var0.delta.covvar(covvar)
  va <- var.delta.covvar(covvar)
  zalpha <- solve.zalpha(nd=ncases,
  					   zbeta = zbeta,
  					   v0 = v0, va = va,
  					   delta = parslist$delta)
  return(zalpha)
}

# Compute the z beta with Obuchowski formula and var(... method=method)
zbeta.obuchowski <- function(roc1, roc2, zalpha, method, ...) {
  delta <- roc1$auc - roc2$auc
  ncases <- length(roc1$cases)
  covvar <- covvar(roc1, roc2, method, ...)
  v0 <- var0.delta.covvar(covvar)
  va <- var.delta.covvar(covvar)
  zbeta <- solve.zbeta(nd=ncases,
  					 zalpha = zalpha,
  					 v0 = v0, va = va,
  					 delta = delta)
  return(zbeta)
}

# Compute the z beta with Obuchowski formula from params
zbeta.obuchowski.params <- function(parslist, zalpha, ncases, kappa) {
	covvar <- list(
		var1 = var.params.obuchowski(parslist$A1, parslist$B1, kappa, parslist$FPR11, parslist$FPR12),
		var2 = var.params.obuchowski(parslist$A2, parslist$B2, kappa, parslist$FPR21, parslist$FPR22),
		cov12 = cov.params.obuchowski(parslist$A1, parslist$B1, parslist$A2, parslist$B2, parslist$rn, parslist$ra, kappa, parslist$FPR11, parslist$FPR12, parslist$FPR21, parslist$FPR22)
	)
	v0 <- var0.delta.covvar(covvar)
	va <- var.delta.covvar(covvar)
	a <- va
	zbeta <- solve.zbeta(nd=ncases,
						 zalpha = zalpha,
						 v0 = v0, va = va,
						 delta = parslist$delta)
	return(zbeta)
}

solve.zbeta <- function(nd, zalpha, v0, va, delta) {
	# Solve for z_\beta in Obuchowski formula:
	# See formula 2 in Obuchowsk & McClish 1997  (2 ROC curves)
	# or formula 2 in Obuchowski et al 2004 (1 ROC curve)
	# The formula is of the form:
	# nd = (z_alpha * sqrt(v0) - z_beta * sqrt(va)) / delta ^ 2
	# Re-organized:
	# z_beta = (sqrt(nd * delta ^ 2) - z_alpha * sqrt(v0)) / sqrt(va)
	# @param nd: number of diseased patients (or abornmal, N_A in Obuchowsk & McClish 1997)
	# @param zalpha: upper \alpha (sig.level) percentile of the standard normal distribution
	# @param v0 the null variance associated with z_alpha
	# @param va: the alternative variance associated with z_beta
	# @param delta: the difference in AUC
	return((sqrt(nd * delta ^ 2) - zalpha * sqrt(v0)) / sqrt(va))
}

solve.nd <- function(zalpha, zbeta, v0, va, delta) {
	# Solve for number of diseased (abnormal) patients in Obuchowski formula:
	# See formula 2 in Obuchowsk & McClish 1997  (2 ROC curves)
	# or formula 2 in Obuchowski et al 2004 (1 ROC curve)
	# nd = (z_alpha * sqrt(v0) - z_beta * sqrt(va)) / delta ^ 2
	# @param zalpha: upper \alpha (sig.level) percentile of the standard normal distribution
	# @param zbeta: upper \beta (power) percentile of the standard normal distribution
	# @param v0 the null variance associated with z_alpha
	# @param va: the alternative variance associated with z_beta
	# @param delta: the difference in AUC
	return((zalpha * sqrt(v0) + zbeta * sqrt(va)) ^ 2 / delta ^ 2)
}

solve.zalpha <- function(nd, zbeta, v0, va, delta) {
	# Solve for z_\alpha in Obuchowski formula:
	# See formula 2 in Obuchowsk & McClish 1997  (2 ROC curves)
	# or formula 2 in Obuchowski et al 2004 (1 ROC curve)
	# The formula is of the form:
	# nd = (z_alpha * sqrt(v0) - z_beta * sqrt(va)) / delta ^ 2
	# Re-organized:
	# z_alpha = (sqrt(nd * delta ^ 2) - z_beta * sqrt(va)) / sqrt(v0)
	# @param nd: number of diseased patients (or abornmal, N_A in Obuchowsk & McClish 1997)
	# @param zbeta: upper \beta (power) percentile of the standard normal distribution
	# @param v0 the null variance associated with z_alpha
	# @param va: the alternative variance associated with z_beta
	# @param delta: the difference in AUC
	return((sqrt(nd * delta ^ 2) - zbeta * sqrt(va)) / sqrt(v0))
}

# Compute var and cov of two ROC curves by bootstrap in a single bootstrap run
covvar <- function(roc1, roc2, method, ...) {
  cov12 <- cov(roc1, roc2, boot.return=TRUE, method=method, ...)
  if (!is.null(attr(cov12, "resampled.values"))) {
    var1 <- var(attr(cov12, "resampled.values")[,1])
    var2 <- var(attr(cov12, "resampled.values")[,2])
    attr(cov12, "resampled.values") <- NULL
  }
  else {
    var1 <- var(roc1, method=method, ...)
    var2 <- var(roc2, method=method, ...)
  }
  ncases <- length(roc1$cases)
  return(list(var1 = var1 * ncases, var2 = var2 * ncases, cov12 = cov12 * ncases))
}
