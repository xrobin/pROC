# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison.
# Copyright (C) 2010-2014 Xavier Robin, Alexandre Hainard, Natacha Turck,
# Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez,
# Markus Müller
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

# DeLong's test, paired or unpaired, used by roc.test.roc
delong.test <- function(calcs) {
  # Input calcs is a list returned by delong.paired.calculations()
  # or delong.unpaired.calculations().

  zscore <- with(calcs, d / sig)

  if (is.nan(zscore) && calcs$d == 0 && calcs$sig == 0) {
    zscore <- 0
  } # special case: no difference between theta's produces a NaN

  return(zscore)
}

# Runs the placements and main calculations for the unpaired DeLong's test.
# The curves are independent, so the variance of the difference is the sum
# of the variances of the two AUCs.
delong.unpaired.calculations <- function(roc1, roc2) {
  nR <- length(roc1$controls)
  mR <- length(roc1$cases)

  nS <- length(roc2$controls)
  mS <- length(roc2$cases)

  VR <- delongPlacements(roc1)
  VS <- delongPlacements(roc2)
  d <- VR$theta - VS$theta

  # A single case or control: the variance is undefined (as in ci_auc_delong)
  if (min(mR, nR, mS, nS) <= 1) {
    return(list("d" = d, "sig" = NA_real_))
  }

  SRX <- sum((VR$X - VR$theta) * (VR$X - VR$theta)) / (mR - 1)
  SSX <- sum((VS$X - VS$theta) * (VS$X - VS$theta)) / (mS - 1)

  SRY <- sum((VR$Y - VR$theta) * (VR$Y - VR$theta)) / (nR - 1)
  SSY <- sum((VS$Y - VS$theta) * (VS$Y - VS$theta)) / (nS - 1)

  SR <- SRX / mR + SRY / nR
  SS <- SSX / mS + SSY / nS

  return(list("d" = d, "sig" = sqrt(SR + SS)))
}

ci_auc_delong <- function(roc, conf.level) {
  YR <- roc$controls # = C2, n, YRj
  XR <- roc$cases # = C1, m, XRi

  n <- length(YR)
  m <- length(XR)

  # If controls or cases have a single observation, we would produce NaNs in SX and SY
  if (m <= 1 || n <= 1) {
    return(rep(NA, 3))
  }

  V <- delongPlacements(roc)

  SX <- sum((V$X - V$theta) * (V$X - V$theta)) / (m - 1)
  SY <- sum((V$Y - V$theta) * (V$Y - V$theta)) / (n - 1)
  S <- SX / m + SY / n
  # S == 0 (AUC of 0 or 1, or all values tied): no interval; the logit interval would be NaN.
  # conf.level == 0: no interval either, and plogis(qlogis(theta)) can be 1 ulp off theta
  if (S == 0 || conf.level == 0) {
    return(rep(V$theta, 3))
  }
  # Normal interval on the logit scale, which stays within (0, 1) (Pepe 2003, p. 107)
  half <- qnorm(1 - (1 - conf.level) / 2) * sqrt(S) / (V$theta * (1 - V$theta))
  ci <- c(plogis(qlogis(V$theta) - half), V$theta, plogis(qlogis(V$theta) + half))

  return(ci)
}

# Runs the placements and main calculations for the paired DeLong's test
#  so that they can be easily used by both the test and CI functions.
delong.paired.calculations <- function(roc1, roc2) {
  calcs <- delong.cov(list(roc1, roc2))
  L <- c(1, -1)
  # A single case or control: cov() and thus sig are NA (as in ci_auc_delong)
  list("d" = sum(L * calcs$theta), "sig" = sqrt(drop(L %*% calcs$S %*% L)))
}

# AUCs (theta) of k paired ROC curves and their covariance matrix S,
# DeLong et al. (1988) p. 840: S = S10 / m + S01 / n
delong.cov <- function(rocs) {
  V <- lapply(rocs, delongPlacements)
  X <- do.call(cbind, lapply(V, `[[`, "X")) # V10, m x k
  Y <- do.call(cbind, lapply(V, `[[`, "Y")) # V01, n x k
  list(
    theta = vapply(V, `[[`, numeric(1), "theta"),
    S = stats::cov(X) / nrow(X) + stats::cov(Y) / nrow(Y)
  )
}

# Validates a list of curves for the k-curve DeLong methods (roc.test.list,
# cov.list). Returns them as roc objects carrying their AUC.
delong_list_rocs <- function(rocs) {
  if (length(rocs) < 2) {
    stop("At least 2 ROC curves are required.")
  }
  rocs <- lapply(rocs, function(roc) {
    if (methods::is(roc, "auc")) {
      # take the roc but keep the auc specifications
      auc <- roc
      roc <- attr(auc, "roc")
      attr(auc, "roc") <- NULL
      roc$auc <- auc
    }
    if (methods::is(roc, "smooth.roc")) {
      stop("DeLong's method is not supported for smoothed ROCs.")
    }
    if (!methods::is(roc, "roc")) {
      stop("All the elements of the list must be 'roc' or 'auc' objects.")
    }
    if (is.null(roc$auc)) {
      roc$auc <- auc(roc)
    }
    if (has.partial.auc(roc)) {
      stop("DeLong's method is not supported for partial AUC.")
    }
    roc
  })
  for (roc in rocs[-1]) {
    if (!are.paired(rocs[[1]], roc)) {
      stop("DeLong's method for a list of ROC curves requires paired ROC curves.")
    }
    if (!identical(unname(rocs[[1]]$response), unname(roc$response))) {
      stop("The ROC curves have different missing values. Build them on the observations complete in all the predictors.")
    }
    if (rocs[[1]]$direction != roc$direction) {
      warning("DeLong's method should not be applied to ROC curves with a different direction.")
    }
  }
  rocs
}

# Calls delongPlacementsCpp safely
# Ensures that the theta value calculated is correct
delongPlacements <- function(roc) {
  if (roc_utils_is_ordered_roc(roc)) {
    roc$cases <- as.numeric(roc$cases)
    roc$controls <- as.numeric(roc$controls)
  }
  placements <- delongPlacementsCpp(roc)

  # Ensure theta equals auc
  auc <- roc$auc / ifelse(roc$percent, 100, 1)
  if (!isTRUE(all.equal(placements$theta, auc))) {
    roc_utils_bug(sprintf("error in calculating DeLong's theta: got %.20f instead of %.20f", placements$theta, auc), roc = roc, placements = placements)
  }

  return(placements)
}
