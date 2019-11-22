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

# Delong's test paired, used by roc.test.roc
delong.paired.test <- function(roc1, roc2) {
  n <- length(roc1$controls)
  m <- length(roc1$cases)

  VR <- delongPlacements(roc1)
  VS <- delongPlacements(roc2)

  SX <- matrix(NA, ncol=2, nrow=2)
  SX[1,1] <- sum((VR$X - VR$theta) * (VR$X - VR$theta))/(m-1)
  SX[1,2] <- sum((VR$X - VR$theta) * (VS$X - VS$theta))/(m-1)
  SX[2,1] <- sum((VS$X - VS$theta) * (VR$X - VR$theta))/(m-1)
  SX[2,2] <- sum((VS$X - VS$theta) * (VS$X - VS$theta))/(m-1)

  SY <- matrix(NA, ncol=2, nrow=2)
  SY[1,1] <- sum((VR$Y - VR$theta) * (VR$Y - VR$theta))/(n-1)
  SY[1,2] <- sum((VR$Y - VR$theta) * (VS$Y - VS$theta))/(n-1)
  SY[2,1] <- sum((VS$Y - VS$theta) * (VR$Y - VR$theta))/(n-1)
  SY[2,2] <- sum((VS$Y - VS$theta) * (VS$Y - VS$theta))/(n-1)

  S <- SX/m + SY/n
  L <- c(1,-1)
  sig <- sqrt(L%*%S%*%L)
  zscore <- (VR$theta-VS$theta)/sig[1]
  if (is.nan(zscore) && VR$theta == VR$theta && sig[1] == 0)
    zscore <- 0 # special case: no difference between theta's produces a NaN
  return(zscore)
}

# Delong's test unpaired, used by roc.test.roc
delong.unpaired.test <- function(roc1, roc2) {

  nR <- length(roc1$controls)
  mR <- length(roc1$cases)

  nS <- length(roc2$controls)
  mS <- length(roc2$cases)
  
  VR <- delongPlacements(roc1)
  VS <- delongPlacements(roc2)

  SRX <- sum((VR$X - VR$theta) * (VR$X - VR$theta))/(mR-1)
  SSX <- sum((VS$X - VS$theta) * (VS$X - VS$theta))/(mS-1)

  SRY <- sum((VR$Y - VR$theta) * (VR$Y - VR$theta))/(nR-1)
  SSY <- sum((VS$Y - VS$theta) * (VS$Y - VS$theta))/(nS-1)

  SR <- SRX/mR + SRY/nR
  SS <- SSX/mS + SSY/nS

  ntotR <- nR + mR
  ntotS <- nS + mS
  SSR <- sqrt((SR) + (SS))
  t <- (VR$theta - VS$theta) / SSR
  df <- ((SR) + (SS))^2 /
    (((SR)^2 / (ntotR-1)) + ((SS)^2 / (ntotS -1 )))

  return(c(t, df))
}

ci.auc.delong <- function(roc, conf.level) {
  YR <- roc$controls # = C2, n, YRj
  XR <- roc$cases # = C1, m, XRi

  n <- length(YR)
  m <- length(XR)
  
  # If controls or cases have a single observation, we would produce NaNs in SX and SY
  if (m <= 1 || n <= 1) {
  	return(rep(NA, 3))
  }

  V <- delongPlacements(roc)

  SX <- sum((V$X - V$theta) * (V$X - V$theta))/(m-1)
  SY <- sum((V$Y - V$theta) * (V$Y - V$theta))/(n-1)
  S <- SX/m + SY/n
  ci <- qnorm(c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2), mean = V$theta, sd = sqrt(S))
  # In some rare cases we have ci[3] > 1 or ci[1] < 0
  ci[ci > 1] <- 1
  ci[ci < 0] <- 0

  # According to Pepe (p. 107), we should probably be doing something like
  # log(roc$auc / (1 - roc$auc)) + pnorm( 1-conf.level/2) * (S / (roc$auc * (1 - roc$auc)))
  # log(roc$auc / (1 - roc$auc)) - pnorm( 1-conf.level/2) * (S / (roc$auc * (1 - roc$auc)))
  # for logit AUC, so that for AUC:
  # exp(log(roc$auc / (1 - roc$auc)) + pnorm( 1-conf.level/2) * (S / (roc$auc * (1 - roc$auc)))) * (1 - roc$auc)
  # exp(log(roc$auc / (1 - roc$auc)) - pnorm( 1-conf.level/2) * (S / (roc$auc * (1 - roc$auc)))) * (1 - roc$auc)
  # However the bounds are very very much smaller (about 10 times) than bootstrap, which seems unrealistic
  # Stay with normal conf interval for now.

  return(ci)
}

# Calls delongPlacementsCpp safely
# Ensures that the theta value calculated is correct
delongPlacements <- function(roc) {
	placements <- delongPlacementsCpp(roc)

	# Ensure theta equals auc
	auc <- roc$auc / ifelse(roc$percent, 100, 1)
	if (! isTRUE(all.equal(placements$theta, auc))) {
		sessionInfo <- sessionInfo()
		save(roc, placements, sessionInfo, file="pROC_bug.RData")
		stop(sprintf("pROC: error in calculating DeLong's theta: got %.20f instead of %.20f. Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", placements$theta, auc, utils:: packageDescription("pROC")$BugReports))
	}
	
	return(placements)
}