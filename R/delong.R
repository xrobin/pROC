# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2013 Xavier Robin, Alexandre Hainard, Natacha Turck,
# Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez,
# Markus Müller and Kazuki Yoshida
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

  delong.warn(roc1, roc2)
  VR <- delong.placements(roc1)
  VS <- delong.placements(roc2)

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
  
  delong.warn(roc1, roc2)
  VR <- delong.placements(roc1)
  VS <- delong.placements(roc2)

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
  mn <- m*n
  
  delong.warn(roc)
  V <- delong.placements(roc)

  SX <- sum((V$X - V$theta) * (V$X - V$theta))/(m-1)
  SY <- sum((V$Y - V$theta) * (V$Y - V$theta))/(n-1)
  S <- SX/m + SY/n
  ci <- qnorm(c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2), mean = V$theta, sd = sqrt(S))
  if (roc$direction == ">") {
    ci <- rev(1 - ci)
  }
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

delong.placements <- function(roc) {
  # returns a list V containing:
  # - theta: the AUC
  # - X: the 10 component
  # - Y: the 01 component
  V <- list()
  Y <- roc$controls
  X <- roc$cases
  n <- length(Y)
  m <- length(X)
  
  # Original computation of MW matrix as given in DeLong et al paper
  # MW <- sapply(1:n, function(j) sapply(1:m, function(i, j) MW.kernel(X[i], Y[j]), j=j))
  # Alternative version by Kazuki Yoshida
  #equal   <- outer(X, Y, "==") * 0.5
  #greater <- outer(X, Y, ">") * 1.0
  MW <- outer(X, Y, "==") * 0.5 + outer(X, Y, ">") * 1.0
  
  V$theta <- sum(MW)/(m*n)
  # Delong-specific computations
  V$X <- sapply(1:m, function(i) {sum(MW[i,])})/n
  V$Y <- sapply(1:n, function(j) {sum(MW[,j])})/m
  return(V)
}

delong.warn <- function(roc, roc2) {
  # Safety function to be called before roc.placements.
  # Will warn and ask for user input. If user supplies any string that starts with N, it will abort the processing
  # Based on the assumption that 3 X*Y matrices will be created, each matrix being X*Y*8 bits at least
  XY <- length(roc$controls) * length(roc$cases)
  if (!missing(roc2)) {
    XY2 <- length(roc2$controls) * length(roc2$cases)
    XY <- max(XY, XY2)
  }
  bytes <- XY * 3 * 8
  if (bytes > 1E9) {
    ans <- readline(sprintf("This could allocate more than %.1f GB of memory. Are you sure you want to continue?\nType N to abort, anything else to continue... ", bytes / 1E9))
    if (grepl("^\\s?n", ans)) {
      stop("Command aborted on user request. Consider using method=\"bootstrap\" instead.", call.=FALSE)
    }
  }
}
