# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2010-2014 Xavier Robin, Alexandre Hainard, Natacha Turck,
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

venkatraman.paired.test <- function(roc1, roc2, boot.n, ties.method="first", progress, parallel) {
  X <- roc1$predictor
  Y <- roc2$predictor
  R <- rank(X, ties.method = ties.method)
  S <- rank(Y, ties.method = ties.method)
  D <- roc1$response # because roc1&roc2 are paired

  E <- venkatraman.paired.stat(R, S, D, roc1$levels)
  EP <- laply(seq_len(boot.n), venkatraman.paired.permutation, R=R, S=S, D=D, levels=roc1$levels, ties.method=ties.method, .progress=progress, .parallel=parallel)
  return(list(E, EP))
}

venkatraman.unpaired.test <- function(roc1, roc2, boot.n, ties.method="first", progress, parallel) {
  X <- roc1$predictor
  Y <- roc2$predictor
  R <- rank(X, ties.method = ties.method)
  S <- rank(Y, ties.method = ties.method)
  D1<- roc1$response
  D2 <- roc2$response
  mp <- (sum(D1 == roc1$levels[2]) + sum(D2 == roc2$levels[2]))/(length(D1) + length(D1)) # mixing proportion, kappa

  E <- venkatraman.unpaired.stat(R, S, D1, D2, roc1$levels, roc2$levels, mp)
  EP <- laply(seq_len(boot.n), venkatraman.unpaired.permutation, R=R, S=S, D1=D1, D2=D2, levels1=roc1$levels, levels2=roc2$levels, mp=mp, ties.method=ties.method, .progress=progress, .parallel=parallel)
  return(list(E, EP))
}

venkatraman.paired.permutation <- function(n, R, S, D, levels, ties.method) {
  # Break ties
  R2 <- R + runif(length(D)) - 0.5 # Add small amount of random but keep same mean
  S2 <- S + runif(length(D)) - 0.5

  # Permutation
  q <- 1 - round(runif(length(D)))
  R3 <- R2 * q + (1 - q) * S
  S3 <- S2 * q + (1 - q) * R

  return(venkatraman.paired.stat(rank(R3, ties.method=ties.method), rank(S3, ties.method=ties.method), D, levels))
}

venkatraman.unpaired.permutation <- function(n, R, S, D1, D2, levels1, levels2, mp, ties.method) {
  # Break ties
  R <- R + runif(length(D1)) - 0.5 # Add small amount of random but keep same mean
  S <- S + runif(length(D2)) - 0.5

  R.controls <- R[D1==levels1[1]]
  R.cases <- R[D1==levels1[2]]
  S.controls <- S[D2==levels2[1]]
  S.cases <- S[D2==levels2[2]]

  # Permutation
  controls <- sample(c(R.controls, S.controls))
  cases <- sample(c(R.cases, S.cases))
  R[D1==levels1[1]] <- controls[1:length(R.controls)]
  S[D2==levels2[1]] <- controls[(length(R.controls)+1):length(controls)]
  R[D1==levels1[2]] <- cases[1:length(R.cases)]
  S[D2==levels2[2]] <- cases[(length(R.cases)+1):length(cases)]

  return(venkatraman.unpaired.stat(rank(R, ties.method=ties.method), rank(S, ties.method=ties.method), D1, D2, levels1, levels2, mp))
}

venkatraman.paired.stat <- function(R, S, D, levels) {
  R.controls <- R[D==levels[1]]
  R.cases <- R[D==levels[2]]
  S.controls <- S[D==levels[1]]
  S.cases <- S[D==levels[2]]
  n <- length(D)

  R.fn <- sapply(1:n, function(x) sum(R.cases <= x))
  R.fp <- sapply(1:n, function(x) sum(R.controls > x))
  S.fn <- sapply(1:n, function(x) sum(S.cases <= x))
  S.fp <- sapply(1:n, function(x) sum(S.controls > x))

  return(sum(abs((S.fn + S.fp) - (R.fn + R.fp))))
}

venkatraman.unpaired.stat <- function(R, S, D1, D2, levels1, levels2, mp) {
  R.controls <- R[D1==levels1[1]]
  R.cases <- R[D1==levels1[2]]
  S.controls <- S[D2==levels2[1]]
  S.cases <- S[D2==levels2[2]]
  n <- length(D1)
  m <- length(D2)

  R.fx <- sapply(1:n, function(x) sum(R.cases <= x)) / length(R.cases)
  R.gx <- sapply(1:n, function(x) sum(R.controls <= x)) / length(R.controls)
  S.fx <- sapply(1:m, function(x) sum(S.cases <= x)) / length(S.cases)
  S.gx <- sapply(1:m, function(x) sum(S.controls <= x)) / length(S.controls)
  R.p <- mp*R.fx + (1 - mp)*R.gx
  S.p <- mp*S.fx + (1 - mp)*S.gx
  R.exp <- mp*R.fx + (1 - mp)*(1-R.gx)
  S.exp <- mp*S.fx + (1 - mp)*(1-S.gx)

  # Do the integration
  x <- sort(c(R.p, S.p))
  R.f <- approxfun(R.p, R.exp)
  S.f <- approxfun(S.p, S.exp)
  f  <- function(x) abs(R.f(x)-S.f(x))
  y <- f(x)
  #trapezoid integration:
  idx <- 2:length(x)
  integral <- sum(((y[idx] + y[idx-1]) * (x[idx] - x[idx-1])) / 2, na.rm=TRUE) # remove NA that can appear in the borders
  return(integral)
}
