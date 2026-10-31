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

# Venkatraman and Begg (1996) define the paired test for continuous data,
# through the ranks of both predictors. The only tie-breaking they use is
# random: the ties created by the permutation are "re-ranked by
# randomisation" (J(.), p. 838), which venkatraman.paired.permutation() does
# with its jitter. Ties in the data themselves are outside the method (their
# section 5); with ties.method = "random" they are broken in one random order
# of the observations, shared by both predictors, so that an observation tied
# in both is ranked the same way in both (identical predictors give E = 0).
venkatraman.paired.test <- function(roc1, roc2, boot.n, ties.method = "random") {
  X <- roc1$predictor
  Y <- roc2$predictor
  if (identical(ties.method, "random")) {
    # No random numbers without ties: untied data give the ranks (and use the
    # random stream) of ties.method = "first"
    o <- if (anyDuplicated(X) || anyDuplicated(Y)) sample.int(length(X)) else seq_along(X)
    R <- S <- numeric(length(X))
    R[o] <- rank(X[o], ties.method = "first")
    S[o] <- rank(Y[o], ties.method = "first")
  } else {
    R <- rank(X, ties.method = ties.method)
    S <- rank(Y, ties.method = ties.method)
  }
  D <- roc1$response # because roc1&roc2 are paired

  E <- venkatraman.paired.stat(R, S, D, roc1$levels)
  EP <- vapply(seq_len(boot.n), venkatraman.paired.permutation, FUN.VALUE = double(1L), R = R, S = S, D = D, levels = roc1$levels, ties.method = ties.method)
  return(list(E, EP))
}

venkatraman.unpaired.test <- function(roc1, roc2, boot.n, ties.method = "random") {
  X <- roc1$predictor
  Y <- roc2$predictor
  R <- venkatraman.rank(X, ties.method)
  S <- venkatraman.rank(Y, ties.method)
  D1 <- roc1$response
  D2 <- roc2$response
  mp <- (sum(D1 %in% roc1$levels[2]) + sum(D2 %in% roc2$levels[2])) / (length(D1) + length(D2)) # mixing proportion, kappa

  E <- venkatraman.unpaired.stat(R, S, D1, D2, roc1$levels, roc2$levels, mp)
  EP <- vapply(seq_len(boot.n), venkatraman.unpaired.permutation, FUN.VALUE = double(1L), R = R, S = S, D1 = D1, D2 = D2, levels1 = roc1$levels, levels2 = roc2$levels, mp = mp, ties.method = ties.method)
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

  return(venkatraman.paired.stat(venkatraman.rank(R3, ties.method), venkatraman.rank(S3, ties.method), D, levels))
}

venkatraman.unpaired.permutation <- function(n, R, S, D1, D2, levels1, levels2, mp, ties.method) {
  # Venkatraman (2000), section 2.2: the reference distribution permutes the
  # values transformed by the mixture CDF of their curve, pooled within each
  # disease group. Raw ranks of curves of different sizes or case mixes are
  # not on a common scale.
  R <- venkatraman.mixture.position(R, D1, levels1, mp)
  S <- venkatraman.mixture.position(S, D2, levels2, mp)

  R.controls <- R[D1 %in% levels1[1]]
  R.cases <- R[D1 %in% levels1[2]]
  S.controls <- S[D2 %in% levels2[1]]
  S.cases <- S[D2 %in% levels2[2]]

  # Permutation
  controls <- sample(c(R.controls, S.controls))
  cases <- sample(c(R.cases, S.cases))
  R[D1 %in% levels1[1]] <- controls[1:length(R.controls)]
  S[D2 %in% levels2[1]] <- controls[(length(R.controls) + 1):length(controls)]
  R[D1 %in% levels1[2]] <- cases[1:length(R.cases)]
  S[D2 %in% levels2[2]] <- cases[(length(R.cases) + 1):length(cases)]

  return(venkatraman.unpaired.stat(venkatraman.rank(R, ties.method), venkatraman.rank(S, ties.method), D1, D2, levels1, levels2, mp))
}

# Position of each observation on the mixture CDF of its curve,
# M(x) = mp * F_cases(x) + (1 - mp) * F_controls(x) (Venkatraman 2000, section
# 2.2). The paper assumes continuous data and uses M itself. Here each
# observation is put at a random point within the step of M at its rank (tied
# ranks share the step of the whole group): a random tie-break like the J(.)
# of Venkatraman and Begg (1996, p. 838), not part of the 2000 paper.
venkatraman.mixture.position <- function(R, D, levels, mp) {
  cases <- sort(R[D %in% levels[2]])
  controls <- sort(R[D %in% levels[1]])
  below <- mp * findInterval(R, cases, left.open = TRUE) / length(cases) +
    (1 - mp) * findInterval(R, controls, left.open = TRUE) / length(controls)
  upto <- mp * findInterval(R, cases) / length(cases) +
    (1 - mp) * findInterval(R, controls) / length(controls)
  return(below + runif(length(R)) * (upto - below))
}

# rank(ties.method = "random") draws random numbers even when nothing is tied:
# rank untied values without them, so that the random stream (and the result
# for a given seed) is the one of ties.method = "first"
venkatraman.rank <- function(x, ties.method) {
  if (identical(ties.method, "random") && !anyDuplicated(x)) {
    ties.method <- "first"
  }
  rank(x, ties.method = ties.method)
}

venkatraman.paired.stat <- function(R, S, D, levels) {
  R.controls <- R[D %in% levels[1]]
  R.cases <- R[D %in% levels[2]]
  S.controls <- S[D %in% levels[1]]
  S.cases <- S[D %in% levels[2]]
  n <- length(D)

  R.fn <- sapply(1:n, function(x) sum(R.cases <= x))
  R.fp <- sapply(1:n, function(x) sum(R.controls > x))
  S.fn <- sapply(1:n, function(x) sum(S.cases <= x))
  S.fp <- sapply(1:n, function(x) sum(S.controls > x))

  return(sum(abs((S.fn + S.fp) - (R.fn + R.fp))))
}

venkatraman.unpaired.stat <- function(R, S, D1, D2, levels1, levels2, mp) {
  R.controls <- R[D1 %in% levels1[1]]
  R.cases <- R[D1 %in% levels1[2]]
  S.controls <- S[D2 %in% levels2[1]]
  S.cases <- S[D2 %in% levels2[2]]
  n <- length(D1)
  m <- length(D2)

  R.fx <- sapply(1:n, function(x) sum(R.cases <= x)) / length(R.cases)
  R.gx <- sapply(1:n, function(x) sum(R.controls <= x)) / length(R.controls)
  S.fx <- sapply(1:m, function(x) sum(S.cases <= x)) / length(S.cases)
  S.gx <- sapply(1:m, function(x) sum(S.controls <= x)) / length(S.controls)
  R.p <- mp * R.fx + (1 - mp) * R.gx
  S.p <- mp * S.fx + (1 - mp) * S.gx
  R.exp <- mp * R.fx + (1 - mp) * (1 - R.gx)
  S.exp <- mp * S.fx + (1 - mp) * (1 - S.gx)

  # Venkatraman (2000), section 2.1: the error curves are estimated at the
  # distinct marker values, start at (p, e) = (0, 1 - mp) (no observation
  # below the threshold), end at (1, mp), and are joined by straight lines.
  # Repeated p values (tied ranks) are the same point.
  R.p <- c(0, R.p)
  R.exp <- c(1 - mp, R.exp)
  S.p <- c(0, S.p)
  S.exp <- c(1 - mp, S.exp)
  R.keep <- !duplicated(R.p)
  S.keep <- !duplicated(S.p)

  # The statistic is the area between them, the integral over [0, 1] of
  # |e_x - e_y| dp. The difference is linear between the points of both curves
  x <- sort(unique(c(R.p, S.p)))
  # rule = 2: the last p of each curve is 1 up to floating-point error
  d <- approx(R.p[R.keep], R.exp[R.keep], x, rule = 2)$y - approx(S.p[S.keep], S.exp[S.keep], x, rule = 2)$y
  dx <- diff(x)
  d0 <- d[-length(d)]
  d1 <- d[-1]
  # If the curves cross inside an interval, integrate both sides of the crossing
  crossing <- d0 * d1 < 0
  integral <- sum(ifelse(crossing,
    (d0^2 + d1^2) / (abs(d0) + abs(d1)) / 2,
    abs(d0 + d1) / 2
  ) * dx)
  return(integral)
}
