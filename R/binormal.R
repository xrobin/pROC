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

# Maximum-likelihood fit of the conventional binormal model, Dorfman DD,
# Alf E (1969) J Math Psychol 6:487-496 (D&A), to the categories of the
# curve (as for the proper fit), started from the least-squares line on
# normal-deviate axes (Metz & Pan 1999, p. 19).
smooth_roc_binormal <- function(roc, n) {
  df <- data.frame(sp = qnorm(roc$specificities * ifelse(roc$percent, 1 / 100, 1)), se = qnorm(roc$sensitivities * ifelse(roc$percent, 1 / 100, 1)))
  df <- df[is.finite(df$sp) & is.finite(df$se), ]
  # Empty levels of ordered predictors give the same ROC point several times,
  # next to each other in curve order: count each point once in the fit
  df <- df[c(TRUE, diff(df$sp) != 0 | diff(df$se) != 0), ]
  # Curves with less than two points and degenerate data are refused as by
  # the least-squares smoothing of pROC 1.19. Degenerate data are fitted
  # exactly at infinite parameters by horizontal and vertical segments (Pan &
  # Metz 1997, Acad Radiol 4:380-389, p. 381).
  if (dim(df)[1] <= 1) { # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
  }
  if (smooth_roc_binormal_degenerate(df)) {
    stop("ROC curve not smoothable (degenerate data).")
  }
  categories <- proper_binormal_categories(roc)
  fit <- binormal_fit(categories$controls, categories$cases, binormal_ls(df))
  if (!fit$converged) {
    stop("ROC curve not smoothable (the binormal fit did not converge).")
  }
  # qnorm(TPF) = a + b * qnorm(FPF), with qnorm(FPF) = -qnorm(sp)
  se <- qnorm(seq(0, 1, 1 / (n - 1)))
  sp <- (fit$a - se) / fit$b

  return(list(
    sensitivities = pnorm(se) * ifelse(roc$percent, 100, 1),
    specificities = pnorm(sp) * ifelse(roc$percent, 100, 1),
    binormal = list(
      a = fit$a, b = fit$b, cutpoints = fit$cutpoints, logLik = fit$logLik,
      categories = data.frame(controls = categories$controls, cases = categories$cases)
    )
  ))
}

# The least-squares line through the points on normal-deviate axes,
# qnorm(sp) = c0 + c1 * qnorm(se): as qnorm(FPF) = -qnorm(sp), its a and b
# in qnorm(TPF) = a + b * qnorm(FPF) are a = -c0 / c1 and b = -1 / c1.
# The points of a ROC curve go down, so b > 0 unless they are degenerate.
binormal_ls <- function(df) {
  binormal_ls_params(lm(sp ~ se, df))
}
binormal_ls_params <- function(model) {
  co <- unname(coefficients(model))
  c(a = -co[1] / co[2], b = -1 / co[2])
}

# a and b of a curve smoothed with method = "binormal". Curves saved by
# pROC 1.19 and earlier hold the least-squares model instead of the fit.
# The $model branch exists only for them: drop it, with the test "binormal
# curves smoothed by pROC 1.19 keep their least-squares fit" in
# test-smooth.R, once they no longer need to be read.
smooth_roc_binormal_params <- function(smooth.roc) {
  if (is.null(smooth.roc$binormal)) {
    return(binormal_ls_params(smooth.roc$model))
  }
  c(a = smooth.roc$binormal$a, b = smooth.roc$binormal$b)
}

# Probabilities of the categories between the cut-points z (D&A Eqs. 1-2):
# F(Z_j) - F(Z_j-1) for the controls, F(b Z_j - a) - F(b Z_j-1 - a) for
# the cases, with Z_0 = -Inf and Z_I = Inf
binormal_probs <- function(a, b, z) {
  list(
    controls = proper_binormal_pdiff(-c(-Inf, z, Inf)),
    cases = proper_binormal_pdiff(-c(-Inf, b * z - a, Inf))
  )
}

# Maximum of the log-likelihood (D&A Eq. 3: sum of r_ij log P_ij) over a,
# b and the cut-points, from start = c(a, b), by the method of scoring
# that D&A use (p. 490, Eq. 7), finished with Newton steps
# (binormal_scoring)
binormal_fit <- function(k, l, start) {
  I <- length(k)
  a <- start[["a"]]
  b <- start[["b"]]
  # Initial cut-points: where the expected cumulative proportion of all the
  # observations under the start curve reaches the observed one. (D&A
  # p. 490 take the controls' cumulative proportions, which the runs of
  # cases would leave tied.)
  pooled <- cumsum(k + l)[-I] / sum(k, l)
  z <- vapply(pooled, function(p) {
    uniroot(function(z) (sum(k) * pnorm(z) + sum(l) * pnorm(b * z - a)) / sum(k, l) - p,
      c(-10, 10),
      extendInt = "upX", tol = 1e-12
    )$root
  }, numeric(1))
  binormal_scoring(c(a, b, z), k, l)
}

# Log-likelihood (D&A Eq. 3) at theta = c(a, b, cut-points); -Inf outside
# the parameter space (b <= 0 or cut-points not increasing)
binormal_loglik <- function(theta, k, l) {
  z <- theta[-(1:2)]
  if (theta[2] <= 0 || any(diff(z) <= 0)) {
    return(-Inf)
  }
  p <- binormal_probs(theta[1], theta[2], z)
  sum(k[k > 0] * log(p$controls[k > 0])) + sum(l[l > 0] * log(p$cases[l > 0]))
}

# Score (D&A Eqs. 4-6) and information at theta = c(a, b, cut-points).
# With observed = FALSE, the information is D&A's A, minus the expected
# second derivatives (p. 490): for each multinomial, the sum over the
# categories of n p' p'^T / p, with p' the derivatives of a category's
# probability. With observed = TRUE, it is minus the second derivatives of
# the log-likelihood: the sum of r (p' p'^T / p^2 - p'' / p). Category j
# depends on the cut-points j - 1 (lo) and j (hi) only, so both are
# tridiagonal in the cut-points, bordered by the rows of a and b: they are
# returned as the diagonal and off-diagonal in the cut-points (t.diag,
# t.off), the cut-points x (a, b) border, and the (a, b) block.
binormal_derivatives <- function(theta, k, l, observed = FALSE) {
  I <- length(k)
  a <- theta[1]
  b <- theta[2]
  z <- theta[-(1:2)]
  p <- binormal_probs(a, b, z)
  # values at the cut-points below (lo) and above (hi) each category; 0
  # beyond the ends
  phi <- c(0, dnorm(z), 0)
  f <- c(0, dnorm(b * z - a), 0)
  f1 <- -c(0, b * z - a, 0) * f # derivative of dnorm at b z - a
  zz <- c(0, z, 0)
  lo <- 1:I
  hi <- 2:(I + 1)
  zero <- numeric(I)
  # first (d) and second (dd) derivatives of each category's probability:
  # controls F(Z_j) - F(Z_j-1) (D&A Eq. 1), cases F(b Z_j - a) - F(b Z_j-1 - a)
  # (D&A Eq. 2)
  groups <- list(
    list(
      r = k, p = p$controls, n = sum(k),
      d = list(a = zero, b = zero, lo = -phi[lo], hi = phi[hi]),
      dd = list(aa = zero, ab = zero, alo = zero, ahi = zero, bb = zero, blo = zero, bhi = zero, lolo = zz[lo] * phi[lo], hihi = -zz[hi] * phi[hi])
    ),
    list(
      r = l, p = p$cases, n = sum(l),
      d = list(a = f[lo] - f[hi], b = f[hi] * zz[hi] - f[lo] * zz[lo], lo = -b * f[lo], hi = b * f[hi]),
      dd = list(
        aa = f1[hi] - f1[lo], ab = f1[lo] * zz[lo] - f1[hi] * zz[hi], alo = b * f1[lo], ahi = -b * f1[hi],
        bb = f1[hi] * zz[hi]^2 - f1[lo] * zz[lo]^2, blo = -(b * f1[lo] * zz[lo] + f[lo]), bhi = b * f1[hi] * zz[hi] + f[hi],
        lolo = -b^2 * f1[lo], hihi = b^2 * f1[hi]
      )
    )
  )
  score <- numeric(I + 1)
  t.diag <- numeric(I - 1)
  t.off <- numeric(I - 2)
  border <- matrix(0, I - 1, 2)
  ab <- matrix(0, 2, 2)
  for (g in groups) {
    s <- ifelse(g$r > 0, g$r / g$p, 0)
    score <- score + c(sum(s * g$d$a), sum(s * g$d$b), (s * g$d$hi)[-I] + (s * g$d$lo)[-1])
    # information of each category between two parameters x and y
    if (observed) {
      w <- ifelse(g$r > 0, g$r / g$p^2, 0)
      info <- function(x, y, xy) w * g$d[[x]] * g$d[[y]] - s * if (is.null(xy)) 0 else g$dd[[xy]]
    } else {
      w <- ifelse(g$p > 0, g$n / g$p, 0)
      info <- function(x, y, xy) w * g$d[[x]] * g$d[[y]]
    }
    t.diag <- t.diag + info("hi", "hi", "hihi")[-I] + info("lo", "lo", "lolo")[-1]
    t.off <- t.off + info("lo", "hi", NULL)[-c(1, I)]
    border <- border + cbind(
      info("a", "hi", "ahi")[-I] + info("a", "lo", "alo")[-1],
      info("b", "hi", "bhi")[-I] + info("b", "lo", "blo")[-1]
    )
    ab.ab <- sum(info("a", "b", "ab"))
    ab <- ab + matrix(c(sum(info("a", "a", "aa")), ab.ab, ab.ab, sum(info("b", "b", "bb"))), 2)
  }
  list(score = score, t.diag = t.diag, t.off = t.off, border = border, ab = ab)
}

# Solve the tridiagonal system T x = y (T symmetric, as diagonal and
# off-diagonal) for each column of y, by Gaussian elimination without
# pivoting: T is positive definite
binormal_tridiagonal_solve <- function(t.diag, t.off, y) {
  y <- as.matrix(y)
  n <- length(t.diag)
  d <- t.diag
  for (i in seq_len(n - 1)) {
    m <- t.off[i] / d[i]
    d[i + 1] <- d[i + 1] - m * t.off[i]
    y[i + 1, ] <- y[i + 1, ] - m * y[i, ]
  }
  y[n, ] <- y[n, ] / d[n]
  for (i in rev(seq_len(n - 1))) {
    y[i, ] <- (y[i, ] - t.off[i] * y[i + 1, ]) / d[i]
  }
  y
}

# The step information^-1 score at theta, solved through the tridiagonal
# block of the information (Schur complement), in O(I); NULL if it fails
binormal_step <- function(theta, k, l, observed) {
  d <- binormal_derivatives(theta, k, l, observed)
  solved <- binormal_tridiagonal_solve(d$t.diag, d$t.off, cbind(d$score[-(1:2)], d$border))
  schur <- d$ab - crossprod(d$border, solved[, -1])
  step.ab <- tryCatch(solve(schur, d$score[1:2] - crossprod(d$border, solved[, 1])), error = function(e) NULL)
  if (is.null(step.ab)) {
    return(NULL)
  }
  step <- c(step.ab, solved[, 1] - solved[, -1] %*% step.ab)
  if (!all(is.finite(step))) {
    return(NULL)
  }
  list(step = step, gain = sum(d$score * step))
}

# Maximize from theta = c(a, b, cut-points). The method of scoring (D&A
# p. 490, Eq. 7), theta + A^-1 score halved until the log-likelihood does
# not decrease, approaches the maximum. It can stall close to it, where
# the expected information differs from the curvature of the
# log-likelihood: Newton steps with the observed information finish, until
# the step is below 1e-10 of the parameters, which then also bounds the
# error of the estimates.
binormal_scoring <- function(theta, k, l) {
  ll <- binormal_loglik(theta, k, l)
  # slack for the Newton steps, which may lower the log-likelihood by its
  # rounding error near the maximum
  tolerance <- 1e-10 * max(1, abs(ll))
  for (iteration in 1:200) {
    st <- binormal_step(theta, k, l, observed = FALSE)
    if (is.null(st) || st$gain <= 1e-6 * max(1, abs(ll))) {
      break # close to the maximum
    }
    for (halving in 0:30) {
      candidate <- theta + st$step / 2^halving
      ll.candidate <- binormal_loglik(candidate, k, l)
      if (ll.candidate >= ll) {
        break
      }
    }
    if (ll.candidate < ll) {
      break
    }
    theta <- candidate
    ll <- ll.candidate
  }
  for (iteration in 1:50) {
    st <- binormal_step(theta, k, l, observed = TRUE)
    if (is.null(st)) {
      break
    }
    for (halving in 0:30) {
      candidate <- theta + st$step / 2^halving
      ll.candidate <- binormal_loglik(candidate, k, l)
      if (ll.candidate >= ll - tolerance) {
        break
      }
    }
    if (!(ll.candidate >= ll - tolerance)) {
      break # not near a maximum
    }
    theta <- candidate
    ll <- ll.candidate
    if (all(abs(st$step) <= 1e-10 * pmax(1, abs(theta)))) {
      return(list(a = theta[1], b = theta[2], cutpoints = theta[-(1:2)], logLik = ll, converged = TRUE))
    }
  }
  list(a = theta[1], b = theta[2], cutpoints = theta[-(1:2)], logLik = ll, converged = FALSE)
}

# Degenerate data (Metz and Pan 1999, p. 20): the points share one
# sensitivity or one specificity, so the fitted line would be horizontal
# or vertical. The rates are exact, so unique() is an exact test.
smooth_roc_binormal_degenerate <- function(df) {
  length(unique(df$se)) == 1 || length(unique(df$sp)) == 1
}

# Maximum-likelihood fit of the "proper" binormal model (PBM). References:
# - Metz CE, Pan X (1999) J Math Psychol 43:1-33 (M&P)
# - Pesce LL, Metz CE (2007) Acad Radiol 14:814-829 (P&M)
smooth_roc_proper <- function(roc, n) {
  load.suggested.package("mvtnorm")
  categories <- proper_binormal_categories(roc)
  k <- categories$controls
  l <- categories$cases
  if (length(k) < 2) {
    stop("ROC curve not smoothable (not enough points).")
  }
  # Perfect performance. M&P p. 27: no curve parameters can be estimated from
  # perfectly separated data. More generally, when every operating point has
  # FPF = 0 or TPF = 1 (no control above it, or no case below it; this also
  # covers ties straddling a separation), the likelihood is maximal only in
  # the limit AUC = 1. As PROPROC (return code 3: AUC = 1, d_a = 0, c = 1),
  # return that limit with a warning: Eq. 36 gives A = 1 at d_a = 0, c = 1
  # (rho = 0, F = 1/4). There is no maximum, hence no logLik or cut-points.
  I <- length(k)
  if (all(rev(cumsum(rev(k)))[-1] == 0 | cumsum(l)[-I] == 0)) {
    warning("Perfect performance: the proper binormal fit is the limit of AUC = 1 (d_a = 0, c = 1).", call. = FALSE)
    return(list(
      sensitivities = rep(1, n) * ifelse(roc$percent, 100, 1),
      specificities = seq(0, 1, length.out = n) * ifelse(roc$percent, 100, 1),
      proper = list(
        da = 0, c = 1, cutpoints = rep(NA_real_, I - 1), logLik = NA_real_,
        categories = data.frame(controls = k, cases = l)
      )
    ))
  }
  fit <- proper_binormal_fit(k, l, proper_binormal_start(k, l))
  if (fit$convergence != 0) {
    stop("ROC curve not smoothable (the proper binormal fit did not converge).")
  }

  # The curve (M&P Eqs. 33-34), at cut-points where FPF and TPF alternately
  # take equally spaced values, as for the fitdistr smoothing
  p <- seq(0, 1, length.out = n + 2)[-c(1, n + 2)]
  odd <- seq_len(n) %% 2 == 1
  v <- c(
    proper_binormal_invert(proper_binormal_fpf, p[odd], fit$da, fit$c),
    proper_binormal_invert(proper_binormal_tpf, p[!odd], fit$da, fit$c)
  )
  v <- sort(v)
  return(list(
    sensitivities = proper_binormal_tpf(v, fit$da, fit$c) * ifelse(roc$percent, 100, 1),
    specificities = (1 - proper_binormal_fpf(v, fit$da, fit$c)) * ifelse(roc$percent, 100, 1),
    proper = list(
      da = fit$da, c = fit$c, cutpoints = fit$cutpoints, logLik = fit$logLik,
      categories = data.frame(controls = k, cases = l)
    )
  ))
}

# Categories of continuous data (M&P pp. 24-25, P&M p. 816): the
# rank-ordered data are split into runs of a single truth state, and tied
# control and case values form a category of their own. From the least to
# the most "positive" category.
proper_binormal_categories <- function(roc) {
  controls <- roc$controls
  cases <- roc$cases
  if (is.ordered(controls) || is.ordered(cases)) {
    controls <- as.integer(controls)
    cases <- as.integer(cases)
  }
  if (roc$direction == ">") {
    controls <- -controls
    cases <- -cases
  }
  values <- sort(unique(c(controls, cases)))
  k <- tabulate(match(controls, values), length(values))
  l <- tabulate(match(cases, values), length(values))
  # 0: tied controls and cases, 1: controls only, 2: cases only
  state <- ifelse(k > 0 & l > 0, 0L, ifelse(k > 0, 1L, 2L))
  run <- cumsum(c(TRUE, state[-1] != state[-length(state)] | state[-1] == 0L))
  list(controls = as.vector(rowsum(k, run)), cases = as.vector(rowsum(l, run)))
}

# FPF and TPF of the PBM (M&P p. 12, Eqs. 33 and 34) at the cut-points v
# (M&P Eq. 35: Inf and -Inf are allowed and give 0 and 1).
# The arguments of the Phi functions in Eqs. 33-34: A is the first one, B the
# second one, without H(c). B is not used for c = 0, where the curly brackets
# vanish (M&P p. 8).
proper_binormal_args <- function(v, da, c, slope, sign) {
  s <- da * sqrt(1 + c^2)
  list(A = -slope * v + sign * s / 2, B = if (c != 0) -slope * v + s / (2 * c))
}
proper_binormal_fpf <- function(v, da, c) {
  proper_binormal_rate(proper_binormal_args(v, da, c, 1 - c, -1), c)
}
proper_binormal_tpf <- function(v, da, c) {
  proper_binormal_rate(proper_binormal_args(v, da, c, 1 + c, 1), c)
}
proper_binormal_rate <- function(args, c) {
  if (c == 0) {
    pnorm(args$A)
  } else if (c < 0) {
    pnorm(args$A) + pnorm(args$B)
  } else { # Phi(B) - H(c), as an upper tail
    pnorm(args$A) - pnorm(args$B, lower.tail = FALSE)
  }
}
# Derivative of FPF or TPF with respect to v
proper_binormal_rate_dv <- function(args, c, slope) {
  d <- dnorm(args$A)
  if (c != 0) {
    d <- d + dnorm(args$B)
  }
  -slope * d
}
# Phi(x[i]) - Phi(x[i + 1]) for decreasing x, from the tail where it does
# not cancel
proper_binormal_pdiff <- function(x) {
  n <- length(x)
  upper <- x[-1] > 0
  lower.tail <- pnorm(x)
  upper.tail <- pnorm(x, lower.tail = FALSE)
  d <- lower.tail[-n] - lower.tail[-1]
  d[upper] <- upper.tail[-1][upper] - upper.tail[-n][upper]
  d
}
# Probabilities of the categories (M&P p. 9, Eqs. 24 and 25) between the
# cut-points 'v' (the bounds of the range included). H(c) cancels out in the
# differences.
proper_binormal_probs <- function(args) {
  p <- proper_binormal_pdiff(args$A)
  if (!is.null(args$B)) {
    p <- p + proper_binormal_pdiff(args$B)
  }
  p
}

# Range of the cut-points v (M&P p. 12, Eq. 35)
proper_binormal_range <- function(da, c) {
  bound <- da * sqrt(1 + c^2) / (4 * c)
  if (c < 0) {
    c(bound, Inf)
  } else if (c > 0) {
    c(-Inf, bound)
  } else {
    c(-Inf, Inf)
  }
}

# Log-likelihood (M&P p. 4, Eq. 2) and its gradient with respect to the
# cut-points v
proper_binormal_loglik <- function(v, da, c, k, l, gradient = FALSE) {
  range <- proper_binormal_range(da, c)
  vv <- c(range[1], v, range[2])
  fpf.args <- proper_binormal_args(vv, da, c, 1 - c, -1)
  tpf.args <- proper_binormal_args(vv, da, c, 1 + c, 1)
  p <- proper_binormal_probs(fpf.args)
  q <- proper_binormal_probs(tpf.args)
  if (any(p[k > 0] <= 0) || any(q[l > 0] <= 0)) {
    return(if (gradient) rep(0, length(v)) else -Inf)
  }
  if (!gradient) {
    return(sum(k[k > 0] * log(p[k > 0])) + sum(l[l > 0] * log(q[l > 0])))
  }
  kp <- k / p
  kp[k == 0] <- 0
  lq <- l / q
  lq[l == 0] <- 0
  I <- length(k)
  inner <- 2:I # v[j] is the upper bound of category j and the lower of j + 1
  dfpf <- proper_binormal_rate_dv(lapply(fpf.args, `[`, inner), c, 1 - c)
  dtpf <- proper_binormal_rate_dv(lapply(tpf.args, `[`, inner), c, 1 + c)
  (kp[-1] - kp[-I]) * dfpf + (lq[-1] - lq[-I]) * dtpf
}

# Initial cut-points (P&M p. 819): each cut-point is the MLE of the data
# dichotomized at it, found with Brent's method, and bounds the next one.
# The maximum is bracketed first: each binomial term of the dichotomized
# log-likelihood is unimodal in v, with its mode where FPF (or TPF) equals
# the observed fraction, so their sum has its maximum between the two modes.
proper_binormal_start_v <- function(da, c, k, l) {
  I <- length(k)
  interval <- proper_binormal_interval(da, c)
  lower <- interval[1]
  upper <- interval[2]
  xlogy <- function(x, y) if (x > 0) x * log(y) else 0
  dichotomized.loglik <- function(v, above.controls, above.cases) {
    fpf <- proper_binormal_fpf(v, da, c)
    tpf <- proper_binormal_tpf(v, da, c)
    ll <- xlogy(above.controls, fpf) + xlogy(sum(k) - above.controls, 1 - fpf) +
      xlogy(above.cases, tpf) + xlogy(sum(l) - above.cases, 1 - tpf)
    # optimize() needs finite values
    if (is.finite(ll)) ll else -.Machine$double.xmax
  }
  v <- numeric(I - 1)
  previous <- lower
  for (j in seq_len(I - 1)) {
    above.controls <- sum(k[-(1:j)])
    above.cases <- sum(l[-(1:j)])
    modes <- c(
      proper_binormal_invert(proper_binormal_fpf, above.controls / sum(k), da, c),
      proper_binormal_invert(proper_binormal_tpf, above.cases / sum(l), da, c)
    )
    bracket <- pmax(previous, range(modes))
    v[j] <- if (bracket[1] == bracket[2]) {
      bracket[1]
    } else {
      optimize(dichotomized.loglik, bracket,
        above.controls = above.controls, above.cases = above.cases,
        maximum = TRUE
      )$maximum
    }
    previous <- v[j]
  }
  # Strictly increasing and strictly inside the range of Eq. 35
  gap <- (upper - lower) * 1e-8
  v <- lower + gap * seq_len(I - 1) + cummax(v - lower)
  if (c > 0) {
    v <- pmin(v, upper - gap * rev(seq_len(I - 1)))
  }
  v
}

# Starting values of (da, c). P&M p. 818: from the least-squares
# conventional binormal fit, with da = 0 when |c| > 2 da (footnote, with c
# then from the relationship between FPF and TPF at da = 0)
proper_binormal_start <- function(k, l) {
  # The empirical operating points (P&M p. 818) at the boundaries of the
  # categories, inside the unit square, from the counts: the same whatever
  # the scale of the curve (percent)
  fpf <- rev(cumsum(rev(k)))[-1] / sum(k)
  tpf <- rev(cumsum(rev(l)))[-1] / sum(l)
  interior <- fpf > 0 & fpf < 1 & tpf > 0 & tpf < 1
  fpf <- fpf[interior]
  tpf <- tpf[interior]
  if (length(fpf) == 0) {
    return(NULL)
  }
  df <- data.frame(sp = qnorm(fpf, lower.tail = FALSE), se = qnorm(tpf))
  # The sign of c: from the least-squares fit; for degenerate data, from its
  # limit: a vertical line (one specificity) has b = Inf, so c = 1, and a
  # horizontal line (one sensitivity) b = 0, so c = -1. A single point has
  # no line: both signs are tried.
  signs <- c(-1, 1)
  if (length(unique(df$sp)) == 1 && length(unique(df$se)) > 1) {
    signs <- 1
  } else if (length(unique(df$se)) == 1 && length(unique(df$sp)) > 1) {
    signs <- -1
  } else if (dim(df)[1] > 1) {
    params <- binormal_ls(df)
    a <- params[["a"]]
    b <- params[["b"]]
    if (b > 0) {
      # M&P p. 9, Eq. 26, and p. 11, Eq. 32
      start <- c(da = sqrt(2) * a / sqrt(1 + b^2), c = (b - 1) / (b + 1))
      if (abs(start[["c"]]) <= 2 * start[["da"]]) {
        return(start)
      }
      signs <- sign(start[["c"]])
    }
  }
  # Eqs. 33-34 at da = 0: FPF = 2 Phi(-(1 - c) v) and TPF = 2 Phi(-(1 + c) v)
  # if c < 0 (v >= 0), and 1 - FPF and 1 - TPF in the same form if c > 0
  # (v <= 0). A least-squares line through the origin of qnorm(TPF / 2) on
  # qnorm(FPF / 2) (or of the same in 1 - TPF and 1 - FPF) has the slope
  # (1 + c) / (1 - c). The fit with c of the expected sign and the smallest
  # residuals is kept (a point on the negative diagonal has no slope: NaN).
  # NULL if there is none.
  best <- NULL
  rss <- Inf
  for (sign in signs) {
    x <- qnorm(if (sign < 0) fpf / 2 else (1 - fpf) / 2)
    y <- qnorm(if (sign < 0) tpf / 2 else (1 - tpf) / 2)
    b <- sum(x * y) / sum(x^2)
    c <- (b - 1) / (b + 1)
    if (isTRUE(sign * c > 0) && sum((y - b * x)^2) < rss) {
      best <- c(da = 0, c = max(-0.99, min(0.99, c)))
      rss <- sum((y - b * x)^2)
    }
  }
  best
}

# AUC of the PBM (M&P p. 12, Eq. 36):
# Phi(da / sqrt(2)) + 2 F(-da / sqrt(2), 0; rho), with rho = -(1 - c^2) / (1 + c^2),
# F being the bivariate normal distribution function, from mvtnorm (TVPACK:
# its algorithm for bivariate semi-infinite regions; the default algorithm is
# randomized). F = 0 at c = 0 (M&P p. 12), where rho = -1.
proper_binormal_auc <- function(da, c) {
  F <- 0
  if (c != 0) {
    load.suggested.package("mvtnorm")
    rho <- -(1 - c^2) / (1 + c^2)
    F <- mvtnorm::pmvnorm(
      upper = c(-da / sqrt(2), 0), corr = matrix(c(1, rho, rho, 1), 2),
      algorithm = mvtnorm::TVPACK()
    )[1]
  }
  pnorm(da / sqrt(2)) + 2 * F
}

# da on the isopleth AUC(da, c) = auc, or 0 if AUC(0, c) >= auc
proper_binormal_isopleth <- function(auc, c) {
  if (proper_binormal_auc(0, c) >= auc) {
    return(0)
  }
  upper <- 1
  while (proper_binormal_auc(upper, c) < auc) {
    if (upper >= 64) {
      return(upper)
    }
    upper <- upper * 2
  }
  uniroot(function(da) proper_binormal_auc(da, c) - auc, c(0, upper), tol = 1e-10)$root
}

# Maximum-likelihood fit. P&M p. 819: all the parameters (da, c and the
# cut-points v) are iterated together (with optim's BFGS here, in place of
# their trust-region Newton method); a step that leaves d_a >= 0, |c| < 1 or
# the range of Eq. 35, or makes the cut-points cross, is rejected (P&M
# pp. 819-820): the objective is Inf there, and BFGS backtracks. da and c are
# searched as |x| and tanh(y). Starts: the least-squares estimate
# (proper_binormal_start, NULL if there is none; P&M p. 818), and a
# symmetrical curve (c = 0) and two curves skewed in opposite directions,
# with the same AUC (P&M p. 819). These are on the isopleth of the empirical
# (trapezoidal) AUC, with the skewed ones chosen as M&P do (pp. 19-20):
# among 20 points equally spaced (in c) on the isopleth, the one of each
# sign with the largest log-likelihood maximized over the cut-points. If the
# starts reach different maxima, 5 more are spread on the isopleth of the
# best one, between its smallest and largest c.
proper_binormal_fit <- function(k, l, start) {
  starts <- c(if (!is.null(start)) list(start), proper_binormal_isopleth_starts(k, l))
  fits <- lapply(unique(starts), proper_binormal_fit_from, k = k, l = l)
  loglik <- vapply(fits, `[[`, numeric(1), "logLik")
  reached <- is.finite(loglik)
  if (sum(reached) > 1 && diff(range(loglik[reached])) > 1e-6 * max(1, abs(max(loglik)))) {
    best <- fits[[which.max(loglik)]]
    cs <- vapply(fits[reached], `[[`, numeric(1), "c")
    cs <- seq(min(cs), max(cs), length.out = 7)[2:6]
    auc <- proper_binormal_auc(best$da, best$c)
    more <- lapply(cs, function(c) c(da = proper_binormal_isopleth(auc, c), c = c))
    fits <- c(fits, lapply(more, proper_binormal_fit_from, k = k, l = l))
    loglik <- vapply(fits, `[[`, numeric(1), "logLik")
  }
  fits[[which.max(loglik)]]
}

proper_binormal_isopleth_starts <- function(k, l) {
  # empirical AUC from the counts: P(case > control) + P(tie) / 2
  cases.above <- rev(cumsum(rev(l))) - l
  auc <- sum(k * (cases.above + l / 2)) / (sum(k) * sum(l))
  # At AUC <= 0.5 the isopleth is the chance line alone, a stationary point
  # of the search: the grid is taken on the isopleth of 1 - AUC instead, and
  # at least on that of the smallest AUC above 0.5 the counts allow (no
  # source covers these data)
  auc <- max(auc, 1 - auc, 0.5 + 1 / (2 * sum(k) * sum(l)))
  # the isopleth reaches da = 0 at +-c.max (the AUC grows with |c| at da = 0)
  c.max <- uniroot(function(c) proper_binormal_auc(0, c) - auc, c(0, 1), tol = 1e-10)$root
  grid <- c.max * seq(-1, 1, length.out = 22)[2:21]
  loglik <- vapply(grid, function(c) {
    proper_binormal_profile(proper_binormal_isopleth(auc, c), c, k, l)
  }, numeric(1))
  skewed <- lapply(c(-1, 1), function(sign) {
    c <- grid[sign * grid > 0][which.max(loglik[sign * grid > 0])]
    c(da = proper_binormal_isopleth(auc, c), c = c)
  })
  c(list(c(da = proper_binormal_isopleth(auc, 0), c = 0)), skewed)
}

# Log-likelihood maximized over the cut-points, at fixed (da, c) (M&P p. 20)
proper_binormal_profile <- function(da, c, k, l) {
  I <- length(k)
  range <- proper_binormal_range(da, c)
  minus.loglik <- function(v) {
    if (v[1] <= range[1] || v[I - 1] >= range[2] || any(diff(v) <= 0)) {
      return(Inf)
    }
    -proper_binormal_loglik(v, da, c, k, l)
  }
  gradient <- function(v) -proper_binormal_loglik(v, da, c, k, l, gradient = TRUE)
  v <- proper_binormal_start_v(da, c, k, l)
  if (!is.finite(minus.loglik(v))) {
    return(-Inf)
  }
  -optim(v, minus.loglik, gradient, method = "BFGS", control = list(maxit = 1000, reltol = 1e-10))$value
}

proper_binormal_fit_from <- function(start, k, l) {
  I <- length(k)
  minus.loglik <- function(par) {
    da <- abs(par[1])
    c <- tanh(par[2])
    v <- par[-(1:2)]
    range <- proper_binormal_range(da, c)
    if (abs(c) >= 1 || v[1] <= range[1] || v[I - 1] >= range[2] || any(diff(v) <= 0)) {
      return(Inf)
    }
    -proper_binormal_loglik(v, da, c, k, l)
  }
  gradient <- function(par) {
    # analytic for the cut-points; central differences for (x, y), one-sided
    # next to a rejected point
    g <- -proper_binormal_loglik(par[-(1:2)], abs(par[1]), tanh(par[2]), k, l, gradient = TRUE)
    f <- minus.loglik(par)
    h <- 1e-6
    numerical <- vapply(1:2, function(i) {
      up <- par
      up[i] <- up[i] + h
      down <- par
      down[i] <- down[i] - h
      f.up <- minus.loglik(up)
      # at x = 0, da = |x| is symmetric: the central difference is 0 and the
      # search could never leave da = 0
      f.down <- if (i == 1 && par[1] == 0) Inf else minus.loglik(down)
      if (is.finite(f.up) && is.finite(f.down)) {
        (f.up - f.down) / (2 * h)
      } else if (is.finite(f.up)) {
        (f.up - f) / h
      } else if (is.finite(f.down)) {
        (f - f.down) / h
      } else {
        0
      }
    }, numeric(1))
    c(numerical, g)
  }
  par <- c(start[["da"]], atanh(start[["c"]]), proper_binormal_start_v(start[["da"]], start[["c"]], k, l))
  if (!is.finite(minus.loglik(par))) { # the start does not fit the data
    return(list(da = start[["da"]], c = start[["c"]], logLik = -Inf, cutpoints = NULL, convergence = 1))
  }
  o <- optim(par, minus.loglik, gradient, method = "BFGS", control = list(maxit = 5000, reltol = 1e-12))
  list(
    da = abs(o$par[1]), c = tanh(o$par[2]), logLik = -o$value,
    cutpoints = o$par[-(1:2)], convergence = o$convergence
  )
}

# Cut-points v where 'rate' (proper_binormal_fpf or _tpf, decreasing in v)
# equals p, within proper_binormal_interval (its bounds for p beyond the
# rates there)
proper_binormal_invert <- function(rate, p, da, c) {
  interval <- proper_binormal_interval(da, c)
  bounds <- rate(interval, da, c)
  vapply(p, function(p) {
    if (p >= bounds[1]) {
      interval[1]
    } else if (p <= bounds[2]) {
      interval[2]
    } else {
      uniroot(function(v) rate(v, da, c) - p, interval, tol = 1e-12)$root
    }
  }, numeric(1))
}

# The range of Eq. 35 within [-far, far], to search the cut-points: beyond
# far, all the arguments of the Phi functions of FPF and TPF in Eqs. 33-34
# are 40 or more away from 0 (with the sign of the limit), so FPF and TPF
# are 0 or 1 to machine precision. The bound of Eq. 35 itself can be very
# far as c tends to 0.
proper_binormal_interval <- function(da, c) {
  range <- proper_binormal_range(da, c)
  far <- (40 + da * sqrt(1 + c^2)) / (1 - abs(c))
  c(max(range[1], -far), min(range[2], far))
}
