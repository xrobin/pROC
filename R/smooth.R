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

smooth <- function(...)
  UseMethod("smooth")

smooth.default <- function(...) {
  stats::smooth(...)
}

smooth.smooth.roc <- function(smooth.roc, ...) {
  roc <- attr(smooth.roc, "roc")
  if (is.null(roc))
    stop("Cannot smooth a ROC curve generated directly with numeric 'density.controls' and 'density.cases'.")
  smooth.roc(roc, ...)
}

smooth.roc <- function(roc, method = c("binormal", "density", "fitdistr", "logcondens", "logcondens.smooth"), n = 512, bw = "nrd0",
                       density = NULL, density.controls = density, density.cases = density,
                       start = NULL, start.controls = start, start.cases = start,
                       reuse.auc = TRUE, reuse.ci = FALSE,
                       ...) {
  method <- match.arg(method)
  
  if (is.ordered(roc$original.predictor) && (method == "density" || method =="fitidstr"))
    stop("ROC curves of ordered predictors can be smoothed only with binormal smoothing.")

  if (method == "binormal") {
    sesp <- smooth.roc.binormal(roc, n)
  }
  else if (method == "fitdistr") {
    sesp <- smooth.roc.fitdistr(roc, n, density.controls, density.cases, start.controls, start.cases, ...)
  }
  else if (method == "density") {
    sesp <- smooth.roc.density(roc, n, density.controls, density.cases, bw, ...)
  }
  else if (method == "logcondens") {
    sesp <- smooth.roc.logcondens(roc, n)
  }
  else if (method == "logcondens.smooth") {
    sesp <- smooth.roc.logcondens.smooth(roc, n)
  }
  else {
    stop(sprintf("Impossible smooth method value '%s'. Please report this bug to <%s>.",
                 method, utils::packageDescription("pROC")$BugReports))
  }

  class(sesp) <- "smooth.roc"
  sesp <- sort(sesp) # sort se and sp
  # anchor SE/SP at 0/100
  sesp$specificities <- c(0, as.vector(sesp$specificities), ifelse(roc$percent, 100, 1))
  sesp$sensitivities <- c(ifelse(roc$percent, 100, 1), as.vector(sesp$sensitivities), 0)
  attr(sesp, "roc") <- roc # keep the original roc. May be useful in CI.
  sesp$percent <- roc$percent # keep some basic roc specifications
  sesp$direction <- roc$direction
  sesp$call <- match.call()
  # keep smoothing arguments (for print and bootstrap)
  sesp$smoothing.args <- list(...)
  sesp$smoothing.args$method <- method
  sesp$smoothing.args$n <- n
  sesp$smoothing.args$start.controls <- start.controls
  sesp$smoothing.args$start.cases <- start.cases
  sesp$smoothing.args$density.controls <- density.controls
  sesp$smoothing.args$density.cases <- density.cases
  sesp$smoothing.args$bw <- bw
  # complete fit.controls/cases if a function was passed as densfun
  if (method == "fitdistr") {
    if (is.null(sesp$fit.controls$densfun)) {
      if (missing(density.controls))
        sesp$fit.controls$densfun <- deparse(substitute(density))
      else
        sesp$fit.controls$densfun <- deparse(substitute(density.controls))
    }
    if (is.null(sesp$fit.cases$densfun)) {
      if (missing(density.cases))
        sesp$fit.cases$densfun <- deparse(substitute(density))
      else
        sesp$fit.cases$densfun <- deparse(substitute(density.cases))
    }
  }

  # if there was an auc and a ci, re-do them
  if (!is.null(roc$auc) && reuse.auc) {
    args <- attributes(roc$auc)
    args$roc <- NULL
    args$smooth.roc <- sesp
    sesp$auc <- do.call("auc.smooth.roc", args)
  }
  if (!is.null(roc$ci) && reuse.ci){
    args <- attributes(roc$ci)
    args$roc <- NULL
    args$smooth.roc <- sesp
    sesp$ci <- do.call(paste(class(roc$ci), "smooth.roc", sep="."), args)
  }

  return(sesp)
}

smooth.roc.density <- function(roc, n, density.controls, density.cases, bw,
                               # catch args for density
                               cut = 3, adjust = 1, kernel = window, window = "gaussian",
                               percent = roc$percent, direction = roc$direction,
                               ...) {
  if (!is.numeric(density.controls) || !is.numeric(density.cases)) {
  	predictor <- c(roc$controls, roc$cases)
    if (is.character(bw)) {
    	bw <- match.fun(paste("bw", bw, sep="."))(predictor)
    }
    bw <- bw * adjust
    from <- min(predictor) - (cut * bw)
    to <- max(predictor) + (cut * bw)
  }
  if (mode(density.controls) == "function") {
    density.controls <- density.controls(roc$controls, n=n, from=from, to=to, bw=bw, kernel=kernel, ...)
    if (! is.numeric(density.controls)) {
      if (is.list(density.controls) && !is.null(density.controls$y) && is.numeric(density.controls$y))
        density.controls <- density.controls$y
      else
        stop("The 'density' function must return a numeric vector or a list with a 'y' item.")
    }
  }
  else if (is.null(density.controls))
    density.controls <- suppressWarnings(density(roc$controls, n=n, from=from, to=to, bw=bw, kernel=kernel, ...))$y
  else if (! is.numeric(density.controls))
    stop("'density.controls' must be either NULL, a function or numeric values of density (over the y axis).")
  if (mode(density.cases) == "function") {
    density.cases <- density.cases(roc$cases, n=n, from=from, to=to, bw=bw, kernel=kernel, ...)
    if (! is.numeric(density.cases)) {
      if (is.list(density.cases) && !is.null(density.cases$y) && is.numeric(density.cases$y))
        density.cases <- density.cases$y
      else
        stop("The 'density' function must return a numeric vector or a list with a 'y' item.")
    }
  }
  else if (is.null(density.cases))
    density.cases <- suppressWarnings(density(roc$cases, n=n, from=from, to=to, bw=bw, kernel=kernel, ...))$y
  else if (! is.numeric(density.cases))
    stop("'density.cases' must be either NULL, a function or numeric values of density (over the y axis).")
  if (length(density.controls) != length(density.cases))
    stop("Length of 'density.controls' and 'density.cases' differ.")

  perfs <- sapply((1:length(density.controls))+.5, roc.utils.perfs.dens, x=(1:length(density.controls))+.5, dens.controls=density.controls, dens.cases=density.cases, direction=direction)

  return(list(sensitivities = perfs[2,] * ifelse(percent, 100, 1),
              specificities = perfs[1,] * ifelse(percent, 100, 1)))
}

smooth.roc.binormal <- function(roc, n) {
  df <- data.frame(sp=qnorm(roc$sp * ifelse(roc$percent, 1/100, 1)), se=qnorm(roc$se * ifelse(roc$percent, 1/100, 1)))
  df <- df[apply(df, 1, function(x) all(is.finite(x))),]
  if (dim(df)[1] <= 1) # ROC curve or with only 1 point
    stop("ROC curve not smoothable (not enough points).")
  model <- lm(sp~se, df)
  if(any(is.na(model$coefficients[2])))
    stop("ROC curve not smoothable (not enough points).")
  se <- qnorm(seq(0, 1, 1/(n-1)))
  sp <- predict(model, data.frame(se))

  return(list(sensitivities = pnorm(se) * ifelse(roc$percent, 100, 1),
              specificities = pnorm(sp) * ifelse(roc$percent, 100, 1),
              model = model))
}

smooth.roc.fitdistr <- function(roc, n, densfun.controls, densfun.cases, start.controls, start.cases, ...) {
  load.suggested.package("MASS")

  densfuns.list <- list(beta = "dbeta", cauchy = "dcauchy", "chi-squared" = "dchisq", exponential = "dexp", f = "df",
                        gamma = "dgamma", geometric = "dgeom", "log-normal" = "dlnorm", lognormal = "dlnorm",
                        logistic = "dlogis", "negative binomial" = "dnbinom", normal = "dnorm", poisson = "dpois",
                        t = "dt", weibull = "dweibull")

  if (is.null(densfun.controls))
    densfun.controls <- "normal"
  else if (is.character(densfun.controls))
    densfun.controls <- match.arg(densfun.controls, names(densfuns.list))

  if (is.null(densfun.cases))
    densfun.cases <- "normal"              
  else if (is.character(densfun.cases))
    densfun.cases <- match.arg(densfun.cases, names(densfuns.list))

  fit.controls <- MASS::fitdistr(roc$controls, densfun.controls, start.controls, ...)
  fit.cases <- MASS::fitdistr(roc$cases, densfun.cases, start.cases, ...)

  # store function name in fitting results
  if (mode(densfun.controls) != "function")
    fit.controls$densfun <- densfun.controls
  if (mode(densfun.cases) != "function")
    fit.cases$densfun <- densfun.cases

  x <- seq(min(c(roc$controls, roc$cases)), max(c(roc$controls, roc$cases)), length.out=n)

  # get the actual function name for do.call
  if (is.character(densfun.controls))
    densfun.controls <- match.fun(densfuns.list[[densfun.controls]])
  if (is.character(densfun.cases))
    densfun.cases <- match.fun(densfuns.list[[densfun.cases]])

  # ... that should be passed to densfun
  if (length(list(...)) > 0 && any(names(list(...)) %in% names(formals(densfun.controls))))
    dots.controls <- list(...)[names(formals(densfun.controls))[match(names(list(...)), names(formals(densfun.controls)))]]
  else
    dots.controls <- list()
  if (length(list(...)) > 0 && any(names(list(...)) %in% names(formals(densfun.cases))))
    dots.cases <- list(...)[names(formals(densfun.cases))[match(names(list(...)), names(formals(densfun.cases)))]]
  else
    dots.cases <- list()
   
  density.controls <- do.call(densfun.controls, c(list(x=x), fit.controls$estimate, dots.controls))
  density.cases <- do.call(densfun.cases, c(list(x=x), fit.cases$estimate, dots.cases))

  perfs <- sapply(x, roc.utils.perfs.dens, x=x, dens.controls=density.controls, dens.cases=density.cases, direction=roc$direction)

  return(list(sensitivities = perfs[2,] * ifelse(roc$percent, 100, 1),
              specificities = perfs[1,] * ifelse(roc$percent, 100, 1),
              fit.controls=fit.controls, fit.cases=fit.cases))
}

smooth.roc.logcondens <- function(roc, n) {
  load.suggested.package("logcondens")

  sp <- seq(0, 1, 1/(n-1))
  logcondens <- logcondens::logConROC(roc$cases, roc$controls, sp)
  se <- logcondens$fROC

  return(list(sensitivities = se * ifelse(roc$percent, 100, 1),
              specificities = (1 - sp) * ifelse(roc$percent, 100, 1),
              logcondens = logcondens))
}

smooth.roc.logcondens.smooth <- function(roc, n) {
  load.suggested.package("logcondens")

  sp <- seq(0, 1, 1/(n-1))
  logcondens <- logcondens::logConROC(roc$cases, roc$controls, sp)
  se <- logcondens$fROC.smooth

  return(list(sensitivities = se * ifelse(roc$percent, 100, 1),
              specificities = (1 - sp) * ifelse(roc$percent, 100, 1),
              logcondens = logcondens))
}
