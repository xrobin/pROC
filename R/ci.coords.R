# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2010, 2011 Xavier Robin, Alexandre Hainard, Natacha Turck,
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

ci.coords <- function(...) {
  UseMethod("ci.coords")
}

ci.coords.formula <- function(formula, data, ...) {
  ci.coords(roc.formula(formula, data, ci=FALSE, ...), ...)
}

ci.coords.default <- function(response, predictor, ...) {
  ci.coords(roc.default(response, predictor, ci=FALSE, ...), ...)
}

ci.coords.smooth.roc <- function(smooth.roc,
											x, 
											input=c("specificity", "sensitivity"), ret=c("specificity", "sensitivity"), 
											best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5),
											conf.level = 0.95,
											boot.n = 2000,
											boot.stratified = TRUE,
											progress = getOption("pROCProgress")$name,
											...
                      ) {
  if (conf.level > 1 | conf.level < 0)
    stop("'conf.level' must be within the interval [0,1].")
  input <- match.arg(input)
  ret <- roc.utils.match.coords.ret.args(ret)
  if (is.character(x)) {
    x <- match.arg(x, c("all", "local maximas", "best"))
    if (x == "all" || x == "local maximas") {
      stop("'all' and 'local maximas' are not available for confidence intervals.")
    }
  }

  # Check if called with density.cases or density.controls
  if (is.null(smooth.roc$smoothing.args) || is.numeric(smooth.roc$smoothing.args$density.cases) || is.numeric(smooth.roc$smoothing.args$density.controls))
    stop("Cannot compute CI of ROC curves smoothed with numeric density.controls and density.cases.")

  # Get the non smoothed roc.
  roc <- attr(smooth.roc, "roc")
  roc$ci <- NULL # remove potential ci in roc to avoid infinite loop with smooth.roc()

  # prepare the calls
  smooth.roc.call <- as.call(c(match.fun("smooth.roc"), smooth.roc$smoothing.args))

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="Coords confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- raply(boot.n, stratified.ci.smooth.coords(roc, x, input, ret, best.method, best.weights, smooth.roc.call), .progress=progress)
  }
  else {
    perfs <- raply(boot.n, nonstratified.ci.smooth.coords(roc, x, input, ret, best.method, best.weights,smooth.roc.call), .progress=progress)
  }

  if (any(is.na(perfs))) {
    warning("NA value(s) produced during bootstrap were ignored.")
    perfs <- perfs[!apply(perfs, 1, function(x) any(is.na(x))),]
  }

  if (length(x) > 1) {
    inputs <- paste(input, x)
    rownames.ret <- paste(rep(inputs, each=length(ret)), ret, sep=": ")
  }
  else {
    rownames.ret <- ret
  }

  ci <- t(apply(perfs, 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
  rownames(ci) <- rownames.ret
  
  class(ci) <- c("ci.coords", "ci", class(ci))
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "roc") <- smooth.roc
  return(ci)
}

ci.coords.roc <- function(roc,
								  x,
								  input=c("threshold", "specificity", "sensitivity"), ret=c("threshold", "specificity", "sensitivity"),
								  best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5),
								  conf.level = 0.95,
								  boot.n = 2000,
								  boot.stratified = TRUE,
								  progress = getOption("pROCProgress")$name,
                      ...
                      ) {
  if (conf.level > 1 | conf.level < 0)
    stop("'conf.level' must be within the interval [0,1].")
  input <- match.arg(input)
  ret <- roc.utils.match.coords.ret.args(ret)
  if (is.character(x)) {
    x <- match.arg(x, c("all", "local maximas", "best"))
    if (x == "all" || x == "local maximas") {
      stop("'all' and 'local maximas' are not available for confidence intervals.")
    }
  }

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="Coords confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- raply(boot.n, stratified.ci.coords(roc, x, input, ret, best.method, best.weights), .progress=progress, .drop=FALSE)
  }
  else {
    perfs <- raply(boot.n, nonstratified.ci.coords(roc, x, input, ret, best.method, best.weights), .progress=progress, .drop=FALSE)
  }

  if (any(is.na(perfs))) {
    warning("NA value(s) produced during bootstrap were ignored.")
    perfs <- perfs[!apply(perfs, 1, function(x) any(is.na(x))),]
  }

  if (length(x) > 1) {
    inputs <- paste(input, x)
    rownames.ret <- paste(rep(inputs, each=length(ret)), ret, sep=": ")
  }
  else {
    rownames.ret <- ret
  }

  ci <- t(apply(perfs, 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
  rownames(ci) <- rownames.ret

  class(ci) <- c("ci.coords", "ci", class(ci))
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "roc") <- roc
  return(ci)
}
