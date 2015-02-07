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

ci.se <- function(...) {
  UseMethod("ci.se")
}

ci.se.formula <- function(formula, data, ...) {
  ci.se(roc.formula(formula, data, ci=FALSE, ...), ...)
}

ci.se.default <- function(response, predictor, ...) {
  ci.se(roc.default(response, predictor, ci=FALSE, ...), ...)
}

ci.se.smooth.roc <- function(smooth.roc,
                      specificities = seq(0, 1, .1) * ifelse(smooth.roc$percent, 100, 1),
                      conf.level = 0.95,
                      boot.n = 2000,
                      boot.stratified = TRUE,
                      progress = getOption("pROCProgress")$name,
                      parallel = FALSE,
                      ...
                      ) {
  if (conf.level > 1 | conf.level < 0)
    stop("'conf.level' must be within the interval [0,1].")
  
  if (roc.utils.is.perfect.curve(smooth.roc)) {
  	warning("ci.se() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }

  # Check if called with density.cases or density.controls
  if (is.null(smooth.roc$smoothing.args) || is.numeric(smooth.roc$smoothing.args$density.cases) || is.numeric(smooth.roc$smoothing.args$density.controls))
    stop("Cannot compute CI of ROC curves smoothed with numeric density.controls and density.cases.")

  # Get the non smoothed roc.
  roc <- attr(smooth.roc, "roc")
  roc$ci <- NULL # remove potential ci in roc to avoid infinite loop with smooth.roc()

  # prepare the calls
  smooth.roc.call <- as.call(c(getS3method("smooth", "roc"), smooth.roc$smoothing.args))

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="SE confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- ldply(1:boot.n, stratified.ci.smooth.se, roc=roc, sp=specificities, smooth.roc.call=smooth.roc.call, .progress=progress, .parallel=parallel)
  }
  else {
    perfs <- ldply(1:boot.n, nonstratified.ci.smooth.se, roc=roc, sp=specificities, smooth.roc.call=smooth.roc.call, .progress=progress, .parallel=parallel)
  }

  if (any(is.na(perfs))) {
    warning("NA value(s) produced during bootstrap were ignored.")
    perfs <- perfs[!apply(perfs, 1, function(x) any(is.na(x))),]
  }

  ci <- t(apply(perfs, 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
  rownames(ci) <- paste(specificities, ifelse(roc$percent, "%", ""), sep="")

  class(ci) <- c("ci.se", "ci", class(ci))
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "specificities") <- specificities
  attr(ci, "roc") <- smooth.roc
  return(ci)
}

ci.se.roc <- function(roc,
                      specificities = seq(0, 1, .1) * ifelse(roc$percent, 100, 1),
                      conf.level = 0.95,
                      boot.n = 2000,
                      boot.stratified = TRUE,
                      progress = getOption("pROCProgress")$name,
                      parallel = FALSE,
                      ...
                      ) {
  if (conf.level > 1 | conf.level < 0)
    stop("'conf.level' must be within the interval [0,1].")
  
  if (roc.utils.is.perfect.curve(roc)) {
  	warning("ci.se() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="SE confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- ldply(1:boot.n, stratified.ci.se, roc=roc, sp=specificities, .progress=progress, .parallel=parallel)
  }
  else {
    perfs <- ldply(1:boot.n, nonstratified.ci.se, roc=roc, sp=specificities, .progress=progress, .parallel=parallel)
  }

  if (any(is.na(perfs))) {
    warning("NA value(s) produced during bootstrap were ignored.")
    perfs <- perfs[!apply(perfs, 1, function(x) any(is.na(x))),]
  }
  ci <- t(apply(perfs, 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
  rownames(ci) <- paste(specificities, ifelse(roc$percent, "%", ""), sep="")
  
  class(ci) <- c("ci.se", "ci", class(ci))
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "specificities") <- specificities
  attr(ci, "roc") <- roc
  return(ci)
}
