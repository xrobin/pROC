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

ci.sp <- function(...) {
  UseMethod("ci.sp")
}

ci.sp.formula <- function(formula, data, ...) {
	data.missing <- missing(data)
	roc.data <- roc.utils.extract.formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = match.call())
	if (length(roc.data$predictor.name) > 1) {
		stop("Only one predictor supported in 'ci.sp'.")
	}
	response <- roc.data$response
	predictor <- roc.data$predictors[, 1]
	ci.sp(roc(response, predictor, ci=FALSE, ...), ...)
}

ci.sp.default <- function(response, predictor, ...) {
	if (methods::is(response, "multiclass.roc") || methods::is(response, "multiclass.auc")) {
		stop("'ci.sp' not available for multiclass ROC curves.")
	}
	roc <- roc.default(response, predictor, ci = FALSE, ...)
	if (methods::is(roc, "smooth.roc")) {
		return(ci.sp(smooth.roc = roc, ...))
	}
	else {
		return(ci.sp(roc = roc, ...))
	}
}

ci.sp.smooth.roc <- function(smooth.roc,
                      sensitivities = seq(0, 1, .1) * ifelse(smooth.roc$percent, 100, 1),
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
  	warning("ci.sp() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }

  # Check if called with density.cases or density.controls
  if (is.null(smooth.roc$smoothing.args) || is.numeric(smooth.roc$smoothing.args$density.cases) || is.numeric(smooth.roc$smoothing.args$density.controls))
    stop("Cannot compute CI of ROC curves smoothed with numeric density.controls and density.cases.")

  # Get the non smoothed roc.
  roc <- attr(smooth.roc, "roc")
  roc$ci <- NULL # remove potential ci in roc to avoid infinite loop with smooth.roc()

  # prepare the calls
  smooth.roc.call <- as.call(c(utils::getS3method("smooth", "roc"), smooth.roc$smoothing.args))

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="SP confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- ldply(1:boot.n, stratified.ci.smooth.sp, roc=roc, se=sensitivities, smooth.roc.call=smooth.roc.call, .progress=progress, .parallel=parallel)
  }
  else {
    perfs <- ldply(1:boot.n, nonstratified.ci.smooth.sp, roc=roc, se=sensitivities, smooth.roc.call=smooth.roc.call, .progress=progress, .parallel=parallel)
  }

  if (any(is.na(perfs))) {
    warning("NA value(s) produced during bootstrap were ignored.")
    perfs <- perfs[!apply(perfs, 1, function(x) any(is.na(x))),]
  }

  ci <- t(apply(perfs, 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
  rownames(ci) <- paste(sensitivities, ifelse(roc$percent, "%", ""), sep="")

  class(ci) <- c("ci.sp", "ci", class(ci))
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "sensitivities") <- sensitivities
  attr(ci, "roc") <- smooth.roc
  return(ci)
}

ci.sp.roc <- function(roc,
                      sensitivities = seq(0, 1, .1) * ifelse(roc$percent, 100, 1),
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
  	warning("ci.sp() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="SP confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- ldply(1:boot.n, stratified.ci.sp, roc=roc, se=sensitivities, .progress=progress, .parallel=parallel)
  }
  else {
    perfs <- ldply(1:boot.n, nonstratified.ci.sp, roc=roc, se=sensitivities, .progress=progress, .parallel=parallel)
  }

  if (any(is.na(perfs))) {
    warning("NA value(s) produced during bootstrap were ignored.")
    perfs <- perfs[!apply(perfs, 1, function(x) any(is.na(x))),]
  }

  ci <- t(apply(perfs, 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
  rownames(ci) <- paste(sensitivities, ifelse(roc$percent, "%", ""), sep="")
  
  class(ci) <- c("ci.sp", "ci", class(ci))
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "sensitivities") <- sensitivities
  attr(ci, "roc") <- roc
  return(ci)
}
