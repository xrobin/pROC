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

ci.thresholds <- function(...) {
  UseMethod("ci.thresholds")
}

ci.thresholds.formula <- function(formula, data, ...) {
	data.missing <- missing(data)
	roc.data <- roc.utils.extract.formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = match.call())
	if (length(roc.data$predictor.name) > 1) {
		stop("Only one predictor supported in 'ci.thresholds'.")
	}
	response <- roc.data$response
	predictor <- roc.data$predictors[, 1]
	ci.thresholds(roc(response, predictor, ci=FALSE, ...), ...)
}

ci.thresholds.default <- function(response, predictor, ...) {
	if (methods::is(response, "multiclass.roc") || methods::is(response, "multiclass.auc")) {
		stop("'ci.thresholds' not available for multiclass ROC curves.")
	}
	ci.thresholds(roc.default(response, predictor, ci=FALSE, ...), ...)
}

ci.thresholds.smooth.roc <- function(smooth.roc, ...)
  stop("'ci.thresholds' is not available for smoothed ROC curves.")

ci.thresholds.roc <- function(roc,
                   conf.level = 0.95,
                   boot.n = 2000,
                   boot.stratified = TRUE,
                   thresholds = "local maximas",
                   progress = getOption("pROCProgress")$name,
                   parallel = FALSE,
                   ...
                   ) {
  if (conf.level > 1 | conf.level < 0)
    stop("'conf.level' must be within the interval [0,1].")
  
  if (roc.utils.is.perfect.curve(roc)) {
  	warning("ci.thresholds() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }

  # Check and prepare thresholds
  if (is.character(thresholds)) {
    if (length(thresholds) != 1)
      stop("'thresholds' of class character must be of length 1.")
    thresholds <- match.arg(thresholds, c("all", "best", "local maximas"))
    thresholds.num <- coords(roc, x=thresholds, input="threshold", ret="threshold", as.matrix = TRUE, transpose = FALSE, ...)[, 1]
    attr(thresholds.num, "coords") <- thresholds
  }
  else if (is.logical(thresholds)) {
    thresholds.num <- roc$thresholds[thresholds]
    attr(thresholds.num, "logical") <- thresholds
  }
  else if (! is.numeric(thresholds)) {
    stop("'thresholds' is not character, logical or numeric.")
  }
  else {
    thresholds.num <- thresholds
  }

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="Thresholds confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- laply(1:boot.n, stratified.ci.thresholds, roc=roc, thresholds=thresholds.num, .progress=progress, .parallel=parallel)
  }
  else {
    perfs <- laply(1:boot.n, nonstratified.ci.thresholds, roc=roc, thresholds=thresholds.num, .progress=progress, .parallel=parallel)
  }

  if (length(thresholds.num) > 1) {
    if (any(is.na(perfs))) {
      warning("NA value(s) produced during bootstrap were ignored.")
      perfs <- perfs[!apply(perfs, 1, function(x) any(is.na(x))),]
    }
    # laply returns a 3d matrix, with dim 1 = bootstrap replicates, dim 2 = SE/SP and dim 3 = thresholds
    # [,1,] = SP and [,2,] = SE
    sp <- t(apply(perfs[,1,], 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
    se <- t(apply(perfs[,2,], 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
  }
  else {
    if (any(is.na(perfs))) {
      warning("NaN value(s) in bootstrap ignored in confidence interval.")
      perfs <- perfs[!apply(perfs, 1, function(x) any(is.na(x))),]
    }
    sp <- as.matrix(t(quantile(perfs[,1], probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2))))
    se <- as.matrix(t(quantile(perfs[,2], probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2))))
  }

  rownames(se) <- rownames(sp) <- thresholds.num

  if (roc$percent) {
    se <- se * 100
    sp <- sp * 100
  }

  ci <- list(specificity = sp, sensitivity = se)
  class(ci) <- c("ci.thresholds", "ci", "list")
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "thresholds") <- thresholds.num
  attr(ci, "roc") <- roc
  return(ci)
}
