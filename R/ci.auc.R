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

ci.auc <- function(...) {
  UseMethod("ci.auc")
}

ci.auc.formula <- function(formula, data, ...) {
	data.missing <- missing(data)
	roc.data <- roc.utils.extract.formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = match.call())
	if (length(roc.data$predictor.name) > 1) {
		stop("Only one predictor supported in 'ci.auc'.")
	}
	response <- roc.data$response
	predictor <- roc.data$predictors[, 1]
	ci.auc.roc(roc.default(response, predictor, ci=FALSE, ...), ...)
}

ci.auc.default <- function(response, predictor, ...) {
	roc <- roc.default(response, predictor, ci = FALSE, ...)
	if (methods::is(roc, "smooth.roc")) {
		return(ci.auc(smooth.roc = roc, ...))
	}
	else {
		return(ci.auc(roc = roc, ...))
	}
}

ci.auc.auc <- function(auc, ...) {
	roc <- attr(auc, "roc")
	roc$auc <- auc
	ci.auc(roc, reuse.auc = TRUE, ...)
}

ci.auc.smooth.roc <- function(smooth.roc,
                   conf.level = 0.95,
                   boot.n = 2000,
                   boot.stratified = TRUE,
                   reuse.auc=TRUE,
                   progress = getOption("pROCProgress")$name,
                   parallel = FALSE,
                   ...
                   ) {
  if (conf.level > 1 | conf.level < 0)
    stop("conf.level must be within the interval [0,1].")
  
  if (roc.utils.is.perfect.curve(smooth.roc)) {
  	warning("ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.")
  }

  # We need an auc
  if (is.null(smooth.roc$auc) | !reuse.auc)
    smooth.roc$auc <- auc(smooth.roc, ...)

  # Check if called with density.cases or density.controls
  if (is.null(smooth.roc$smoothing.args) || is.numeric(smooth.roc$smoothing.args$density.cases) || is.numeric(smooth.roc$smoothing.args$density.controls))
    stop("Cannot compute CI of ROC curves smoothed with numeric density.controls and density.cases.")

  # Get the non smoothed roc.
  roc <- attr(smooth.roc, "roc")
  roc$ci <- NULL # remove potential ci in roc to avoid infinite loop with smooth.roc()

  # do all the computations in fraction, re-transform in percent later if necessary
  percent <- smooth.roc$percent
  smooth.roc$percent <- FALSE
  roc$percent <- FALSE
  oldauc <- smooth.roc$auc
  if (percent) {
    attr(smooth.roc$auc, "percent") <- FALSE
    if (! identical(attr(smooth.roc$auc, "partial.auc"), FALSE)) {
      attr(smooth.roc$auc, "partial.auc") <- attr(smooth.roc$auc, "partial.auc") / 100
    }
  }

  # prepare the calls
  smooth.roc.call <- as.call(c(utils::getS3method("smooth", "roc"), smooth.roc$smoothing.args))
  auc.args <- attributes(smooth.roc$auc)[grep("partial.auc", names(attributes(smooth.roc$auc)))]
  auc.args$allow.invalid.partial.auc.correct <- TRUE
  auc.call <- as.call(c(utils::getS3method("auc", "smooth.roc"), auc.args))

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="AUC confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    aucs <- unlist(llply(1:boot.n, stratified.ci.smooth.auc, roc=roc, smooth.roc.call=smooth.roc.call, auc.call=auc.call, .progress=progress, .parallel=parallel))
  }
  else {
    aucs <- unlist(llply(1:boot.n, nonstratified.ci.smooth.auc, roc=roc, smooth.roc.call=smooth.roc.call, auc.call=auc.call, .progress=progress, .parallel=parallel))
  }

  if (sum(is.na(aucs)) > 0) {
    warning("NA value(s) produced during bootstrap were ignored.")
    aucs <- aucs[!is.na(aucs)]
  }
  # TODO: Maybe apply a correction (it's in the Tibshirani?) What do Carpenter-Bithell say about that?
  # Prepare the return value
  ci <- quantile(aucs, c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2))
  if (percent) {
    ci <- ci * 100
    aucs <- aucs * 100
  }
  attr(ci, "conf.level") <- conf.level
  attr(ci, "method") <- "bootstrap"
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "auc") <- oldauc
  class(ci) <- c("ci.auc", "ci", class(ci))
  return(ci)  
}

ci.auc.roc <- function(roc,
                   conf.level = 0.95,
                   method=c("delong", "bootstrap"),
                   boot.n = 2000,
                   boot.stratified = TRUE,
                   reuse.auc=TRUE,
                   progress = getOption("pROCProgress")$name,
                   parallel = FALSE,
                   ...
                   ) {
  if (conf.level > 1 | conf.level < 0)
    stop("conf.level must be within the interval [0,1].")
  
  if (roc.utils.is.perfect.curve(roc)) {
  	warning("ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.")
  }

  # We need an auc
  if (is.null(roc$auc) | !reuse.auc)
    roc$auc <- auc(roc, ...)

  # do all the computations in fraction, re-transform in percent later if necessary
  percent <- roc$percent
  oldauc <- roc$auc
  if (percent) {
  	roc <- roc.utils.unpercent(roc)
  }

  # Check the method
  if (missing(method) | is.null(method)) {
    # determine method if missing
    if (has.partial.auc(roc)) {
      # partial auc: go for bootstrap
      method <- "bootstrap"
    }
    else if ("smooth.roc" %in% class(roc)) {
      # smoothing: bootstrap
      method <- "bootstrap"
    }
    else {
      method <- "delong"
    }
  }
  else {
    method <- match.arg(method, c("delong", "bootstrap"))
    # delong NA to pAUC: warn + change
    if (has.partial.auc(roc) && method == "delong") {
      stop("DeLong method is not supported for partial AUC. Use method=\"bootstrap\" instead.")
    }
    else if ("smooth.roc" %in% class(roc)) {
      stop("DeLong method is not supported for smoothed ROCs. Use method=\"bootstrap\" instead.")
    }
  }

  if (method == "delong")
    ci <- ci.auc.delong(roc, conf.level)
  else
    ci <- ci.auc.bootstrap(roc, conf.level, boot.n, boot.stratified, progress, parallel, ...)

  if (percent) {
    ci <- ci * 100
  }
  attr(ci, "conf.level") <- conf.level
  attr(ci, "method") <- method
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "auc") <- oldauc
  class(ci) <- c("ci.auc", "ci", class(ci))
  return(ci)
}

ci.auc.multiclass.roc <- function(multiclass.roc, ...) {
	stop("CI of a multiclass ROC curve not implemented")
}

ci.auc.multiclass.auc <- function(multiclass.auc, ...) {
	stop("CI of a multiclass AUC not implemented")
}
