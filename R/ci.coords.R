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

ci.coords <- function(...) {
  UseMethod("ci.coords")
}

ci.coords.formula <- function(formula, data, ...) {
	# Get the data. Use standard code from survival::coxph as suggested by Terry Therneau
	Call <- match.call()
	indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(Call), nomatch=0)
	if (indx[1] == 0) {
		stop("A formula argument is required")
	}
	# Keep the standard arguments and run them in model.frame
	temp <- Call[c(1,indx)]  
	temp[[1]] <- as.name('model.frame')
	m <- eval(temp, parent.frame())
	
	if (!is.null(model.weights(m))) stop("weights are not supported")
	
	response <- model.response(m)
	predictor <- m[[attr(terms(formula), "term.labels")]]
  	ci.coords(roc(response, predictor, ci=FALSE, ...), ...)
}

ci.coords.default <- function(response, predictor, ...) {
	if (methods::is(response, "multiclass.roc") || methods::is(response, "multiclass.auc")) {
		stop("'ci.coords' not available for multiclass ROC curves.")
	}
  ci.coords(roc.default(response, predictor, ci=FALSE, ...), ...)
}

ci.coords.smooth.roc <- function(smooth.roc,
											x, 
											input=c("specificity", "sensitivity"), ret=c("specificity", "sensitivity"), 
											best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5),
											best.policy = c("stop", "omit", "random"),
											conf.level = 0.95,
											boot.n = 2000,
											boot.stratified = TRUE,
											progress = getOption("pROCProgress")$name,
											...
                      ) {
  if (conf.level > 1 | conf.level < 0)
    stop("'conf.level' must be within the interval [0,1].")
  
  if (roc.utils.is.perfect.curve(smooth.roc)) {
  	warning("ci.coords() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
 
  input <- match.arg(input)
  ret <- roc.utils.match.coords.ret.args(ret)
  best.policy <- match.arg(best.policy)
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
  smooth.roc.call <- as.call(c(utils::getS3method("smooth", "roc"), smooth.roc$smoothing.args))

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="Coords confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- raply(boot.n, stratified.ci.smooth.coords(roc, x, input, ret, best.method, best.weights, smooth.roc.call, best.policy), .progress=progress)
  }
  else {
    perfs <- raply(boot.n, nonstratified.ci.smooth.coords(roc, x, input, ret, best.method, best.weights,smooth.roc.call, best.policy), .progress=progress)
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

  if (length(x) > 1 && length(ret) > 1) {
  	quant.perfs <- apply(perfs, c(2, 3), quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2))
  	ci <- matrix(quant.perfs, ncol=3, byrow=TRUE)
  	colnames(ci) <- dimnames(quant.perfs)[[1]]
  }
  else {
  	ci <- t(apply(perfs, 2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
  }
  
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
								  best.policy = c("stop", "omit", "random"),
								  conf.level = 0.95,
								  boot.n = 2000,
								  boot.stratified = TRUE,
								  progress = getOption("pROCProgress")$name,
                      ...
                      ) {
  if (conf.level > 1 | conf.level < 0)
    stop("'conf.level' must be within the interval [0,1].")
  
  if (roc.utils.is.perfect.curve(roc)) {
  	warning("ci.coords() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
 
  input <- match.arg(input)
  
  if (missing(ret) && input != "threshold") {
  	# Don't show NA thresholds by default
  	ret <- roc.utils.match.coords.ret.args(ret, threshold = FALSE)
  }
  else {
  	ret <- roc.utils.match.coords.ret.args(ret)
  }
  
  best.policy <- match.arg(best.policy)
  if (is.character(x)) {
    x <- match.arg(x, c("all", "local maximas", "best"))
    if (x == "all" || x == "local maximas") {
      stop("'all' and 'local maximas' are not available for confidence intervals.")
    }
  }
  
  if ("threshold" %in% ret && ! (identical(x, "best") || input == "threshold")) {
  	stop("'threshold' output is only supported for best ROC point ('x = \"best\"') or if \"threshold\" was given as input.")
  }

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="Coords confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    perfs <- rlply(boot.n, stratified.ci.coords(roc, x, input, ret, best.method, best.weights, best.policy), .progress=progress)
  }
  else {
    perfs <- raply(boot.n, nonstratified.ci.coords(roc, x, input, ret, best.method, best.weights, best.policy), .progress=progress)
  }

  if (any(which.ones <- sapply(perfs, function(x) all(is.na(x))))) {
  	if (all(which.ones)) {
  		warning("All bootstrap iterations produced NA values only.")
  	}
  	else {
  		how.many <- sum(which.ones)
  		warning(sprintf("%s NA value(s) produced during bootstrap were ignored.", how.many))
  	}
  }

  if (length(x) > 1 || length(ret) > 1) {
    ci <- t(apply(sapply(perfs, c), 1, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2), na.rm=TRUE))
  }
  else { # 1 x and 1 ret
    # If x == "best" coords may return multiple best thresholds
    # Be very conservative and take the most extreme ones
    ci.max <- quantile(sapply(perfs, max), probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2), na.rm=TRUE)
    ci.min <- quantile(sapply(perfs, min), probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2), na.rm=TRUE)
    ci <- t(ci.max * c(0, 0.5, 1) + ci.min * c(1, 0.5, 0))
  }
  rownames(ci) <- paste(rep(paste(input, x), each=length(ret)), ret, sep=": ")

  class(ci) <- c("ci.coords", "ci", class(ci))
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "roc") <- roc
  return(ci)
}

# Function to be called when "best" threshold returned more than 1 column
# Will follow the action defined by best.policy
# For instance:
#   if (x == "best" && ncol(res) != 1) {
# return(enfore.best.policy(res, best.policy))
# }
enfore.best.policy <- function(res, best.policy) {
	if (best.policy == "stop") {
		stop("More than one \"best\" threshold was found, aborting. Change 'best.policy' to alter this behavior.")
	}
	else if (best.policy == "omit") {
		res[, 1] <- NA
		return(res[, 1, drop = FALSE])
	}
	else {
		return(res[, sample(seq_len(ncol(res)), size = 1), drop = FALSE])
	}
}

