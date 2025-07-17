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
	data.missing <- missing(data)
	roc.data <- roc_utils_extract_formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = match.call())
	if (length(roc.data$predictor.name) > 1) {
		stop("Only one predictor supported in 'ci.coords'.")
	}
	response <- roc.data$response
	predictor <- roc.data$predictors[, 1]
  	ci.coords(roc(response, predictor, ci=FALSE, ...), ...)
}

ci.coords.default <- function(response, predictor, ...) {
	if (methods::is(response, "multiclass.roc") || methods::is(response, "multiclass.auc")) {
		stop("'ci.coords' not available for multiclass ROC curves.")
	}
	roc <- roc.default(response, predictor, ci = FALSE, ...)
	if (methods::is(roc, "smooth.roc")) {
		return(ci.coords(smooth.roc = roc, ...))
	}
	else {
		return(ci.coords(roc = roc, ...))
	}
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
  
  if (roc_utils_is_perfect_curve(smooth.roc)) {
  	warning("ci.coords() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
	
  input <- roc_utils_match_coords_input_args(input)
  ret <- roc_utils_match_coords_ret_args(ret)
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

  if(inherits(progress, "list"))
    progress <- roc_utils_get_progress_bar(progress, title="Coords confidence interval", label="Bootstrap in progress...", ...)

  smooth_coords_fun <- if (boot.stratified) stratified.ci.smooth.coords else nonstratified.ci.smooth.coords
  # Replicate with simplify=FALSE returns a list of length boot.n
  perfs <- replicate(boot.n, smooth_coords_fun(roc, x, input, ret, best.method, best.weights, smooth.roc.call, best.policy), simplify=FALSE)
  # Reshape into an array of length(x) x length(ret) x boot.n suited for summary
  perfs_array <- array(unlist(perfs),
  					 dim = c(length(x), length(ret), boot.n),
  					 dimnames = list(x, ret, NULL))
  
  if (any(which.ones <- sapply(perfs, function(x) all(is.na(x))))) {
  	if (all(which.ones)) {
  		warning("All bootstrap iterations produced NA values only.")
  	}
  	else {
  		how.many <- sum(which.ones)
  		warning(sprintf("%s NA value(s) produced during bootstrap were ignored.", how.many))
  	}
  }

  summarized.perfs <- apply(perfs_array, 1:2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2), na.rm=TRUE)
  ci <- lapply(ret, function(x) t(summarized.perfs[,,x]))
  names(ci) <- ret

  class(ci) <- c("ci.coords", "ci", class(ci))
  attr(ci, "input") <- input
  attr(ci, "x") <- x
  attr(ci, "ret") <- ret
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "roc") <- smooth.roc
  return(ci)
}

ci.coords.roc <- function(roc,
								  x,
								  input="threshold", ret=c("threshold", "specificity", "sensitivity"),
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
  
  if (roc_utils_is_perfect_curve(roc)) {
  	warning("ci.coords() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
	
  input <- roc_utils_match_coords_input_args(input)
  
  if (missing(ret) && input != "threshold") {
  	# Don't show NA thresholds by default
  	ret <- roc_utils_match_coords_ret_args(ret, threshold = FALSE)
  }
  else {
  	ret <- roc_utils_match_coords_ret_args(ret)
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

  if(inherits(progress, "list"))
    progress <- roc_utils_get_progress_bar(progress, title="Coords confidence interval", label="Bootstrap in progress...", ...)

  coords_fun <- if (boot.stratified) stratified.ci.coords else nonstratified.ci.coords
  # Replicate with simplify=FALSE returns a list of length boot.n
  perfs <- replicate(boot.n, coords_fun(roc, x, input, ret, best.method, best.weights, best.policy), simplify=FALSE)
  # Reshape into an array of length(x) x length(ret) x boot.n suited for summary
  perfs_array <- array(unlist(perfs),
  					 dim = c(length(x), length(ret), boot.n),
  					 dimnames = list(x, ret, NULL))
  
  if (any(which.ones <- sapply(perfs, function(x) all(is.na(x))))) {
  	if (all(which.ones)) {
  		warning("All bootstrap iterations produced NA values only.")
  	}
  	else {
  		how.many <- sum(which.ones)
  		warning(sprintf("%s NA value(s) produced during bootstrap were ignored.", how.many))
  	}
  }
  
  summarized.perfs <- apply(perfs_array, 1:2, quantile, probs=c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2), na.rm=TRUE)
  ci <- lapply(ret, function(x) t(summarized.perfs[,,x]))
  names(ci) <- ret

  class(ci) <- c("ci.coords", "ci", class(ci))
  attr(ci, "input") <- input
  attr(ci, "x") <- x
  attr(ci, "ret") <- ret
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "roc") <- roc
  return(ci)
}

# Function to be called when "best" threshold returned more than 1 column
# Will follow the action defined by best.policy
# For instance:
#   if (x == "best" && nrow(res) != 1) {
# return(enforce.best.policy(res, best.policy))
# }
enforce.best.policy <- function(res, best.policy) {
	if (best.policy == "stop") {
		stop("More than one \"best\" threshold was found, aborting. Change 'best.policy' to alter this behavior.")
	}
	else if (best.policy == "omit") {
		res[1, ] <- NA
		return(res[1, drop = FALSE])
	}
	else {
		return(res[sample(seq_len(nrow(res)), size = 1), , drop = FALSE])
	}
}

