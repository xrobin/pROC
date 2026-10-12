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
    call = match.call()
  )
  if (length(roc.data$predictor.names) > 1) {
    stop("Only one predictor supported in 'ci.coords'.")
  }
  response <- roc.data$response
  predictor <- roc.data$predictors[, 1]
  roc <- roc(response, predictor, ci = FALSE, ...)
  if (!roc_utils_is_curve(roc)) {
    return(roc) # NA or NaN from roc(): propagate it
  }
  if (methods::is(roc, "smooth.roc")) {
    return(ci.coords(smooth.roc = roc, ...))
  } else {
    return(ci.coords(roc = roc, ...))
  }
}

ci.coords.default <- function(response, predictor, ...) {
  if (inherits(response, c("multiclass.roc", "multiclass.auc", "mv.multiclass.roc", "mv.multiclass.auc"))) {
    stop("'ci.coords' not available for multiclass ROC curves.")
  }
  roc <- roc.default(response, predictor, ci = FALSE, ...)
  if (!roc_utils_is_curve(roc)) {
    return(roc) # NA or NaN from roc(): propagate it
  }
  if (methods::is(roc, "smooth.roc")) {
    return(ci.coords(smooth.roc = roc, ...))
  } else {
    return(ci.coords(roc = roc, ...))
  }
}

ci.coords.smooth.roc <- function(smooth.roc,
                                 x,
                                 input = "specificity", ret = c("specificity", "sensitivity"),
                                 best.method = c("youden", "closest.topleft"), best.weights = c(1, 0.5),
                                 best.policy = c("unique.stop", "stop", "omit", "random"),
                                 conf.level = 0.95,
                                 boot.n = 2000,
                                 boot.stratified = TRUE,
                                 progress = getOption("pROCProgress", interactive()),
                                 cl = NULL,
                                 ...) {
  if (!roc_utils_is_single_number(conf.level) || conf.level > 1 || conf.level < 0) {
    stop("'conf.level' must be a single number within the interval [0,1].")
  }

  if (roc_utils_is_perfect_curve(smooth.roc)) {
    warning("ci.coords() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
  progress <- roc_utils_normalise_progress(progress)

  # A smoothed curve has no thresholds, hence threshold = FALSE, as in
  # coords.smooth.roc().
  input <- roc_utils_match_coords_input_args(input, threshold = FALSE)
  ret <- roc_utils_match_coords_ret_args(ret, threshold = FALSE)
  best.policy <- match.arg(best.policy)
  special <- coords_special_x(x, roc = smooth.roc)
  if (!is.na(special)) {
    x <- special
    if (x == "all" || x == "local maximas") {
      stop("'all' and 'local maximas' are not available for confidence intervals.")
    }
  }

  roc_utils_stop_if_numeric_density(smooth.roc, "CI of")

  # Get the non smoothed roc.
  roc <- attr(smooth.roc, "roc")
  auc.args <- NULL
  if (!is.null(smooth.roc$auc) && !identical(attr(smooth.roc$auc, "partial.auc"), FALSE)) {
    auc.args <- attributes(smooth.roc$auc)[c("partial.auc", "partial.auc.focus")]
  }

  # prepare the calls
  smoothing.args <- smooth.roc$smoothing.args
  resmooth <- function(roc) roc_utils_resmooth(roc, smoothing.args)

  boot.n <- roc_utils_check_boot_n(boot.n)
  perfs <- bootstrap.replicates(boot.n, bootstrap.smooth.coords,
    roc = roc, stratified = boot.stratified, x = x, input = input, ret = ret,
    best.method = best.method, best.weights = best.weights,
    resmooth = resmooth, best.policy = best.policy, auc.args = auc.args,
    simplify = "list", progress = progress, cl = cl
  )
  # Reshape into an array of length(x) x length(ret) x boot.n suited for summary
  perfs_array <- array(unlist(perfs),
    dim = c(length(x), length(ret), boot.n),
    dimnames = list(x, ret, NULL)
  )

  if (any(which.ones <- sapply(perfs, function(x) all(is.na(x))))) {
    if (all(which.ones)) {
      warning("All bootstrap iterations produced NA values only.")
    } else {
      how.many <- sum(which.ones)
      warning(sprintf("%s NA value(s) produced during bootstrap were ignored.", how.many))
    }
  }

  summarized.perfs <- apply(perfs_array, 1:2, quantile, probs = c(0 + (1 - conf.level) / 2, .5, 1 - (1 - conf.level) / 2), na.rm = TRUE)
  ci <- lapply(ret, function(x) t(summarized.perfs[, , x]))
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
                          input = "threshold", ret = c("threshold", "specificity", "sensitivity"),
                          best.method = c("youden", "closest.topleft"), best.weights = c(1, 0.5),
                          best.policy = c("unique.stop", "stop", "omit", "random"),
                          conf.level = 0.95,
                          boot.n = 2000,
                          boot.stratified = TRUE,
                          progress = getOption("pROCProgress", interactive()),
                          cl = NULL,
                          ...) {
  if (!roc_utils_is_single_number(conf.level) || conf.level > 1 || conf.level < 0) {
    stop("'conf.level' must be a single number within the interval [0,1].")
  }

  if (roc_utils_is_perfect_curve(roc)) {
    warning("ci.coords() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
  progress <- roc_utils_normalise_progress(progress)

  input <- roc_utils_match_coords_input_args(input)

  best.policy <- match.arg(best.policy)
  special <- coords_special_x(x, roc = roc)
  if (!is.na(special)) {
    x <- special
    if (x == "all" || x == "local maximas") {
      stop("'all' and 'local maximas' are not available for confidence intervals.")
    }
  }

  if (missing(ret) && (input != "threshold" || roc_utils_is_ordered_roc(roc))) {
    # Don't show NA thresholds by default (including ordered curves, where
    # a threshold interval cannot be computed).
    ret <- roc_utils_match_coords_ret_args(setdiff(ret, "threshold"), threshold = FALSE)
  } else if (identical(ret, "all") && !(identical(x, "best") || input == "threshold")) {
    # No threshold interval for a non-threshold input (see the check below):
    # "all" means all the other coordinates, as for smoothed curves
    ret <- roc_utils_match_coords_ret_args(ret, threshold = FALSE)
  } else {
    ret <- roc_utils_match_coords_ret_args(ret)
    if (roc_utils_is_ordered_roc(roc) && "threshold" %in% ret) {
      warning("A confidence interval for 'threshold' is not available for ordered ROC curves.")
    }
  }

  if ("threshold" %in% ret && !(identical(x, "best") || input == "threshold")) {
    stop("'threshold' output is only supported for best ROC point ('x = \"best\"') or if \"threshold\" was given as input.")
  }

  boot.n <- roc_utils_check_boot_n(boot.n)
  perfs <- bootstrap.replicates(boot.n, bootstrap.coords,
    roc = roc, stratified = boot.stratified, x = x, input = input, ret = ret,
    best.method = best.method, best.weights = best.weights,
    best.policy = best.policy, simplify = "list", progress = progress, cl = cl
  )
  perfs <- lapply(perfs, ci_coords_numeric_columns)
  # Reshape into an array of length(x) x length(ret) x boot.n suited for summary
  perfs_array <- array(unlist(perfs),
    dim = c(length(x), length(ret), boot.n),
    dimnames = list(x, ret, NULL)
  )

  if (any(which.ones <- sapply(perfs, function(x) all(is.na(x))))) {
    if (all(which.ones)) {
      warning("All bootstrap iterations produced NA values only.")
    } else {
      how.many <- sum(which.ones)
      warning(sprintf("%s NA value(s) produced during bootstrap were ignored.", how.many))
    }
  }

  summarized.perfs <- apply(perfs_array, 1:2, quantile, probs = c(0 + (1 - conf.level) / 2, .5, 1 - (1 - conf.level) / 2), na.rm = TRUE)
  ci <- lapply(ret, function(x) t(summarized.perfs[, , x]))
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
  if (best.policy == "unique.stop") {
    # Rows that are identical in all the values summarized into the CI are the
    # same point (ie. thresholds at empty levels of an ordered predictor).
    distinct <- !duplicated(ci_coords_numeric_columns(res))
    if (sum(distinct) == 1) {
      return(res[distinct, , drop = FALSE])
    }
    stop("More than one distinct \"best\" threshold was found, aborting. Change 'best.policy' to alter this behavior.")
  } else if (best.policy == "stop") {
    stop("More than one \"best\" threshold was found, aborting. Change 'best.policy' to alter this behavior.")
  } else if (best.policy == "omit") {
    res[1, ] <- NA
    return(res[1, , drop = FALSE])
  } else {
    return(res[sample(seq_len(nrow(res)), size = 1), , drop = FALSE])
  }
}

# Replace non-numeric columns (ie. ordered thresholds) with NA, as they
# cannot be summarized into a CI.
ci_coords_numeric_columns <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      col
    } else {
      rep(NA_real_, length(col))
    }
  })
  df
}
