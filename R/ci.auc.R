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

ci.auc.formula <- function(formula, data, ..., subset, na.action) {
  data.missing <- missing(data)
  roc.data <- roc_utils_extract_formula(formula, data, ...,
    data.missing = data.missing,
    call = match.call()
  )
  if (length(roc.data$predictor.names) > 1) {
    stop("Only one predictor supported in 'ci.auc'.")
  }
  response <- roc.data$response
  predictor <- roc.data$predictors[, 1]
  roc <- roc_utils_roc_ci(response, predictor, ci = FALSE, ...)
  if (!roc_utils_is_curve(roc)) {
    return(roc) # NA or NaN from roc(): propagate it
  }
  if (methods::is(roc, "smooth.roc")) {
    return(ci.auc(smooth.roc = roc, ...))
  } else {
    return(ci.auc(roc = roc, ...))
  }
}

ci.auc.default <- function(response, predictor, ...) {
  if (inherits(response, "mv.multiclass.roc")) {
    stop("CI of a multiclass ROC curve not implemented")
  } else if (inherits(response, "mv.multiclass.auc")) {
    stop("CI of a multiclass AUC not implemented")
  }
  roc <- roc_utils_roc_ci(response, predictor, ci = FALSE, ...)
  if (!roc_utils_is_curve(roc)) {
    return(roc) # NA or NaN from roc(): propagate it
  }
  if (methods::is(roc, "smooth.roc")) {
    return(ci.auc(smooth.roc = roc, ...))
  } else {
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
                              reuse.auc = TRUE,
                              progress = getOption("pROCProgress", interactive()),
                              parallel = FALSE,
                              cl = NULL,
                              sides = c("two.sided", "lower", "upper"),
                              ...,
                              method = NULL,
                              scale = c("auc", "logit")) {
  scale.given <- !missing(scale)
  scale <- match.arg(scale)
  if (!roc_utils_is_single_number(conf.level) || conf.level > 1 || conf.level < 0) {
    stop("conf.level must be a single number within the interval [0,1].")
  }
  sides <- match.arg(sides)

  if (roc_utils_is_perfect_curve(smooth.roc)) {
    warning("ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.")
  }

  progress <- roc_utils_normalise_progress(progress)
  roc_utils_warn_deprecated_parallel(parallel)

  # We need an auc
  if (is.null(smooth.roc$auc) | !reuse.auc) {
    # Named: smooth = TRUE in ... (from ci.auc.default/formula) would
    # partially match the smooth.roc argument of auc.smooth.roc
    smooth.roc$auc <- auc(smooth.roc = smooth.roc, ...)
  }

  roc_utils_stop_if_numeric_density(smooth.roc, "CI of")

  # Check the method
  if (is.null(method)) {
    method <- if (roc_utils_has_ml(smooth.roc)) "ml" else "bootstrap"
  } else {
    method <- match.arg(method, c("delong", "bootstrap", "ml"))
    if (method == "delong") {
      warning("DeLong's method is not supported for smoothed ROCs. Using bootstrap instead.")
      method <- "bootstrap"
    } else if (method == "ml") {
      roc_utils_stop_if_no_ml(smooth.roc)
    }
  }
  roc_utils_warn_scale_bootstrap(method, scale.given)

  # Get the non smoothed roc.
  roc <- attr(smooth.roc, "roc")

  # do all the computations in fraction, re-transform in percent later if necessary
  percent <- smooth.roc$percent
  smooth.roc$percent <- FALSE
  roc$percent <- FALSE
  oldauc <- smooth.roc$auc
  if (percent) {
    attr(smooth.roc$auc, "percent") <- FALSE
    if (!identical(attr(smooth.roc$auc, "partial.auc"), FALSE)) {
      attr(smooth.roc$auc, "partial.auc") <- attr(smooth.roc$auc, "partial.auc") / 100
    }
  }

  if (method == "ml") {
    ci <- roc_utils_auc_ci(as.numeric(oldauc) / ifelse(percent, 100, 1), sqrt(var_roc_ml(smooth.roc)), conf.level, scale, sides)
  } else {
    boot.n <- roc_utils_check_boot_n(boot.n)
    aucs <- roc_utils_drop_na_replicates(
      bootstrap.replicates(boot.n, bootstrap.replicate,
        roc = roc, stratified = boot.stratified,
        statistic = bootstrap.auc, auc.args = bootstrap.auc.args(smooth.roc$auc), na = NA_real_,
        smoothing.args = smooth.roc$smoothing.args,
        FUN.VALUE = NA_real_, progress = progress, cl = cl
      )
    )
    ci <- roc_utils_ci_quantile(aucs, conf.level, sides, roc_utils_auc_range(smooth.roc$auc))
  }
  if (percent) {
    ci <- ci * 100
  }
  attr(ci, "conf.level") <- conf.level
  attr(ci, "sides") <- sides
  attr(ci, "method") <- method
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  # The bootstrap percentiles have no scale
  if (method != "bootstrap") {
    attr(ci, "scale") <- scale
  }
  attr(ci, "auc") <- oldauc
  class(ci) <- c("ci.auc", "ci", class(ci))
  return(ci)
}

ci.auc.roc <- function(roc,
                       conf.level = 0.95,
                       method = c("delong", "bootstrap"),
                       boot.n = 2000,
                       boot.stratified = TRUE,
                       reuse.auc = TRUE,
                       progress = getOption("pROCProgress", interactive()),
                       parallel = FALSE,
                       cl = NULL,
                       sides = c("two.sided", "lower", "upper"),
                       ...,
                       scale = c("auc", "logit")) {
  scale.given <- !missing(scale)
  scale <- match.arg(scale)
  if (!roc_utils_is_single_number(conf.level) || conf.level > 1 || conf.level < 0) {
    stop("conf.level must be a single number within the interval [0,1].")
  }
  sides <- match.arg(sides)

  if (roc_utils_is_perfect_curve(roc)) {
    warning("ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.")
  }
  progress <- roc_utils_normalise_progress(progress)
  roc_utils_warn_deprecated_parallel(parallel)

  # We need an auc
  if (is.null(roc$auc) | !reuse.auc) {
    roc$auc <- auc(roc, ...)
  }

  # do all the computations in fraction, re-transform in percent later if necessary
  percent <- roc$percent
  oldauc <- roc$auc
  if (percent) {
    roc <- roc_utils_unpercent(roc)
  }

  # Check the method
  if (missing(method) | is.null(method)) {
    # determine method if missing
    if (has.partial.auc(roc)) {
      # partial auc: go for bootstrap
      method <- "bootstrap"
    } else if ("smooth.roc" %in% class(roc)) {
      # smoothing: bootstrap
      method <- "bootstrap"
    } else if (!is.null(roc$weights)) {
      # observation weights: bootstrap
      method <- "bootstrap"
    } else {
      method <- "delong"
    }
  } else {
    method <- match.arg(method, c("delong", "bootstrap", "ml"))
    if (method == "ml") {
      roc_utils_stop_if_no_ml(roc)
    }
    # delong NA to pAUC: warn + change
    if (has.partial.auc(roc) && method == "delong") {
      stop(roc_utils_no_delong("partial AUC"))
    } else if ("smooth.roc" %in% class(roc)) {
      stop(roc_utils_no_delong("smoothed ROCs"))
    }
  }
  roc_utils_warn_scale_bootstrap(method, scale.given)

  if (method == "delong") {
    ci <- ci_auc_delong(roc, conf.level, scale, sides)
  } else {
    boot.n <- roc_utils_check_boot_n(boot.n)
    ci <- ci_auc_bootstrap(roc, conf.level, sides, boot.n, boot.stratified, progress = progress, cl = cl, ...)
  }

  if (percent) {
    ci <- ci * 100
  }
  attr(ci, "conf.level") <- conf.level
  attr(ci, "sides") <- sides
  attr(ci, "method") <- method
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  # The bootstrap percentiles have no scale
  if (method != "bootstrap") {
    attr(ci, "scale") <- scale
  }
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
