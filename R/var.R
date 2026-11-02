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

var <- function(...) {
  UseMethod("var")
}

var.default <- function(...) {
  stats::var(...)
}

var.auc <- function(auc, ...) {
  # Change roc from an auc to a roc object but keep the auc specifications
  roc <- auc
  attr(auc, "roc") <- NULL
  roc <- attr(roc, "roc")
  roc$auc <- auc
  # Pass to var.roc
  var.roc(roc, ...)
}

var.smooth.roc <- function(smooth.roc, ...) {
  var.roc(smooth.roc, ...) # just pass to var.roc that will do the job
}

var.roc <- function(roc,
                    method = c("delong", "bootstrap", "obuchowski"),
                    boot.n = 2000,
                    boot.stratified = TRUE,
                    reuse.auc = TRUE,
                    progress = getOption("pROCProgress", interactive()),
                    parallel = FALSE,
                    cl = NULL,
                    ...) {
  # We need an auc
  if (is.null(roc$auc) | !reuse.auc) {
    roc$auc <- auc(roc, ...)
  }

  if (roc_utils_is_perfect_curve(roc)) {
    warning("var() of a ROC curve with AUC == 1 is always 0 and can be misleading.")
  }
  progress <- roc_utils_normalise_progress(progress)
  roc_utils_warn_deprecated_parallel(parallel)

  # do all the computations in fraction, re-transform in percent later
  percent <- roc$percent
  if (percent) {
    roc <- roc_utils_unpercent(roc)
  }

  # Check the method
  if (missing(method) | is.null(method)) {
    # determine method if missing
    if (has.partial.auc(roc)) {
      # partial auc: go for bootstrap
      method <- "bootstrap"
    } else if (inherits(roc, "smooth.roc")) {
      # smoothing: bootstrap
      method <- "bootstrap"
    } else {
      method <- "delong"
    }
  } else {
    method <- match.arg(method)
    # delong NA to pAUC: warn + change
    if (method == "delong") {
      if (has.partial.auc(roc)) {
        stop("DeLong method is not supported for partial AUC. Use method=\"bootstrap\" instead.")
      } else if ("smooth.roc" %in% class(roc)) {
        stop("DeLong method is not supported for smoothed ROCs. Use method=\"bootstrap\" instead.")
      }
    } else if (method == "obuchowski") {
      if ("smooth.roc" %in% class(roc)) {
        stop("Obuchowski method is not supported for smoothed ROCs. Use method=\"bootstrap\" instead.")
      }
      if (has.partial.auc(roc) && attr(roc$auc, "partial.auc.focus") == "sensitivity") {
        stop("Obuchowski method is not supported for partial AUC on sensitivity region. Use method=\"bootstrap\" instead.")
      }
    }
  }

  if (method == "delong") {
    n <- length(roc$controls)
    m <- length(roc$cases)
    V <- delongPlacements(roc)
    var <- var(V$Y) / n + var(V$X) / m
  } else if (method == "obuchowski") {
    var <- var_roc_obuchowski(roc) / length(roc$cases)
  } else {
    var <- var_roc_bootstrap(roc, boot.n, boot.stratified, progress = progress, cl = cl, ...)
  }

  if (percent) {
    var <- var * (100^2)
  }
  return(var)
}

var_roc_bootstrap <- function(roc, boot.n, boot.stratified, progress = FALSE, cl = NULL, ...) {
  ## Smoothed ROC curve variance
  if (inherits(roc, "smooth.roc")) {
    roc_utils_stop_if_numeric_density(roc, "the variance of")
    non.smoothed.roc <- attr(roc, "roc")
    non.smoothed.roc$percent <- FALSE # as we did earlier for the smoothed.roc
    aucs <- bootstrap.replicates(boot.n, bootstrap.replicate,
      roc = non.smoothed.roc, stratified = boot.stratified,
      statistic = bootstrap.auc, auc.args = bootstrap.auc.args(roc$auc), na = NA_real_,
      smoothing.args = roc$smoothing.args,
      FUN.VALUE = NA_real_, progress = progress, cl = cl
    )
  }
  ## Non smoothed ROC curves variance
  else {
    aucs <- bootstrap.replicates(boot.n, bootstrap.replicate,
      roc = roc, stratified = boot.stratified,
      statistic = bootstrap.auc, auc.args = bootstrap.auc.args(roc$auc), na = NA_real_,
      FUN.VALUE = NA_real_, progress = progress, cl = cl
    )
  }

  var(roc_utils_drop_na_replicates(aucs))
}
