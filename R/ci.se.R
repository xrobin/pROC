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
  data.missing <- missing(data)
  roc.data <- roc_utils_extract_formula(formula, data, ...,
    data.missing = data.missing,
    call = match.call()
  )
  if (length(roc.data$predictor.names) > 1) {
    stop("Only one predictor supported in 'ci.se'.")
  }
  response <- roc.data$response
  predictor <- roc.data$predictors[, 1]
  roc <- roc(response, predictor, ci = FALSE, ...)
  if (!roc_utils_is_curve(roc)) {
    return(roc) # NA or NaN from roc(): propagate it
  }
  if (methods::is(roc, "smooth.roc")) {
    return(ci.se(smooth.roc = roc, ...))
  } else {
    return(ci.se(roc = roc, ...))
  }
}

ci.se.default <- function(response, predictor, ...) {
  if (inherits(response, c("multiclass.roc", "multiclass.auc", "mv.multiclass.roc", "mv.multiclass.auc"))) {
    stop("'ci.se' not available for multiclass ROC curves.")
  }
  roc <- roc.default(response, predictor, ci = FALSE, ...)
  if (!roc_utils_is_curve(roc)) {
    return(roc) # NA or NaN from roc(): propagate it
  }
  if (methods::is(roc, "smooth.roc")) {
    return(ci.se(smooth.roc = roc, ...))
  } else {
    return(ci.se(roc = roc, ...))
  }
}

ci.se.smooth.roc <- function(smooth.roc,
                             specificities = seq(0, 1, .1) * ifelse(smooth.roc$percent, 100, 1),
                             conf.level = 0.95,
                             boot.n = 2000,
                             boot.stratified = TRUE,
                             progress = getOption("pROCProgress", interactive()),
                             cl = NULL,
                             parallel = FALSE,
                             ...) {
  if (!roc_utils_is_single_number(conf.level) || conf.level > 1 || conf.level < 0) {
    stop("'conf.level' must be a single number within the interval [0,1].")
  }

  if (roc_utils_is_perfect_curve(smooth.roc)) {
    warning("ci.se() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
  progress <- roc_utils_normalise_progress(progress)
  roc_utils_warn_deprecated_parallel(parallel)

  roc_utils_stop_if_numeric_density(smooth.roc, "CI of")

  # Get the non smoothed roc.
  roc <- attr(smooth.roc, "roc")

  # prepare the calls
  smoothing.args <- smooth.roc$smoothing.args
  resmooth <- function(roc) roc_utils_resmooth(roc, smoothing.args)

  boot.n <- roc_utils_check_boot_n(boot.n)
  perfs <- bootstrap.replicates(boot.n, bootstrap.smooth.se,
    roc = roc, stratified = boot.stratified, sp = specificities,
    resmooth = resmooth, simplify = "rows", progress = progress, cl = cl
  )

  perfs <- roc_utils_drop_na_replicates(perfs, margin = 1L, ignore = is.na(specificities))

  ci <- t(apply(perfs, 2, quantile, probs = c(0 + (1 - conf.level) / 2, .5, 1 - (1 - conf.level) / 2), na.rm = TRUE))
  rownames(ci) <- paste(specificities, ifelse(roc$percent, "%", ""), sep = "")

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
                      progress = getOption("pROCProgress", interactive()),
                      cl = NULL,
                      parallel = FALSE,
                      ...) {
  if (!roc_utils_is_single_number(conf.level) || conf.level > 1 || conf.level < 0) {
    stop("'conf.level' must be a single number within the interval [0,1].")
  }

  if (roc_utils_is_perfect_curve(roc)) {
    warning("ci.se() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
  progress <- roc_utils_normalise_progress(progress)
  roc_utils_warn_deprecated_parallel(parallel)

  boot.n <- roc_utils_check_boot_n(boot.n)
  perfs <- bootstrap.replicates(boot.n, bootstrap.se,
    roc = roc, stratified = boot.stratified, sp = specificities, simplify = "rows", progress = progress, cl = cl
  )

  perfs <- roc_utils_drop_na_replicates(perfs, margin = 1L, ignore = is.na(specificities))
  ci <- t(apply(perfs, 2, quantile, probs = c(0 + (1 - conf.level) / 2, .5, 1 - (1 - conf.level) / 2), na.rm = TRUE))
  rownames(ci) <- paste(specificities, ifelse(roc$percent, "%", ""), sep = "")

  class(ci) <- c("ci.se", "ci", class(ci))
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "specificities") <- specificities
  attr(ci, "roc") <- roc
  return(ci)
}
