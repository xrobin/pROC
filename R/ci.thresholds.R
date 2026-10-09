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
  roc.data <- roc_utils_extract_formula(formula, data, ...,
    data.missing = data.missing,
    call = match.call()
  )
  if (length(roc.data$predictor.name) > 1) {
    stop("Only one predictor supported in 'ci.thresholds'.")
  }
  response <- roc.data$response
  predictor <- roc.data$predictors[, 1]
  ci.thresholds(roc(response, predictor, ci = FALSE, ...), ...)
}

ci.thresholds.default <- function(response, predictor, ...) {
  if (methods::is(response, "multiclass.roc") || methods::is(response, "multiclass.auc")) {
    stop("'ci.thresholds' not available for multiclass ROC curves.")
  }
  ci.thresholds(roc.default(response, predictor, ci = FALSE, ...), ...)
}

ci.thresholds.smooth.roc <- function(smooth.roc, ...) {
  stop("'ci.thresholds' is not available for smoothed ROC curves.")
}

ci.thresholds.roc <- function(roc,
                              conf.level = 0.95,
                              boot.n = 2000,
                              boot.stratified = TRUE,
                              thresholds = "local maximas",
                              progress = getOption("pROCProgress", interactive()),
                              cl = NULL,
                              parallel = FALSE,
                              ...) {
  if (conf.level > 1 | conf.level < 0) {
    stop("'conf.level' must be within the interval [0,1].")
  }

  if (roc_utils_is_perfect_curve(roc)) {
    warning("ci.thresholds() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
  progress <- roc_utils_normalise_progress(progress)
  roc_utils_warn_deprecated_parallel(parallel)

  # Check and prepare thresholds
  special <- coords_special_x(thresholds, roc = roc, keywords = c("all", "best", "local maximas"))
  if (!is.na(special)) {
    thresholds <- special
    thresholds.num <- coords(roc, x = thresholds, input = "threshold", ret = "threshold", ...)[, 1]
    attr(thresholds.num, "coords") <- thresholds
  } else if (is.logical(thresholds)) {
    thresholds.num <- roc$thresholds[thresholds]
    attr(thresholds.num, "logical") <- thresholds
  } else if (is.numeric(thresholds) || is.ordered(thresholds) || is.character(thresholds)) {
    ordered_roc <- roc_utils_is_ordered_roc(roc)
    if (is.numeric(thresholds) && ordered_roc) {
      stop("Numeric 'thresholds' is not supported for ordered ROC curves. Pass level labels as character or ordered.")
    }
    if ((is.character(thresholds) || is.ordered(thresholds)) && !ordered_roc) {
      if (is.character(thresholds) && length(thresholds) != 1) {
        stop("'thresholds' of class character must be of length 1.")
      }
      stop("'thresholds' must be a numeric or one of \"all\", \"local maximas\", \"best\".")
    }
    thresholds.num <- thresholds
    if (ordered_roc) {
      # Resolve to a T-scale ordered vector so the bootstrap loop can rank
      # by as.integer() rather than by label.
      thr_idx <- roc_utils_thr_idx(roc, thresholds.num)
      thresholds.num <- roc$thresholds[thr_idx]
    }
  } else {
    stop("'thresholds' is not character, logical, ordered or numeric.")
  }

  # One 2 x length(thresholds) matrix per replicate, stacked into a
  # 2 x length(thresholds) x boot.n array.
  perfs <- bootstrap.replicates(boot.n, bootstrap.thresholds,
    roc = roc, stratified = boot.stratified, thresholds = thresholds.num,
    simplify = "columns", progress = progress, cl = cl
  )
  # Drop the replicates with NA (a class can vanish from a non-stratified
  # resample): one column per replicate, then back to the array
  perfs.dim <- dim(perfs)
  perfs <- roc_utils_drop_na_replicates(matrix(perfs, ncol = perfs.dim[3]), margin = 2L)
  perfs <- array(perfs, dim = c(perfs.dim[1:2], ncol(perfs)))

  probs <- c(0 + (1 - conf.level) / 2, .5, 1 - (1 - conf.level) / 2)
  # output is length(probs) x 2 x length(thresholds.num)
  perf_quantiles <- apply(perfs, 1:2, quantile, probs = probs)

  sp <- t(perf_quantiles[, 1L, ])
  se <- t(perf_quantiles[, 2L, ])

  rownames(se) <- rownames(sp) <- as.character(thresholds.num)

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
