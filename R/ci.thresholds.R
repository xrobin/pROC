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
                   progress = NULL,
                   parallel = FALSE,
                   ...
                   ) {
  if (conf.level > 1 | conf.level < 0)
    stop("'conf.level' must be within the interval [0,1].")
  
  if (roc_utils_is_perfect_curve(roc)) {
  	warning("ci.thresholds() of a ROC curve with AUC == 1 is always a null interval and can be misleading.")
  }
  if ( ! is.null(progress)) {
    warning("Progress bars are deprecated in pROC 1.19. Ignoring 'progress' argument")
  }

  # Check and prepare thresholds
  if (is.character(thresholds)) {
    if (length(thresholds) != 1)
      stop("'thresholds' of class character must be of length 1.")
    thresholds <- match.arg(thresholds, c("all", "best", "local maximas"))
    thresholds.num <- coords(roc, x=thresholds, input="threshold", ret="threshold", ...)[, 1]
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

  perfs_shape <- matrix(NA_real_, nrow=2L, ncol=length(thresholds.num))
  bootstrap_fun <- if (boot.stratified) stratified.ci.thresholds else nonstratified.ci.thresholds
  perfs <- vapply(seq_len(boot.n), bootstrap_fun, FUN.VALUE=perfs_shape, roc=roc, thresholds=thresholds.num)

  probs <- c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)
  # output is length(probs) x 2 x length(thresholds.num)
  perf_quantiles <- apply(perfs, 1:2, quantile, probs=probs)

  sp <- t(perf_quantiles[,1L,])
  se <- t(perf_quantiles[,2L,])

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
