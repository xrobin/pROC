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

auc <- function(...) {
  UseMethod("auc")
}

auc.formula <- function(formula, data, ...) {
	data.missing <- missing(data)
	roc.data <- roc.utils.extract.formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = match.call())
	if (length(roc.data$predictor.name) > 1) {
		stop("Only one predictor supported in 'auc'.")
	}
	response <- roc.data$response
	predictor <- roc.data$predictors[, 1]
	
	roc(response, predictor, auc = TRUE, ...)$auc
}

auc.default <- function(response, predictor, ...) {
  roc.default(response, predictor, auc=TRUE, ...)$auc
}

auc.smooth.roc <- function(smooth.roc, ...) {
  auc.roc(smooth.roc, ...) # force usage of auc.roc: compatible
}

auc.multiclass.roc <- function(multiclass.roc, ...) {
  sum <- sum(sapply(multiclass.roc$rocs, auc, ...))
  count <- length(multiclass.roc$levels)
  # Hand & Till formula:
  auc <- (2 * sum) / (count * (count - 1))

  # Prepare auc object
  auc <- as.vector(auc) # remove potential pre-existing attributes
  attr(auc, "percent") <- multiclass.roc$percent
  attr(auc, "roc") <- multiclass.roc
  # Get partial auc details from first computed auc
  # TODO: find a better way to recover partial.auc!
  aucs <- lapply(multiclass.roc$rocs, auc, ...) # keep individual AUCs in a list for later
  attr(auc, "partial.auc") <- attr(aucs[[1]], "partial.auc")
  if (!identical(attr(aucs[[1]], "partial.auc"), FALSE)) {
    attr(auc, "partial.auc.focus") <- attr(aucs[[1]], "partial.auc.focus")
    attr(auc, "partial.auc.correct") <- attr(aucs[[1]], "partial.auc.correct") 
  }
  class(auc) <- c("multiclass.auc", "numeric")
  return(auc)
}

auc.mv.multiclass.roc <- function(mv.multiclass.roc, ...) {
	aucs <- lapply(mv.multiclass.roc$rocs, function(x) list(auc(x[[1]], ...), auc(x[[2]], ...)))
	A.ij.total <- sum(sapply(aucs, function(x) mean(unlist(x))))
	c <- length(mv.multiclass.roc$levels)
	auc <- 2 / (c * (c-1)) * A.ij.total
	
	# Prepare auc object
	auc <- as.vector(auc) # remove potential pre-existing attributes
	attr(auc, "percent") <- mv.multiclass.roc$percent
	attr(auc, "roc") <- mv.multiclass.roc
	
	# Get partial auc details from first computed auc
	attr(auc, "partial.auc") <- attr(aucs[[1]][[1]], "partial.auc")
	if (!identical(attr(aucs[[1]], "partial.auc"), FALSE)) {
		attr(auc, "partial.auc.focus") <- attr(aucs[[1]][[1]], "partial.auc.focus")
		attr(auc, "partial.auc.correct") <- attr(aucs[[1]][[1]], "partial.auc.correct") 
	}
	class(auc) <- c("mv.multiclass.auc", "numeric")
	return(auc)
}

auc.roc <- function(roc,
                    # Partial auc definition
                    partial.auc=FALSE, # false (consider total area) or numeric length 2: boundaries of the AUC to consider, between 0 and 1, or 0 and 100 if percent is TRUE
                    partial.auc.focus=c("specificity", "sensitivity"), # if partial.auc is not FALSE: do the boundaries
                    partial.auc.correct=FALSE,
					allow.invalid.partial.auc.correct = FALSE,
                    ... # unused required to allow roc passing arguments to plot or ci.
                    ) {
  if (!identical(partial.auc, FALSE)) {
    partial.auc.focus <- match.arg(partial.auc.focus)
  }

  percent <- roc$percent
  
  # Validate partial.auc
  if (! identical(partial.auc, FALSE) & !(is.numeric(partial.auc) && length(partial.auc)==2))
    stop("partial.auc must be either FALSE or a numeric vector of length 2")
  
  # Ensure partial.auc is sorted with partial.auc[1] >= partial.auc[2]
  partial.auc <- sort(partial.auc, decreasing=TRUE)
  # Get and sort the sensitivities and specificities
  roc <- sort(roc)
  se <- roc$sensitivities
  sp <- roc$specificities

  # Full area if partial.auc is FALSE
  if (identical(partial.auc, FALSE)) {
    if (methods::is(roc, "smooth.roc") && ! is.null(roc$smoothing.args) && roc$smoothing.args$method == "binormal") {
      coefs <- coefficients(roc$model)
      auc <- unname(pnorm(coefs[1] / sqrt(1+coefs[2]^2)) * ifelse(percent, 100^2, 1))
    }
    else {
      diffs.x <- sp[-1] - sp[-length(sp)]
      means.vert <- (se[-1] + se[-length(se)])/2
      auc <- sum(means.vert * diffs.x)
    }
  }
  # Partial area
  else {
    if (partial.auc.focus == "sensitivity") {
      # if we focus on SE, just swap and invert x and y and the computations for SP will work
      x <- rev(se)
      y <- rev(sp)
    }
    else {
      x <- sp
      y <- se
    }
    
    # find the SEs and SPs in the interval
    x.inc <- x[x <= partial.auc[1] & x >= partial.auc[2]]
    y.inc <- y[x <= partial.auc[1] & x >= partial.auc[2]]
    # compute the AUC strictly in the interval
    diffs.x <- x.inc[-1] - x.inc[-length(x.inc)]
    means.vert <- (y.inc[-1] + y.inc[-length(y.inc)])/2
    auc <- sum(means.vert * diffs.x)
    # add the borders:
    if (length(x.inc) == 0) { # special case: the whole AUC is between 2 se/sp points. Need to interpolate from both
      diff.horiz <- partial.auc[1] - partial.auc[2]
      # determine indices
      idx.hi <- match(FALSE, x < partial.auc[1])
      idx.lo <- idx.hi - 1
      # proportions
      proportion.hi <- (x[idx.hi] - partial.auc[1]) / (x[idx.hi] - x[idx.lo])
      proportion.lo <- (partial.auc[2] - x[idx.lo]) / (x[idx.hi] - x[idx.lo])
      # interpolated y's
      y.hi <- y[idx.hi] + proportion.hi * (y[idx.lo] - y[idx.hi])
      y.lo <- y[idx.lo] - proportion.lo * (y[idx.lo] - y[idx.hi])
      # compute AUC
      mean.vert <- (y.hi + y.lo)/2
      auc <- mean.vert*diff.horiz
    }
    else { # if the upper limit is not exactly present in SPs, interpolate
      if (!(partial.auc[1] %in% x.inc)) {
        # find the limit indices
        idx.out <- match(FALSE, x < partial.auc[1])
        idx.in <- idx.out - 1
        # interpolate y
        proportion <- (partial.auc[1] - x[idx.out]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.out] + proportion * (y[idx.in] - y[idx.out])
        # add to AUC
        auc <- auc + (partial.auc[1] - x[idx.in]) * (y[idx.in] + y.interpolated)/2
      }
      if (!(partial.auc[2] %in% x.inc)) { # if the lower limit is not exactly present in SPs, interpolate
        # find the limit indices in and out
        #idx.out <- length(x) - match(TRUE, rev(x) < partial.auc[2]) + 1
        idx.out <- match(TRUE, x > partial.auc[2]) - 1
        idx.in <- idx.out + 1
        # interpolate y
        proportion <- (x[idx.in] - partial.auc[2]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.in] + proportion * (y[idx.out] - y[idx.in])
        # add to AUC
        auc <- auc + (x[idx.in] - partial.auc[2]) * (y[idx.in] + y.interpolated)/2
      }
    }
  }

  # In percent, we have 100*100 = 10,000 as maximum area, so we need to divide by a factor 100
  if (percent)
    auc <- auc/100

  # Correction according to McClish DC, 1989
  if (all(!identical(partial.auc, FALSE), partial.auc.correct)) { # only for pAUC
    min <- roc.utils.min.partial.auc(partial.auc, percent)
    max <- roc.utils.max.partial.auc(partial.auc, percent)
    # The correction is defined only when auc >= min
    if (!allow.invalid.partial.auc.correct && auc < min) {
    	warning("Partial AUC correction not defined for ROC curves below the diagonal.")
    	auc <- NA
    }
    else if (percent) {
      auc <- (100+((auc-min)*100/(max-min)))/2 # McClish formula adapted for %
    }
    else {
      auc <- (1+((auc-min)/(max-min)))/2 # original formula by McClish
    }
  }
  # Prepare the AUC to return with attributes
  auc <- as.vector(auc) # remove potential pre-existing attributes
  attr(auc, "partial.auc") <- partial.auc
  attr(auc, "percent") <- percent
  attr(auc, "roc") <- roc
  if (!identical(partial.auc, FALSE)) {
    attr(auc, "partial.auc.focus") <- partial.auc.focus
    attr(auc, "partial.auc.correct") <- partial.auc.correct
  }
  class(auc) <- c("auc", class(auc))
  return(auc)
}
