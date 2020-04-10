# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2010-2019 Xavier Robin, Matthias Doering, 
# Alexandre Hainard, Natacha Turck, Natalia Tiberti, 
# Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller
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

multiclass.roc <- function(...)
  UseMethod("multiclass.roc")

multiclass.roc.formula <- function(formula, data, ...) {
	data.missing <- missing(data)
	call <- match.call()
	roc.data <- roc.utils.extract.formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = call)
	response <- roc.data$response
	predictors <- roc.data$predictors
	if (ncol(predictors) == 1) {
		predictors <- predictors[, 1]
	}
	multiclass.roc <- multiclass.roc.default(response, predictors, ...)
	multiclass.roc$call <- call
	if (! data.missing) {
		multiclass.roc$data <- data
	}
	return(multiclass.roc)
}

multiclass.roc.univariate <- function(response, predictor,
                                   levels=base::levels(as.factor(response)),
                                   percent=FALSE, # Must sensitivities, specificities and AUC be reported in percent? Note that if TRUE, and you want a partial area, you must pass it in percent also (partial.area=c(100, 80))
                                   direction,
                                   # what computation must be done
                                   #auc=TRUE, # call auc.roc on the current object
                                   #ci=FALSE, # call ci.roc on the current object
                                   ...) {
  multiclass.roc <- list(
                         response = response,
                         predictor = predictor,
                         percent = percent)
  class(multiclass.roc) <- "multiclass.roc"
  if (is.factor(response) && any(names(table(response))[table(response) == 0] %in% levels)) {
    missing.levels <- names(table(response))[table(response) == 0]
    missing.levels.requested <- missing.levels[missing.levels %in% levels]
    warning(paste("No observation for response level(s):", paste(missing.levels.requested, collapse=", ")))
    levels <- levels[!(levels %in% missing.levels.requested)]
  }
  multiclass.roc$levels <- levels
  
  rocs <- utils::combn(levels, 2, function(X, response, predictor, percent, ...) {
    roc(response, predictor, levels=X, percent=percent, auc=FALSE, ci=FALSE, ...)
  }, simplify=FALSE, response=response, predictor=predictor, percent=percent, direction=direction, ...)

  multiclass.roc$rocs <- rocs

  # Makes no sense to turn auc off, so remove this option
  #if (auc)
    multiclass.roc$auc <- auc.multiclass.roc(multiclass.roc, ...)
  # CI is not implemented yet.
  #if (ci)
  #  multiclass.roc$ci <- ci.multiclass.roc(multiclass.roc, ...)

  return(multiclass.roc)
}

compute.pair.AUC <- function(pred.matrix, i, j, ref.outcome, levels, percent, direction, ... ) {
    # computes A(i|j), the probability that a randomly 
    # chosen member of class j has a lower estimated probability (or score) 
    # of belonging to class i than a randomly chosen member of class i
    pred.i <- pred.matrix[which(ref.outcome == i), i] # p(G = i) assigned to class i observations
    pred.j <- pred.matrix[which(ref.outcome == j), i] # p(G = i) assigned to class j observations
    classes <- factor(c(rep(i, length(pred.i)), rep(j, length(pred.j))))
    # override levels argument by new levels
    levels <- unique(classes)
    predictor <- c(pred.i, pred.j)
    auc <- roc(classes, predictor, levels = levels, percent = percent, auc = FALSE, ci = FALSE, direction = direction, ...)
    return(auc)
}

multiclass.roc.multivariate <- function(response, predictor, levels, percent, direction, ...) {
    # Reference: "A Simple Generalisation of the Area Under the ROC 
    # Curve for Multiple Class Classification Problems" (Hand and Till, 2001)
    if (!methods::is(predictor, "matrix") && !methods::is(predictor, "data.frame")) {
        stop("Please provide a matrix or data frame via 'predictor'.")
    }
    if (nrow(predictor) != length(response)) {
        stop("Number of rows in 'predictor' does not agree with 'response'");
    }
	if (direction == "auto") {
		stop("'direction=\"auto\"' not available for multivariate multiclass.roc")
	}
	if (is.factor(response) && any(names(table(response))[table(response) == 0] %in% levels)) {
		missing.levels <- names(table(response))[table(response) == 0]
		missing.levels.requested <- missing.levels[missing.levels %in% levels]
		warning(paste("No observation for response level(s):", paste(missing.levels.requested, collapse=", ")))
		levels <- levels[!(levels %in% missing.levels.requested)]
	}
	
    # check whether the columns of the prediction matrix agree with the factors in 'response'
    m <- match(colnames(predictor), levels)
    missing.classes <- levels[setdiff(seq_along(levels), m)]
    levels <- colnames(predictor)[!is.na(m)]
    if (length(levels) == 1) {
        stop("For a single decision value, please provide 'predictor' as a vector.")
    } else if (length(levels) == 0) {
        stop("The column names of 'predictor' could not be matched to the levels of 'response'.")
    }
    if (length(missing.classes) != 0) {
        out.classes <- paste0(missing.classes, collapse = ",")
        if (length(missing.classes) == length(levels)) {
            # no decision values found
            stop(paste0("Could not find any decision values in 'predictor' matching the 'response' levels.",
                 " Could not find the following classes: ", out.classes, ". Check your column names!"))
        } else {
            # some decision values not found
            warning("You did not provide decision values for the following classes: ", out.classes, ".")
        }
    }
    additional.classes <- colnames(predictor)[which(is.na(m))]
    if (length(additional.classes) != 0) {
        out.classes <- paste0(additional.classes, collapse = ",")
        warning("The following classes were not found in 'response': ", out.classes, ".")
    }
    multiclass.roc <- list(
                         response = response,
                         predictor = predictor,
                         percent = percent)
    class(multiclass.roc) <- "mv.multiclass.roc"
    multiclass.roc$levels <- levels
    rocs <- utils::combn(levels, 2, function(x, predictor, response, levels, percent, direction, ...) {
        A1 <- compute.pair.AUC(predictor, x[1], x[2], response, levels, percent, direction, ...)
        A2 <- compute.pair.AUC(predictor, x[2], x[1], response, levels, percent, direction, ...)
        # merging A1 and A2 is infeasible as auc() would not be well-defined
        A <- list(A1, A2) 
        return(A)
    }, simplify = FALSE, predictor = predictor, response = response, 
        levels = levels, percent = percent, direction, ...)
    pairs <- unlist(lapply(utils::combn(levels, 2, simplify = FALSE), 
                           function(x) paste(x, collapse = "/")))
    names(rocs) <- pairs
    multiclass.roc$rocs <- rocs
    
    multiclass.roc$auc <- auc.mv.multiclass.roc(multiclass.roc, ...)
    
    return(multiclass.roc)
}

multiclass.roc.default <- function(response, predictor,
                                   levels = base::levels(as.factor(response)),
                                   percent = FALSE, # Must sensitivities, specificities and AUC be reported in percent? Note that if TRUE, and you want a partial area, you must pass it in percent also (partial.area=c(100, 80)),
                                   direction = c("auto", "<", ">"),
                                   ...) {
	# We need at least two levels in response
	if (length(unique(response)) < 2) {
		stop("'response' must have at least two levels")
	}
	
    # implements the approach from Hand & Till (2001)
    if (methods::is(predictor, "matrix") || methods::is(predictor, "data.frame")) {
        # for decision values for multiple classes (e.g. probabilities of individual classes)
    	if (missing("direction")) {
            # need to have uni-directional decision values for consistency
            direction <- ">" 
    	}
    	else {
    		direction <- match.arg(direction)
    	}
        mc.roc <- multiclass.roc.multivariate(response, predictor, levels, percent, direction, ...)
    } else {
        # for a single decision value for separating the classes
    	direction <- match.arg(direction)
    	mc.roc <- multiclass.roc.univariate(response, predictor, levels, percent, direction, ...)
    }
	mc.roc$call <- match.call()
	return(mc.roc)
}
