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

multiclass.roc <- function(...)
  UseMethod("multiclass.roc")

multiclass.roc.formula <- function(formula, data, ...) {
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  term.labels <- attr(attr(m, "terms"), "term.labels")
  response <- model.extract(m, "response")

  if (length(response) == 0) {
    stop("Error in the formula: a response is required in a formula of type response~predictor.")
  }

  multiclass.roc <- multiclass.roc.default(response, m[[term.labels]], ...)
  multiclass.roc$call <- match.call()
  return(multiclass.roc)
}

multiclass.roc.univariate <- function(response, predictor,
                                   levels=base::levels(as.factor(response)),
                                   percent=FALSE, # Must sensitivities, specificities and AUC be reported in percent? Note that if TRUE, and you want a partial area, you must pass it in percent also (partial.area=c(100, 80))
                                   # what computation must be done
                                   #auc=TRUE, # call auc.roc on the current object
                                   #ci=FALSE, # call ci.roc on the current object
                                   ...) {
  multiclass.roc <- list(
                         response = response,
                         predictor = predictor,
                         percent = percent,
                         call=match.call())
  class(multiclass.roc) <- "multiclass.roc"

  if ("ordered" %in% class(response) && any(names(table(response))[table(response) == 0] %in% levels)) {
    missing.levels <- names(table(response))[table(response) == 0]
    missing.levels.requested <- missing.levels[missing.levels %in% levels]
    warning(paste("No observation for response level(s):", paste(missing.levels.requested, collapse=", ")))
    levels <- levels[!(levels %in% missing.levels.requested)]
  }
  multiclass.roc$levels <- levels
  
  rocs <- combn(levels, 2, function(X, response, predictor, percent, ...) {
    roc(response, predictor, levels=X, percent=percent, auc=FALSE, ci=FALSE, ...)
  }, simplify=FALSE, response=response, predictor=predictor, percent=percent, ...)

  multiclass.roc$rocs <- rocs

  # Makes no sense to turn auc off, so remove this option
  #if (auc)
    multiclass.roc$auc <- auc.multiclass.roc(multiclass.roc, ...)
  # CI is not implemented yet.
  #if (ci)
  #  multiclass.roc$ci <- ci.multiclass.roc(multiclass.roc, ...)

  return(multiclass.roc)
}

compute.pair.AUC <- function(pred.matrix, i, j, ref.outcome, levels, percent, ... ) {
    # computes A(i|j), the probability that a randomly 
    # chosen member of class j has a lower estimated probability (or score) 
    # of belonging to class i than a randomly chosen member of class i

    pred.i <- pred.matrix[which(ref.outcome == i), i] # p(G = i) assigned to class i observations
    pred.j <- pred.matrix[which(ref.outcome == j), i] # p(G = i) assigned to class j observations
    classes <- factor(c(rep(i, length(pred.i)), rep(j, length(pred.j))))
    # override levels argument by new levels
    levels <- unique(classes)
    predictor <- c(pred.i, pred.j)
    auc <- roc(classes, predictor, levels = levels, percent = percent, auc = FALSE, ci = FALSE, ...)
    return(auc)
}

compute.new.AUC <- function(pred.matrix, i, j, ref.outcome, levels, percent, ... ) {
    # computes A(i|j), the probability that a randomly 
    # chosen member of class j has a lower estimated probability (or score) 
    # of belonging to class i than a randomly chosen member of class i

    pred.i <- pred.matrix[which(ref.outcome == i), i] # p(G = i) assigned to class i observations
    pred.j <- pred.matrix[which(ref.outcome == j), i] # p(G = i) assigned to class j observations
    classes <- factor(c(rep(i, length(pred.i)), rep(j, length(pred.j))))
    # override levels argument by new levels
    levels <- unique(classes)
    predictor <- c(pred.i, pred.j)
    auc <- roc(classes, predictor, levels = levels, percent = percent, auc = FALSE, ci = FALSE, ...)
    return(auc)
}

multiclass.roc.multivariate <- function(response, predictor, levels, percent, ...) {
    # Reference: "A Simple Generalisation of the Area Under the ROC 
    # Curve for Multiple Class Classification Problems" (Hand and Till, 2001)
    if (!is(predictor, "matrix")) {
        stop("Please provide a matrix via 'predictor'.")
    }
    if (nrow(predictor) != length(response)) {
        stop("Number of rows in 'predictor' does not agree with 'response'");
    }
    # check whether the columns of the prediction matrix agree with the factors in 'response'
    m <- match(colnames(predictor), levels)
    levels <- colnames(predictor)[!is.na(m)]
    if (length(levels) == 1) {
        stop("For a single decision value, please provide 'predictor' as a vector.")
    } else if (length(levels) == 0) {
        stop("The column names of 'predictor' could not be matched to the levels of 'response'.")
    }
    missing.classes <- levels[setdiff(seq_along(levels), m)]
    if (length(missing.classes) != 0) {
        out.classes <- paste0(missing.classes, collapse = ",")
        if (length(missing.classes) == length(levels)) {
            # no decision values found
            stop("Could not find any decision values in 'predictor' matching the 'response' levels. Could not find the following classes: ", out.classes, ". Check your column names!")
        } else {
            # some decision values found
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
                         percent = percent,
                         call=match.call())
    class(multiclass.roc) <- "multiclass.roc"
    multiclass.roc$levels <- levels
    rocs <- utils::combn(levels, 2, function(x, predictor, response, levels, percent, ...) {
        A1 <- compute.pair.AUC(predictor, x[1], x[2], response, levels, percent, ...)
        A2 <- compute.pair.AUC(predictor, x[2], x[1], response, levels, percent, ...)
        #print(auc(A1))
        #print(auc(A2))
        # merging A1 and A2 is infeasible as auc() would not be well-defined
        A <- list(A1, A2) 
        return(A)
    }, simplify = FALSE, predictor = predictor, response = response, 
        levels = levels, percent = percent, ...)
    pairs <- unlist(lapply(combn(levels, 2, simplify = FALSE), 
                           function(x) paste(x, collapse = "/")))
    names(rocs) <- pairs
    multiclass.roc$rocs <- rocs
    # can't use 'auc.multiclass.roc' because mean of ROCs needs to be computed before ...
    #A.total <- auc.multiclass.roc(multiclass.roc)
    mean.AUCs <- unlist(lapply(rocs, function(x) mean(auc(x[[1]]), auc(x[[2]]))))
    A.ij.total <- sum(mean.AUCs)
    c <- length(levels)
    M <- 2 / (c * (c-1)) * A.ij.total 
    multiclass.roc$auc <- M # manually store auc; TODO: improve functioning of auc function
    multiclass.roc$pair_aucs <- mean.AUCs # manually store pairwise AUCs
    return(multiclass.roc)
}

multiclass.roc.default <- function(response, predictor,
                                   levels=base::levels(as.factor(response)),
                                   percent=FALSE, # Must sensitivities, specificities and AUC be reported in percent? Note that if TRUE, and you want a partial area, you must pass it in percent also (partial.area=c(100, 80))
                                   # what computation must be done
                                   #auc=TRUE, # call auc.roc on the current object
                                   #ci=FALSE, # call ci.roc on the current object
                                   ...) {
    # implements the approach from Hand & Till (2001)

    if (is(predictor, "matrix")) {
        # for decision values for multiple classes (e.g. probabilities of individual classes)
        return(multiclass.roc.multivariate(response, predictor, levels, percent, ...))
    } else {
        # for a single decision value for separating the classes
        return(multiclass.roc.univariate(response, predictor, levels, percent, ...))
    }
}

# Previous code included for testing: results are different!

compute.A.conditional <- function(pred.matrix, i, j, ref.outcome) {
    # computes A(i|j), the probability that a randomly 
    # chosen member of class j has a lower estimated probability (or score) 
    # of belonging to class i than a randomly chosen member of class i

    i.idx <- which(ref.outcome == i)
    j.idx <- which(ref.outcome == j)
    pred.i <- pred.matrix[i.idx, i] # p(G = i) assigned to class i observations
    pred.j <- pred.matrix[j.idx, i] # p(G = i) assigned to class j observations
    all.preds <- c(pred.i, pred.j)
    classes <- c(rep(i, length(pred.i)), rep(j, length(pred.j)))
    o <- order(all.preds)
    classes.o <- classes[o]
    # Si: sum of ranks from class i observations
    Si <- sum(which(classes.o == i))
    ni <- length(i.idx)
    nj <- length(j.idx)
    # calculate A(i|j)
    A <- (Si - ((ni * (ni + 1))/2)) / (ni * nj)
    return(A)
}

multiclass.roc.multivariate.custom <- function(response, predictor) {
    # Reference: "A Simple Generalisation of the Area Under the ROC 
    # Curve for Multiple Class Classification Problems" (Hand and Till, 2001)
    if (!is(predictor, "matrix")) {
        stop("Please provide a matrix via 'predictor'.")
    }
    if (nrow(predictor) != length(response)) {
        stop("Number of rows in 'predictor' does not agree with 'response'");
    }
    # check whether the columns of the prediction matrix agree with the factors in 'response'
    m <- match(colnames(predictor), levels(response))
    labels <- colnames(predictor)[!is.na(m)]
    if (length(labels) == 1) {
        stop("For a single decision value, please provide 'predictor' as a vector.")
    } else if (length(labels) == 0) {
        stop("The column names of 'predictor' could not be matched to the levels of 'response'.")
    }
    missing.classes <- levels(response)[setdiff(seq_along(levels(response)), m)]
    if (length(missing.classes) != 0) {
        out.classes <- paste0(missing.classes, collapse = ",")
        if (length(missing.classes) == length(levels(response))) {
            # no decision values found
            stop("Could not find any decision values in 'predictor' matching the 'response' levels. Could not find the following classes: ", out.classes, ". Check your column names!")
        } else {
            # some decision values found
            warning("You did not provide decision values for the following classes: ", out.classes, ".")
        }
    }
    additional.classes <- colnames(predictor)[which(is.na(m))]
    if (length(additional.classes) != 0) {
        out.classes <- paste0(additional.classes, collapse = ",")
        warning("The following classes were not found in 'response': ", out.classes, ".")
    }
    A.ij.cond <- utils::combn(labels, 2, function(x, predictor, response) {x
        i <- x[1]
        j <- x[2]
        A.ij <- compute.A.conditional(predictor, i, j, response)
        A.ji <- compute.A.conditional(predictor, j, i, response)
        pair <- paste0(i, "/", j)
        return(c(A.ij, A.ji))
    }, simplify = FALSE, predictor = predictor, response = response)
    c <- length(labels)
    pairs <- unlist(lapply(combn(labels, 2, simplify = FALSE), function(x) paste(x, collapse = "/")))
    A.ij.joint <- unlist(lapply(A.ij.cond, mean)) # <=> A(i,j)
    names(A.ij.joint) <- pairs
    A.ij.total <- sum(unlist(A.ij.joint))
    M <- 2 / (c * (c-1)) * A.ij.total 
    attr(M, "pair_AUCs") <- A.ij.joint
    return(M)
}
