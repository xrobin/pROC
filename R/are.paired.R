# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2010, 2011 Xavier Robin, Alexandre Hainard, Natacha Turck,
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

are.paired <- function(...) {
  UseMethod("are.paired")
}

are.paired.auc <- function(roc1, roc2, ...) {
  return(are.paired.roc(roc1, roc2, ...))
}

are.paired.smooth.roc <- function(roc1, roc2, ...) {
  return(are.paired.roc(roc1, roc2, ...))
}

are.paired.roc <- function(roc1, roc2,
                           return.paired.rocs=FALSE,
                           reuse.auc = TRUE, reuse.ci = FALSE, reuse.smooth=TRUE,
                           ...) {
  # check return.paired.rocs
  if (! is.logical(return.paired.rocs) || length(return.paired.rocs) != 1)
    stop("'return.paired.rocs' must be either TRUE or FALSE.")
  # Recover base ROC curves (not auc or smoothed)
  if ("auc" %in% class(roc1))
    roc1 <- attr(roc1, "roc")
  if ("auc" %in% class(roc2))
    roc2 <- attr(roc2, "roc")
  if ("smooth.roc" %in% class(roc1)) {
    oroc1 <- roc1
    roc1 <- attr(roc1, "roc")
  }
  if ("smooth.roc" %in% class(roc2)) {
    oroc2 <- roc2
    roc2 <- attr(roc2, "roc")
  }
  # check if responses of roc 1 and 2 are identical
  if (identical(roc1$response, roc2$response)) {
    retval <- TRUE
    if (return.paired.rocs) {
      attr(retval, "roc1") <- roc1
      attr(retval, "roc2") <- roc2
    }
    return(retval)
  }
  else {
    # Make sure the difference is not only due to missing values
    # restore original response (with NAs)
    res1 <- roc1$response
    if (any(attr(res1, "na.action"))) {
      res1.with.nas <- rep(NA, length(res1) + length(attr(res1, "na.action")))
      res1.with.nas[-attr(res1, "na.action")] <- res1
      attributes(res1.with.nas) <- attributes(res1)
      attr(res1.with.nas, "na.action") <- NULL
    }
    else {
      res1.with.nas <- roc1$response
    }
    res2 <- roc2$response
    if (any(attr(res2, "na.action"))) {
      res2.with.nas <- rep(NA, length(res2) + length(attr(res2, "na.action")))
      res2.with.nas[-attr(res2, "na.action")] <- res2
      attributes(res2.with.nas) <- attributes(res2)
      attr(res2.with.nas, "na.action") <- NULL
    }
    else {
      res2.with.nas <- roc2$response
    }
    # if lengthes with NAs restored differ, not paired
    if (length(res1.with.nas) != length(res2.with.nas))
      return(FALSE)
    # re-remove NAs, but on both responses
    nas <- is.na(res1.with.nas) | is.na(res2.with.nas)
    res1.without.nas <- res1.with.nas[!nas]
    res2.without.nas <- res2.with.nas[!nas]
    if (! return.paired.rocs) {
      return(identical(res1.without.nas, res2.without.nas))
    }
    else { # return.paired.rocs == TRUE, so we must return the paired rocs
      # re-check identity
      if (identical(res1.without.nas, res2.without.nas)) {
        retval <- TRUE
        # If identical, we still need to remove the NAs and re-compute the ROC curves
        if (any(attr(res1, "na.action"))) {
          predictor1.with.nas <- rep(NA, length(res1.with.nas))
          predictor1.with.nas[-attr(res1, "na.action")] <- roc1$predictor
          attributes(predictor1.with.nas) <- attributes(roc1$predictor)
          attr(predictor1.with.nas, "na.action") <- NULL
        }
        else {
          predictor1.with.nas <- roc1$predictor
        }
        roc1.without.nas <- roc(res1.without.nas, predictor1.with.nas[!nas], levels=roc1$levels, percent=roc1$percent, na.rm=FALSE, direction=roc1$direction, auc=FALSE)
        
        if (any(attr(res2, "na.action"))) {
          predictor2.with.nas <- rep(NA, length(res2.with.nas))
          predictor2.with.nas[-attr(res2, "na.action")] <- roc2$predictor
          attributes(predictor2.with.nas) <- attributes(roc2$predictor)
          attr(predictor2.with.nas, "na.action") <- NULL
        }
        else {
          predictor2.with.nas <- roc2$predictor
        }
        roc2.without.nas <- roc(res2.with.nas[!nas], predictor2.with.nas[!nas], levels=roc2$levels, percent=roc2$percent, na.rm=FALSE, direction=roc2$direction, auc=FALSE)

        # Re-use auc/ci/smooth for roc1
        if (exists("oroc1") && reuse.smooth) {
          print("Smoothing ROC1")
          args <- oroc1$smoothing.args
          args$roc <- roc1.without.nas
          roc1.without.nas <- do.call("smooth.roc", args)
          roc1.without.nas$call$roc <- as.name("roc1.without.nas")
        }
        if (!is.null(roc1$auc) && reuse.auc) {
          args <- attributes(roc1$auc)
          args$roc <- roc1.without.nas
          roc1.without.nas$auc <- do.call("auc.roc", args)
        }
        if (!is.null(roc1$ci) && reuse.ci) {
          args <- attributes(roc1$ci)
          args$roc <- NULL
          roc1.without.nas$ci <- do.call(class(roc1$ci), c(roc=list(roc1.without.nas), args))
        }
        # Re-use auc/ci/smooth for roc2
        if (exists("oroc2") && reuse.smooth) {
          args <- roc2$smoothing.args
          args$roc <- roc2.without.nas
          roc2.without.nas <- do.call("smooth.roc", args)
          roc2.without.nas$call$roc <- as.name("roc2.without.nas")
        }
        if (!is.null(roc2$auc) && reuse.auc) {
          args <- attributes(roc2$auc)
          args$roc <- roc2.without.nas
          roc2.without.nas$auc <- do.call("auc.roc", args)
        }
        if (!is.null(roc2$ci) && reuse.ci) {
          args <- attributes(roc2$ci)
          args$roc <- NULL
          roc2.without.nas$ci <- do.call(class(roc2$ci), c(roc=list(roc2.without.nas), args))
        }

        # attach the paired rocs and return
        attr(retval, "roc1") <- roc1.without.nas
        attr(retval, "roc2") <- roc2.without.nas
        return(retval)
      }
      else { # ! identical(res1.without.nas, res2.without.nas)
        return(FALSE)
      }
    }
  }
}
