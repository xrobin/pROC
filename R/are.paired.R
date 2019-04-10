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
  # Check if the levels are the same. Otherwise it is not paired.
  if (!identical(roc1$levels, roc2$levels))
    return(FALSE)
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
    # compare original response (with NAs and response not in levels)
    if (identical(roc1$original.response, roc2$original.response)) {
      retval <- TRUE
      if (! return.paired.rocs)
        return(retval)
      # else prepare paired ROCs
      idx.exclude <- is.na(roc1$original.predictor) | is.na(roc2$original.predictor) | is.na(roc1$original.response) | ! roc1$original.response %in% roc1$levels
      roc1.paired <- roc(roc1$original.response[!idx.exclude], roc1$original.predictor[!idx.exclude], levels=roc1$levels, percent=roc1$percent, na.rm=FALSE, direction=roc1$direction, auc=FALSE)
      roc2.paired <- roc(roc2$original.response[!idx.exclude], roc2$original.predictor[!idx.exclude], levels=roc2$levels, percent=roc2$percent, na.rm=FALSE, direction=roc2$direction, auc=FALSE)
      # Re-use auc/ci/smooth for roc1
      if (exists("oroc1") && reuse.smooth) {
        args <- oroc1$smoothing.args
        args$roc <- roc1.paired
        roc1.paired <- do.call("smooth.roc", args)
        roc1.paired$call$roc <- as.name("roc1.paired")
      }
      if (!is.null(roc1$auc) && reuse.auc) {
        args <- attributes(roc1$auc)
        args$roc <- roc1.paired
        roc1.paired$auc <- do.call("auc.roc", args)
      }
      if (!is.null(roc1$ci) && reuse.ci) {
        args <- attributes(roc1$ci)
        args$roc <- NULL
        roc1.paired$ci <- do.call(class(roc1$ci)[1], c(roc=list(roc1.paired), args))
      }
      # Re-use auc/ci/smooth for roc2
      if (exists("oroc2") && reuse.smooth) {
        args <- oroc2$smoothing.args
        args$roc <- roc2.paired
        roc2.paired <- do.call("smooth.roc", args)
        roc2.paired$call$roc <- as.name("roc2.paired")
      }
      if (!is.null(roc2$auc) && reuse.auc) {
        args <- attributes(roc2$auc)
        args$roc <- roc2.paired
        roc2.paired$auc <- do.call("auc.roc", args)
      }
      if (!is.null(roc2$ci) && reuse.ci) {
        args <- attributes(roc2$ci)
        args$roc <- NULL
        roc2.paired$ci <- do.call(class(roc2$ci)[1], c(roc=list(roc2.paired), args))
      }

      # Attach ROCs and return value
      attr(retval, "roc1") <- roc1.paired
      attr(retval, "roc2") <- roc2.paired
      return(retval)
    }
    else {
      return(FALSE)
    }
  }
}
