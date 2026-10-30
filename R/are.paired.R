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
                           return.paired.rocs = FALSE,
                           reuse.auc = TRUE, reuse.ci = FALSE, reuse.smooth = TRUE,
                           ...) {
  # check return.paired.rocs
  if (!is.logical(return.paired.rocs) || length(return.paired.rocs) != 1) {
    stop("'return.paired.rocs' must be either TRUE or FALSE.")
  }
  # Recover base ROC curves (not auc or smoothed)
  if ("auc" %in% class(roc1)) {
    roc1 <- attr(roc1, "roc")
  }
  if ("auc" %in% class(roc2)) {
    roc2 <- attr(roc2, "roc")
  }
  oroc1 <- oroc2 <- NULL # the original smoothed curves, if any
  if ("smooth.roc" %in% class(roc1)) {
    oroc1 <- roc1
    roc1 <- attr(roc1, "roc")
  }
  if ("smooth.roc" %in% class(roc2)) {
    oroc2 <- roc2
    roc2 <- attr(roc2, "roc")
  }
  # Check if the levels are the same. Otherwise it is not paired.
  # Compare them as character: levels given explicitly are stored as given
  # (numeric, factor...) while the default levels are character
  if (!identical(as.character(roc1$levels), as.character(roc2$levels))) {
    return(FALSE)
  }
  # check if responses of roc 1 and 2 are identical. Ignore the names (a
  # formula gives the row names, as model.response() does), but not the
  # na.action attribute: two curves that lost different observations must be
  # paired on the observations they share, below.
  if (identical(unname(roc1$response), unname(roc2$response))) {
    retval <- TRUE
    if (return.paired.rocs) {
      # Return the smoothed curves if the input was smoothed
      attr(retval, "roc1") <- if (!is.null(oroc1) && reuse.smooth) oroc1 else roc1
      attr(retval, "roc2") <- if (!is.null(oroc2) && reuse.smooth) oroc2 else roc2
    }
    return(retval)
  } else {
    # Make sure the difference is not only due to missing values
    # compare original response (with NAs and response not in levels)
    if (identical(unname(roc1$original.response), unname(roc2$original.response))) {
      retval <- TRUE
      if (!return.paired.rocs) {
        return(retval)
      }
      # else prepare paired ROCs
      idx.exclude <- is.na(roc1$original.predictor) | is.na(roc2$original.predictor) | is.na(roc1$original.response) | !roc1$original.response %in% roc1$levels
      roc1.paired <- roc(roc1$original.response[!idx.exclude], roc1$original.predictor[!idx.exclude], levels = roc1$levels, percent = roc1$percent, na.rm = FALSE, direction = roc1$direction, auc = FALSE)
      roc2.paired <- roc(roc2$original.response[!idx.exclude], roc2$original.predictor[!idx.exclude], levels = roc2$levels, percent = roc2$percent, na.rm = FALSE, direction = roc2$direction, auc = FALSE)
      # Print the data as the input curves, not as the expressions above
      for (field in c("call", "predictor.name", "response.name")) {
        roc1.paired[[field]] <- roc1[[field]]
        roc2.paired[[field]] <- roc2[[field]]
      }
      # Re-use auc/ci/smooth for roc1
      # Take auc and ci from the smoothed curve if roc1.paired is smoothed
      source.roc1 <- roc1
      if (!is.null(oroc1) && reuse.smooth) {
        source.roc1 <- oroc1
        roc1.paired <- roc_utils_resmooth(roc1.paired, oroc1$smoothing.args)
        roc1.paired$call <- oroc1$call
      }
      if (!is.null(source.roc1$auc) && reuse.auc) {
        args <- attributes(source.roc1$auc)
        args$roc <- roc1.paired
        roc1.paired$auc <- do.call("auc.roc", args)
      }
      if (!is.null(source.roc1$ci) && reuse.ci) {
        args <- attributes(source.roc1$ci)
        args$roc <- NULL
        roc1.paired$ci <- do.call(class(source.roc1$ci)[1], c(list(roc1.paired), args))
      }
      # Re-use auc/ci/smooth for roc2
      # Take auc and ci from the smoothed curve if roc2.paired is smoothed
      source.roc2 <- roc2
      if (!is.null(oroc2) && reuse.smooth) {
        source.roc2 <- oroc2
        roc2.paired <- roc_utils_resmooth(roc2.paired, oroc2$smoothing.args)
        roc2.paired$call <- oroc2$call
      }
      if (!is.null(source.roc2$auc) && reuse.auc) {
        args <- attributes(source.roc2$auc)
        args$roc <- roc2.paired
        roc2.paired$auc <- do.call("auc.roc", args)
      }
      if (!is.null(source.roc2$ci) && reuse.ci) {
        args <- attributes(source.roc2$ci)
        args$roc <- NULL
        roc2.paired$ci <- do.call(class(source.roc2$ci)[1], c(list(roc2.paired), args))
      }

      # Attach ROCs and return value
      attr(retval, "roc1") <- roc1.paired
      attr(retval, "roc2") <- roc2.paired
      return(retval)
    } else {
      return(FALSE)
    }
  }
}

# The pairing decision of roc.test, cov and power.roc.test: paired = NULL
# detects it with are.paired(), TRUE requires it, FALSE treats the curves as
# unpaired. Returns are.paired()'s result (TRUE with the re-paired curves in
# attributes "roc1" and "roc2"), or FALSE.
roc_utils_paired <- function(roc1, roc2, paired) {
  if (!is.null(paired) && !paired) {
    return(FALSE)
  }
  rocs.are.paired <- are.paired(roc1, roc2, return.paired.rocs = TRUE, reuse.auc = TRUE, reuse.ci = FALSE, reuse.smooth = TRUE)
  if (!is.null(paired) && !rocs.are.paired) {
    stop(simpleError("'paired = TRUE' cannot be applied to unpaired ROC curves.", sys.call(-1)))
  }
  rocs.are.paired
}
