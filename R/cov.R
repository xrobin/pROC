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

cov <- function(...) {
  UseMethod("cov")
}

cov.default <- function(...) {
  stats::cov(...)
}

cov.auc <- function(roc1, roc2, ...) {
  # Change roc1 from an auc to a roc object but keep the auc specifications
  auc1 <- roc1
  attr(auc1, "roc") <- NULL
  roc1 <- attr(roc1, "roc")
  roc1$auc <- auc1
  # Pass to cov.roc
  return(cov.roc(roc1, roc2, ...))
}

cov.smooth.roc <- function(roc1, roc2, ...) {
  cov.roc(roc1, roc2, ...)
}

cov.roc <- function(roc1, roc2,
                         method=c("delong", "bootstrap", "obuchowski"),
                         reuse.auc=TRUE,
                         boot.n=2000, boot.stratified=TRUE, boot.return=FALSE,
                         progress=getOption("pROCProgress")$name,
                         parallel = FALSE,
                         ...) {
  # If roc2 is an auc, take the roc but keep the auc specifications
  if (methods::is(roc2, "auc")) {
    auc2 <- roc2
    attr(auc2, "roc") <- NULL
    roc2 <- attr(roc2, "roc")
    roc2$auc <- auc2
  }
  
  if (roc.utils.is.perfect.curve(roc1) && roc.utils.is.perfect.curve(roc2)) {
  	warning("cov() of two ROC curves with AUC == 1 is always 0 and can be misleading.")
  }

  # store which objects are smoothed, and how
  smoothing.args <- list()
  if ("smooth.roc" %in% class(roc1)) {
    smoothing.args$roc1 <- roc1$smoothing.args
    smoothing.args$roc1$smooth <- TRUE
    roc1 <- attr(roc1, "roc")
    #oroc1$auc <- roc1$auc
  }
  else {
    smoothing.args$roc1 <- list(smooth=FALSE)
  }
  if ("smooth.roc" %in% class(roc2)) {
    smoothing.args$roc2 <- roc2$smoothing.args
    smoothing.args$roc2$smooth <- TRUE
    roc2 <- attr(roc2, "roc")
    #oroc2$auc <- roc2$auc
  }
  else {
    smoothing.args$roc2 <- list(smooth=FALSE)
  }

  # then determine whether the rocs are paired or not
  rocs.are.paired <- are.paired(roc1, roc2, return.paired.rocs=FALSE, reuse.auc=TRUE, reuse.ci=FALSE, reuse.smooth=TRUE)
  if (! rocs.are.paired) {
    message("ROC curves are unpaired.")
    return(0)
  }    

  # check that the AUC was computed, or do it now
  if (is.null(roc1$auc) | !reuse.auc) {
    if (smoothing.args$roc1$smooth) {
      roc1$auc <- auc(smooth.roc=do.call("smooth.roc", c(list(roc=roc1), smoothing.args$roc1)), ...)
      # remove partial.auc.* arguments that are now in roc1$auc and that will mess later processing
      # (formal argument "partial.auc(.*)" matched by multiple actual arguments)
      # This removal should be safe because we always use smoothing.args with roc1 in the following processing,
      # however it is a potential source of bugs.
      smoothing.args$roc1$partial.auc <- NULL
      smoothing.args$roc1$partial.auc.correct <- NULL
      smoothing.args$roc1$partial.auc.focus <- NULL
    }
    else
      roc1$auc <- auc(roc1, ...)
  }
  if (is.null(roc2$auc) | !reuse.auc) {
    if (smoothing.args$roc2$smooth) {
      roc2$auc <- auc(smooth.roc=do.call("smooth.roc", c(list(roc=roc2), smoothing.args$roc2)), ...)
      # remove partial.auc.* arguments that are now in roc1$auc and that will mess later processing
      # (formal argument "partial.auc(.*)" matched by multiple actual arguments)
      # This removal should be safe because we always use smoothing.args with roc2 in the following processing,
      # however it is a potential source of bugs.
      smoothing.args$roc2$partial.auc <- NULL
      smoothing.args$roc2$partial.auc.correct <- NULL
      smoothing.args$roc2$partial.auc.focus <- NULL
    }
    else
      roc2$auc <- auc(roc2, ...)
  }
    
  # check that the same region was requested in auc. Otherwise, issue a warning
  if (!identical(attributes(roc1$auc)[names(attributes(roc1$auc))!="roc"], attributes(roc2$auc)[names(attributes(roc2$auc))!="roc"]))
    warning("Different AUC specifications in the ROC curves. Enforcing the inconsistency, but unexpected results may be produced.")
  # check that the same smoothing params were requested in auc. Otherwise, issue a warning
  if (!identical(smoothing.args$roc1, smoothing.args$roc2))
    warning("Different smoothing parameters in the ROC curves. Enforcing the inconsistency, but unexpected results may be produced.")

  # Check the method
  if (missing(method) | is.null(method)) {
    # determine method if missing
    if (has.partial.auc(roc1)) {
      # partial auc: go for bootstrap
      method <- "bootstrap"
    }
    else if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth) {
      # smoothing in one or both: bootstrap
      method <- "bootstrap"
    }
    else if (roc1$direction != roc2$direction) {
      # delong doesn't work well with opposite directions (will report high significance if roc1$auc and roc2$auc are similar and high)
      method <- "bootstrap"
    }
    else {
      method <- "delong"
    }
  }
  else {
    method <- match.arg(method)
    if (method == "delong") {
      # delong NA to pAUC: warn + change
      if (has.partial.auc(roc1) || has.partial.auc(roc2)) {
      	stop("DeLong method is not supported for partial AUC. Use method=\"bootstrap\" instead.")
      }
      if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth) {
      	stop("DeLong method is not supported for smoothed ROCs. Use method=\"bootstrap\" instead.")
      }
      if (roc1$direction != roc2$direction)
        warning("DeLong method should not be applied to ROC curves with a different direction.")
    }
    else if (method == "obuchowski") {
      if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth) {
        stop("Obuchowski method is not supported for smoothed ROCs. Use method=\"bootstrap\" instead.")
      }
      if ((has.partial.auc(roc1) && attr(roc1$auc, "partial.auc.focus") == "sensitivity")
          || (has.partial.auc(roc2) && attr(roc2$auc, "partial.auc.focus") == "sensitivity")) {
        stop("Obuchowski method is not supported for partial AUC on sensitivity region. Use method=\"bootstrap\" instead.")
      }
      if (roc1$direction != roc2$direction)
        warning("Obuchowski method should not be applied to ROC curves with a different direction.")
    }
  }
  
  if (method == "delong") {
    n <- length(roc1$controls)
    m <- length(roc1$cases)

    V1 <- delongPlacements(roc1)
    var1 <- var(V1$Y) / n + var(V1$X) / m

    V2 <- delongPlacements(roc2)
    var2 <- var(V2$Y) / n + var(V2$X) / m

    cov <- cov(V2$X, V1$X) / m + cov(V2$Y, V1$Y) / n

    if (roc1$percent) {
      cov <- cov * (100^2)
    }
  }
  
  else if (method == "obuchowski") {
    cov <- cov.roc.obuchowski(roc1, roc2) / length(roc1$cases)

    if (roc1$percent) {
      cov <- cov * (100^2)
    }
  }
  else { # method == "bootstrap"
    # Check if called with density.cases or density.controls
    if (is.null(smoothing.args) || is.numeric(smoothing.args$density.cases) || is.numeric(smoothing.args$density.controls))
      stop("Cannot compute the covariance of ROC curves smoothed with numeric density.controls and density.cases.")

    if(class(progress) != "list") {
      progress <- roc.utils.get.progress.bar(progress, title="Bootstrap covariance", label="Bootstrap in progress...", ...)
    }
    
    cov <- bootstrap.cov(roc1, roc2, boot.n, boot.stratified, boot.return, smoothing.args, progress, parallel)
  }

  return(cov)
}
