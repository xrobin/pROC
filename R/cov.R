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

cov <- function(...) {
  UseMethod("cov")
}

cov.default <- function(...) {
  stats::cov(...)
}

cov.auc <- function(roc1, roc2, ...) {
  return(cov.roc(attr(roc1, "roc"), roc2, ...))
}

cov.smooth.roc <- function(roc1, roc2, ...) {
  cov.roc(roc1, roc2, ...)
}

cov.roc <- function(roc1, roc2,
                         method=c("delong", "bootstrap"),
                         reuse.auc=TRUE,
                         boot.n=2000, boot.stratified=TRUE, boot.return=FALSE,
                         progress=getOption("pROCProgress")$name,
                         ...) {
  if ("auc" %in% class(roc2))
    roc2 <- attr(roc2, "roc")

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
    if (is.numeric(attr(roc1$auc, "partial.auc")) && length(attr(roc1$auc, "partial.auc") == 2)) {
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
      if ((is.numeric(attr(roc1$auc, "partial.auc")) && length(attr(roc1$auc, "partial.auc") == 2)) || (is.numeric(attr(roc2$auc, "partial.auc")) && length(attr(roc2$auc, "partial.auc") == 2))) {
        warning("Using DeLong for partial AUC is not supported. Using bootstrap instead.")
        method <- "bootstrap"
      }
      if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth) {
        warning("Using DeLong for smoothed ROCs is not supported. Using bootstrap instead.")
        method <- "bootstrap"
      }
      if (roc1$direction != roc2$direction)
        warning("DeLong should not be applied to ROC curves with a different direction.")
    }
  }
  
  if (method == "delong") {
    n <- length(roc1$controls)
    m <- length(roc1$cases)

    V1 <- roc.utils.delong.placements(roc1)
    var1 <- var(V1$Y) / n + var(V1$X) / m

    V2 <- roc.utils.delong.placements(roc2)
    var2 <- var(V2$Y) / n + var(V2$X) / m

    cov <- cov(V2$X, V1$X) / m + cov(V2$Y, V1$Y) / n

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

    cov <- bootstrap.cov(roc1, roc2, boot.n, boot.stratified, boot.return, smoothing.args, progress)
  }

  return(cov)
}

bootstrap.cov <- function(roc1, roc2, boot.n, boot.stratified, boot.return, smoothing.args, progress) {

  # rename method into smooth.method for roc
  smoothing.args$roc1$smooth.method <- smoothing.args$roc1$method
  smoothing.args$roc1$method <- NULL
  smoothing.args$roc2$smooth.method <- smoothing.args$roc2$method
  smoothing.args$roc2$method <- NULL

  # Prepare arguments for later calls to roc
  auc1skeleton <- attributes(roc1$auc)
  auc1skeleton$roc <- NULL
  auc1skeleton$direction <- roc1$direction
  auc1skeleton$class <- NULL
  auc1skeleton <- c(auc1skeleton, smoothing.args$roc1)
  auc2skeleton <- attributes(roc2$auc)
  auc2skeleton$roc <- NULL
  auc2skeleton$direction <- roc2$direction
  auc2skeleton$class <- NULL
  auc2skeleton <- c(auc2skeleton, smoothing.args$roc2)

  # Some attributes may be duplicated in AUC skeletons and will mess the boostrap later on when we do.call().
  # If this condition happen, it probably means we have a bug elsewhere.
  # Rather than making a complicated processing to remove the duplicates,
  # just throw an error and let us solve the bug when a user reports it.
  duplicated.auc1skeleton <- duplicated(names(auc1skeleton))
  duplicated.auc2skeleton <- duplicated(names(auc2skeleton))
  if (any(duplicated.auc1skeleton))
    stop(sprintf("duplicated argument(s) in AUC1 skeleton: \"%s\". Please report this bug to the package maintainer %s", paste(names(auc1skeleton)[duplicated(names(auc1skeleton))], collapse=", "), packageDescription("pROC")$Maintainer))
  if (any(duplicated.auc2skeleton))
    stop(sprintf("duplicated argument(s) in AUC2 skeleton: \"%s\". Please report this bug to the package maintainer %s", paste(names(auc2skeleton)[duplicated(names(auc2skeleton))], collapse=", "), packageDescription("pROC")$Maintainer))

  if (boot.stratified) { # precompute sorted responses if stratified
    response.roc1 <- factor(c(rep(roc1$levels[1], length(roc1$controls)), rep(roc1$levels[2], length(roc1$cases))), levels=roc1$levels)
    response.roc2 <- factor(c(rep(roc2$levels[1], length(roc2$controls)), rep(roc2$levels[2], length(roc2$cases))), levels=roc2$levels)
    auc1skeleton$response <- response.roc1
    auc2skeleton$response <- response.roc2
    resampled.values <- raply(boot.n, stratified.bootstrap.test(roc1, roc2, "boot", NULL, TRUE, auc1skeleton, auc2skeleton), .progress=progress) # stratified.bootstrap.test: returns resampled values just as we need for cov, so re-use it!
  }
  else {
    resampled.values <- raply(boot.n, nonstratified.bootstrap.test(roc1, roc2, "boot", NULL, TRUE, auc1skeleton, auc2skeleton), .progress=progress)
  }

  # are there NA values?
  if ((num.NAs <- sum(apply(resampled.values, 1, is.na))) > 0) {
    warning(sprintf("%i NA value(s) produced during bootstrap were ignored.", num.NAs))
    resampled.values <- resampled.values[!apply(resampled.values, 1, function(x) any(is.na(x))),]
  }

  cov <- stats::cov(resampled.values[,1], resampled.values[,2])
  if (boot.return) {
    attr(cov, "resampled.values") <- resampled.values
  }
  return(cov)
}
