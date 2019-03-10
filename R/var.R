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

var <- function(...)
  UseMethod("var")

var.default <- function(...) {
  stats::var(...)
}

var.auc <- function(auc, ...) {
  # Change roc from an auc to a roc object but keep the auc specifications
  roc <- auc
  attr(auc, "roc") <- NULL
  roc <- attr(roc, "roc")
  roc$auc <- auc
  # Pass to var.roc
  var.roc(roc, ...)
}

var.smooth.roc <- function(smooth.roc, ...) {
  var.roc(smooth.roc, ...) # just pass to var.roc that will do the job
}

var.roc <- function(roc,
                    method=c("delong", "bootstrap", "obuchowski"),
                    boot.n = 2000,
                    boot.stratified = TRUE,
                    reuse.auc=TRUE,
                    progress = getOption("pROCProgress")$name,
                    parallel = FALSE,
                    ...) {
  # We need an auc
  if (is.null(roc$auc) | !reuse.auc)
    roc$auc <- auc(roc, ...)
  
  if (roc.utils.is.perfect.curve(roc)) {
  	warning("var() of a ROC curve with AUC == 1 is always 0 and can be misleading.")
  }

  # do all the computations in fraction, re-transform in percent later
  percent <- roc$percent
  if (percent) {
  	roc <- roc.utils.unpercent(roc)
  }

  # Check the method
  if (missing(method) | is.null(method)) {
    # determine method if missing
    if (has.partial.auc(roc)) {
      # partial auc: go for bootstrap
      method <- "bootstrap"
    }
    else if (class(roc) == "smooth.roc") {
      # smoothing: bootstrap
      method <- "bootstrap"
    }
    else {
      method <- "delong"
    }
  }
  else {
    method <- match.arg(method)
    # delong NA to pAUC: warn + change
    if (method == "delong") {
      if (has.partial.auc(roc)) {
      	stop("DeLong method is not supported for partial AUC. Use method=\"bootstrap\" instead.")
      }
      else if ("smooth.roc" %in% class(roc)) {
      	stop("DeLong method is not supported for smoothed ROCs. Use method=\"bootstrap\" instead.")
      }
    }

    else if (method == "obuchowski") {
      if ("smooth.roc" %in% class(roc)) {
        stop("Using Obuchowski for smoothed ROCs is not supported. Using bootstrap instead.")
      }
      if (has.partial.auc(roc) && attr(roc$auc, "partial.auc.focus") == "sensitivity") {
        stop("Using Obuchowski for partial AUC on sensitivity region is not supported. Using bootstrap instead.")
      }
    }
  }
  
  if (method == "delong") {
    n <- length(roc$controls)
    m <- length(roc$cases)  
    V <- delongPlacements(roc)
    var <- var(V$Y) / n + var(V$X) / m
  }
  else if (method == "obuchowski") {
    var <- var.roc.obuchowski(roc) / length(roc$cases)
  }
  else {
    var <- var.roc.bootstrap(roc, boot.n, boot.stratified, progress, parallel, ...)
  }
  
  if (percent) {
    var <- var * (100^2)
  }
  return(var)
}

var.roc.bootstrap <- function(roc, boot.n, boot.stratified, progress, parallel, ...) {
  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="AUC variance", label="Bootstrap in progress...", ...)

  ## Smoothed ROC curve variance
  if (class(roc) == "smooth.roc") {
    smoothing.args <- roc$smoothing.args
    smoothing.args$smooth <- TRUE
    non.smoothed.roc <- attr(roc, "roc")
    non.smoothed.roc$percent <- FALSE # as we did earlier for the smoothed.roc
    smooth.roc.call <- as.call(c(utils::getS3method("smooth", "roc"), roc$smoothing.args))
    auc.args <- attributes(roc$auc)[grep("partial.auc", names(attributes(roc$auc)))]
    auc.args$allow.invalid.partial.auc.correct <- TRUE
    auc.call <- as.call(c(utils::getS3method("auc", "smooth.roc"), auc.args))

    if (boot.stratified) {
      aucs <- unlist(llply(1:boot.n, stratified.ci.smooth.auc, roc=non.smoothed.roc, smooth.roc.call=smooth.roc.call, auc.call=auc.call, .progress=progress, .parallel=parallel))
    }
    else {
      aucs <- unlist(llply(1:boot.n, nonstratified.ci.smooth.auc, roc=non.smoothed.roc, smooth.roc.call=smooth.roc.call, auc.call=auc.call, .progress=progress, .parallel=parallel))
    }
  }
  ## Non smoothed ROC curves variance
  else {
    if (boot.stratified) {
      aucs <- unlist(llply(1:boot.n, stratified.ci.auc, roc=roc, .progress=progress, .parallel=parallel)) # ci.auc: returns aucs just as we need for var, so re-use it!
    }
    else {
      aucs <- unlist(llply(1:boot.n, nonstratified.ci.auc, roc=roc, .progress=progress, .parallel=parallel))
    }
  }

  if ((num.NAs <- sum(is.na(aucs))) > 0) {
    warning(sprintf("%i NA value(s) produced during bootstrap were ignored.", num.NAs))
    aucs <- aucs[!is.na(aucs)]
  }
  return(var(aucs))
}
