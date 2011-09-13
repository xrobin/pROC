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

var <- function(...)
  UseMethod("var")

var.default <- function(...) {
  stats::var(...)
}

var.auc <- function(auc, ...) {
  var.roc(attr(auc, "roc"), ...)
}

var.smooth.roc <- function(smooth.roc, ...) {
  var.roc(smooth.roc, ...) # just pass to var.roc that will do the job
}

var.roc <- function(roc,
                    method=c("delong", "bootstrap"),
                    boot.n = 2000,
                    boot.stratified = TRUE,
                    reuse.auc=TRUE,
                    progress = getOption("pROCProgress")$name,
                    ...) {
  # We need an auc
  if (is.null(roc$auc) | !reuse.auc)
    roc$auc <- auc(roc, ...)

  # do all the computations in fraction, re-transform in percent later
  percent <- roc$percent
  roc$percent <- FALSE
  oldauc <- roc$auc
  if (percent) {
    attr(roc$auc, "percent") <- FALSE
    if (! identical(attr(roc$auc, "partial.auc"), FALSE))
      attr(roc$auc, "partial.auc") <- attr(roc$auc, "partial.auc") / 100
  }

  # Check the method
  if (missing(method) | is.null(method)) {
    # determine method if missing
    if (is.numeric(attr(roc$auc, "partial.auc")) && length(attr(roc$auc, "partial.auc") == 2)) {
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
    method <- match.arg(method, c("delong", "bootstrap"))
    # delong NA to pAUC: warn + change
    if (is.numeric(attr(roc$auc, "partial.auc")) && length(attr(roc$auc, "partial.auc") == 2) && method == "delong") {
      warning("Using DeLong for partial AUC is not supported. Using bootstrap instead.")
      method <- "bootstrap"
    }
    else if (class(roc) == "smooth.roc") {
      warning("Using DeLong's test for smoothed ROCs is not supported. Using bootstrap instead.")
      method <- "bootstrap"
    }
  }
  
  if (method == "delong") {
    n <- length(roc$controls)
    m <- length(roc$cases)  
    V <- roc.utils.delong.placements(roc)
    var <- var(V$Y) / n + var(V$X) / m
  }
  else {
    var <- var.roc.bootstrap(roc, boot.n, boot.stratified, progress, ...)
  }
  
  if (percent) {
    var <- var * (100^2)
  }
  return(var)
}

var.roc.bootstrap <- function(roc, boot.n, boot.stratified, progress, ...) {

  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="AUC variance", label="Bootstrap in progress...", ...)

  ## Smoothed ROC curve variance
  if (class(roc) == "smooth.roc") {
    smoothing.args <- roc$smoothing.args
    smoothing.args$smooth <- TRUE
    non.smoothed.roc <- attr(roc, "roc")
    non.smoothed.roc$percent <- FALSE # as we did earlier for the smoothed.roc
    smooth.roc.call <- as.call(c(match.fun("smooth.roc"), roc$smoothing.args))
    auc.args <- attributes(roc$auc)[grep("partial.auc", names(attributes(roc$auc)))]
    auc.call <- as.call(c(match.fun("auc.smooth.roc"), auc.args))

    if (boot.stratified) {
      aucs <- unlist(rlply(boot.n, stratified.ci.smooth.auc(non.smoothed.roc, smooth.roc.call, auc.call), .progress=progress))
    }
    else {
      aucs <- unlist(rlply(boot.n, nonstratified.ci.smooth.auc(non.smoothed.roc, smooth.roc.call, auc.call), .progress=progress))
    }
  }
  ## Non smoothed ROC curves variance
  else {
    if (boot.stratified) {
      aucs <- unlist(rlply(boot.n, stratified.ci.auc(roc), .progress=progress)) # ci.auc: returns aucs just as we need for var, so re-use it!
    }
    else {
      aucs <- unlist(rlply(boot.n, nonstratified.ci.auc(roc), .progress=progress))
    }
  }

  if ((num.NAs <- sum(is.na(aucs))) > 0) {
    warning(sprintf("%i NA value(s) produced during bootstrap were ignored.", num.NAs))
    aucs <- aucs[!is.na(aucs)]
  }
  return(var(aucs))
}


## SD

sd <- function(...)
  UseMethod("sd")

sd.default <- function(...) {
  stats::sd(...)
}

sd.auc <- function(auc, ...) {
  sqrt(var.roc(attr(auc, "roc"), ...))
}

sd.smooth.roc <- function(smooth.roc, ...) {
  sqrt(var.roc(smooth.roc, ...)) 
}

sd.roc <- function(roc, ...) {
  sqrt(var.roc(roc, ...)) 
}
