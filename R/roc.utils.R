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

# Helper functions for the ROC curves. These functions should not be called directly as they peform very specific tasks and do nearly no argument validity checks. Not documented in RD and not exported.

# returns a list of sensitivities (se) and specificities (sp) for the given data. Robust algorithm
roc.utils.perfs.all.safe <- function(thresholds, controls, cases, direction) {
  perf.matrix <- sapply(thresholds, roc.utils.perfs, controls=controls, cases=cases, direction=direction)
  #stopifnot(identical(roc.utils.perfs.all.fast(thresholds, controls, cases, direction), list(se=perf.matrix[2,], sp=perf.matrix[1,])))
  return(list(se=perf.matrix[2,], sp=perf.matrix[1,]))
}


roc.utils.perfs.all.test <- function(thresholds, controls, cases, direction) {
	perfs.safe <- roc.utils.perfs.all.safe(thresholds=thresholds, controls=controls, cases=cases, direction=direction)
	perfs.fast <- roc.utils.perfs.all.fast(thresholds=thresholds, controls=controls, cases=cases, direction=direction)
	perfs.C <- rocUtilsPerfsAllC(thresholds=thresholds, controls=controls, cases=cases, direction=direction)
	if (! (identical(perfs.safe, perfs.fast) && identical(perfs.safe, perfs.C))) {
		pROCpackageDescription <- packageDescription("pROC")
		save(thresholds, controls, cases, direction, pROCpackageDescription, file="pROC_bug.RData")
		stop(sprintf("Bug in pROC: algorithms returned different values. Diagnostic data saved in pROC_bug.RData. Please report this bug to the package maintainer %s", packageDescription("pROC")$Maintainer))
	}
	return(perfs.safe)
}


# returns a list of sensitivities (se) and specificities (sp) for the given data. Fast algorithm
roc.utils.perfs.all.fast <- function(thresholds, controls, cases, direction) {
  ncontrols <- length(controls)
  ncases <- length(cases)
  predictor <- c(controls, cases)
  response <- c(rep(0, length(controls)), rep(1, length(cases)))
  decr <- direction=="<"
  predictor.order <- order(predictor, decreasing=decr)
  predictor.sorted <- predictor[predictor.order]
  response.sorted <- response[predictor.order]
  
  tp <- cumsum(response.sorted==1)
  fp <- cumsum(response.sorted==0)
  se <- tp / ncases
  sp <- (ncontrols - fp) / ncontrols
  # filter duplicate thresholds
  dups.pred <- rev(duplicated(rev(predictor.sorted)))
  dups.sesp <- duplicated(matrix(c(se, sp), ncol=2), MARGIN=1)
  dups <- dups.pred | dups.sesp
  # Make sure we have the right length
  if (sum(!dups) != length(thresholds) - 1) {
    stop(sprintf("Bug in pROC: fast algorithm computed an incorrect number of sensitivities and specificities. Please report this bug to the package maintainer %s", packageDescription("pROC")$Maintainer))
  }
  if (direction == "<") {
    se <- rev(c(0, se[!dups]))
    sp <- rev(c(1, sp[!dups]))
  }
  else {
    se <- c(0, se[!dups])
    sp <- c(1, sp[!dups])
  }
  return(list(se=se, sp=sp))
}

# As roc.utils.perfs.all but returns an "old-style" matrix (pre-fast-algo-compatible)
#roc.utils.perfs.all.matrix <- function(...) {
#  perfs <- roc.utils.perfs.all(...)
#  return(matrix(c(perfs$sp, perfs$se), nrow=2, byrow=TRUE))
#}

# returns a vector with two elements, sensitivity and specificity, given the threshold at which to evaluate the performance, the values of controls and cases and the direction of the comparison, a character '>' or '<' as controls CMP cases
# sp <- roc.utils.perfs(...)[1,]
# se <- roc.utils.perfs(...)[2,]
roc.utils.perfs <- function(threshold, controls, cases, direction) {
  if (direction == '>') {
    tp <- sum(cases <= threshold)
    tn <- sum(controls > threshold)
  }
  else if (direction == '<') {
    tp <- sum(cases >= threshold)
    tn <- sum(controls < threshold)
  }
  # return(c(sp, se))
  return(c(sp=tn/length(controls), se=tp/length(cases)))
}

# as roc.utils.perfs, but for densities
roc.utils.perfs.dens <- function(threshold, x, dens.controls, dens.cases, direction) {
  if (direction == '>') {
    tp <- sum(dens.cases[x <= threshold])
    tn <- sum(dens.controls[x > threshold])
  }
  else if (direction == '<') {
    tp <- sum(dens.cases[x >= threshold])
    tn <- sum(dens.controls[x < threshold])
  }
  # return(c(sp, se))
  return(c(sp=tn/sum(dens.controls), se=tp/sum(dens.cases)))
}

# return the thresholds to evaluate in the ROC curve, given the 'predictor' values. Returns all unique values of 'predictor' plus 2 extreme values
roc.utils.thresholds <- function(predictor) {
  thresholds <- sort(unique(predictor))
  return((c(-Inf, thresholds) + c(thresholds, +Inf))/2)
}

# Find all the local maximas of the ROC curve. Returns a logical vector
roc.utils.max.thresholds.idx <- function(thresholds, sp, se) {
  reversed <- FALSE
  if (is.unsorted(sp)) {
    # make sure SP is sorted increasingly, and sort thresholds accordingly
    thresholds <- rev(thresholds)
    sp <- rev(sp)
    se <- rev(se)
    reversed <- TRUE
  }
  # TODO: find whether the duplicate check is still needed.
  # Should have been fixed by passing only c(controls, cases)
  # instead of whole 'predictor' to roc.utils.thresholds in roc.default
  # but are there other potential issues like that?
  dup <- duplicated(data.frame(sp, se))
  thresholds <- thresholds[!dup]
  sp <- sp[!dup]
  se <- se[!dup]
  # Find the local maximas
  if (length(thresholds) == 1) {
    local.maximas <- TRUE # let's consider that if there is only 1 point, we should print it.
  }
  else if (length(thresholds) == 2) {
    local.maximas <- c(se[1] > se[2], sp[2] > sp[1])
  }
  else {
    local.maximas <- se[1] > se[2]
    for (i in 2:(length(thresholds)-1)) {
      if (sp[i] > sp[i-1] & se[i] > se[i+1])
        local.maximas <- c(local.maximas, TRUE)
      else if (sp[i] > sp[i-1] & se[i] == 0)
        local.maximas <- c(local.maximas, TRUE)
      else if (se[i] > se[i-1] & sp[i] == 1)
        local.maximas <- c(local.maximas, TRUE)
      else
        local.maximas <- c(local.maximas, FALSE)
    }
    local.maximas <- c(local.maximas, sp[length(thresholds)] > sp[length(thresholds)-1])
  }
  if (any(dup)) {
    lms <- rep(FALSE, length(dup))
    lms[!dup] <- local.maximas
    local.maximas <- lms
  }
  if (reversed)
    rev(local.maximas)

  # Remove +-Inf at the limits of the curve
  #local.maximas <- local.maximas & is.finite(thresholds)
  # Question: why did I do that? It breaks coords.roc when partial.auc contains only the extreme point

  return(local.maximas)
}

# Define which progress bar to use
roc.utils.get.progress.bar <- function(name = getOption("pROCProgress")$name, title = "Bootstrap", label = "", width = getOption("pROCProgress")$width, char = getOption("pROCProgress")$char, style = getOption("pROCProgress")$style, ...) {
  if (name == "tk") { # load tcltk if possible
    if (!requireNamespace("tcltk")) {
      # If tcltk cannot be loaded fall back to default text progress bar
      name <- "text"
      style <- 3
      char <- "="
      width <- NA
      warning("Package tcltk required with progress='tk' but could not be loaded. Falling back to text progress bar.")
    }
  }
  if (name == "none")
    progress_none()
  else if (name == "text") {
    # Put some default values if user only passed a name
    if (missing(style) && missing(char) && missing(width) && getOption("pROCProgress")$name != "text") {
      style <- 3
      char <- "="
      width <- NA
    }
    progress_text(char=char, style=style, width=width)
  }
  else if (name == "tk" || name == "win")
    match.fun(paste("progress", name, sep = "_"))(title=title, label=label, width=width)
  else # in the special case someone made a progress_other function
    match.fun(paste("progress", name, sep = "_"))(title=title, label=label, width=width, char=char, style=style)
}

# sort roc curve. Make sure specificities are increasing.
sort.roc <- function(roc) {
  if (is.unsorted(roc$specificities)) {
    roc$sensitivities <- rev(roc$sensitivities)
    roc$specificities <- rev(roc$specificities)
    roc$thresholds <- rev(roc$thresholds)
  }
  return(roc)
}

# sort smoothed roc curve. Make sure specificities are increasing.
sort.smooth.roc <- function(roc) {
  if (is.unsorted(roc$specificities)) {
    roc$sensitivities <- rev(roc$sensitivities)
    roc$specificities <- rev(roc$specificities)
  }
  return(roc)
}

# Arguments which can be returned by coords
roc.utils.match.coords.ret.args <- function(x) {

  valid.ret.args <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
  x <- replace(x, x=="t", "threshold")
  x <- replace(x, x=="npe", "1-npv")
  x <- replace(x, x=="ppe", "1-ppv")
  match.arg(x, valid.ret.args, several.ok=TRUE)
}

# Compute the min/max for partial AUC
# ... with an auc
roc.utils.min.partial.auc.auc <- function(auc) {
  roc.utils.min.partial.auc(attr(auc, "partial.auc"), attr(auc, "percent"))
}

roc.utils.max.partial.auc.auc <- function(roc) {
  roc.utils.max.partial.auc(attr(auc, "partial.auc"), attr(auc, "percent"))
}

# ... with partial.auc/percent
roc.utils.min.partial.auc <- function(partial.auc, percent) {
  if (!identical(partial.auc, FALSE)) {
    min <- sum(ifelse(percent, 100, 1)-partial.auc)*abs(diff(partial.auc))/2/ifelse(percent, 100, 1)
  }
  else {
    min <- 0.5 * ifelse(percent, 100, 1)
  }
  return(min)
}

roc.utils.max.partial.auc <- function(partial.auc, percent) {
  if (!identical(partial.auc, FALSE)) {
    max <- abs(diff(partial.auc))
  }
  else {
    max <- 1 * ifelse(percent, 100, 1)
  }
  return(max)
}

# Checks if the 
# Input: roc object
# Output: boolean, true the curve reaches 100%/100%, false otherwise
roc.utils.is.perfect.curve <- function(roc) {
	best.point <- max(roc$sensitivities + roc$specificities) / ifelse(roc$percent, 100, 1)
	return(abs(best.point - 2) < .Machine$double.eps ^ 0.5) # or best.point == 2, with numerical tolerance
}
