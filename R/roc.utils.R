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
		sessionInfo <- sessionInfo()
		save(thresholds, controls, cases, direction, sessionInfo, file="pROC_bug.RData")
		stop(sprintf("pROC: algorithms returned different values. Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", utils::packageDescription("pROC")$BugReports))
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
  dups.sesp <- duplicated(se) & duplicated(sp)
  dups <- dups.pred | dups.sesp
  # Make sure we have the right length
  if (sum(!dups) != length(thresholds) - 1) {
  	sessionInfo <- sessionInfo()
  	save(thresholds, controls, cases, direction, sessionInfo, file="pROC_bug.RData")
    stop(sprintf("pROC: fast algorithm computed an incorrect number of sensitivities and specificities. Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", utils::packageDescription("pROC")$BugReports))
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
roc.utils.thresholds <- function(predictor, direction) {
  unique.candidates <- sort(unique(predictor))
  thresholds1 <- (c(-Inf, unique.candidates) + c(unique.candidates, +Inf))/2
  thresholds2 <- (c(-Inf, unique.candidates)/2 + c(unique.candidates, +Inf)/2)
  thresholds <- ifelse(abs(thresholds1) > 1e100, thresholds2, thresholds1)
  if (any(ties <- thresholds %in% predictor)) {
  	# If we get here, some thresholds are identical to the predictor
  	# This is caused by near numeric ties that caused the mean to equal
  	# one of the candidate
  	# We need to make sure we select the right threshold more carefully
  	if (direction == '>') {
  		# We have:
  		# tp <- sum(cases <= threshold)
  		# tn <- sum(controls > threshold)
  		# We need to make sure the selected threshold
  		# Corresponds to the lowest observation of the predictor
  		# Identify problematic thresholds
  		# rows <- which(ties)
  		for (tie.idx in which(ties)) {
  			if (thresholds[tie.idx] == unique.candidates[tie.idx - 1]) {
  				# We're already good, nothing to do
  			}
  			else if (thresholds[tie.idx] == unique.candidates[tie.idx]) {
  				thresholds[tie.idx] <- unique.candidates[tie.idx - 1]
  			}
  			else {
  				sessionInfo <- sessionInfo()
  				save(predictor, direction, sessionInfo, file="pROC_bug.RData")
  				stop(sprintf("Couldn't fix near ties in thresholds: %s, %s, %s, %s. Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", thresholds[tie.idx], unique.candidates[tie.idx - 1], unique.candidates[tie.idx], direction, utils::packageDescription("pROC")$BugReports))
  			}
  		}
  	}
  	else if (direction == '<') {
  		# We have:
  		# tp <- sum(cases >= threshold)
  		# tn <- sum(controls < threshold)
  		# We need to make sure the selected threshold
  		# Corresponds to the highest observation of the predictor
  		# Identify the problematic thresholds:
  		# rows <- which(apply(o, 1, any))
  		for (tie.idx in which(ties)) {
  			if (thresholds[tie.idx] == unique.candidates[tie.idx - 1]) {
  				# Easy to fix: should be unique.candidates[tie.idx]
  				thresholds[tie.idx] <- unique.candidates[tie.idx]
  			} else if (thresholds[tie.idx] == unique.candidates[tie.idx]) {
  				# We're already good, nothing to do
  			}
  			else {
  				sessionInfo <- sessionInfo()
  				save(predictor, direction, sessionInfo, file="pROC_bug.RData")
  				stop(sprintf("Couldn't fix near ties in thresholds: %s, %s, %s, %s. Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", thresholds[tie.idx], unique.candidates[tie.idx - 1], unique.candidates[tie.idx], direction, utils::packageDescription("pROC")$BugReports))
  			}
  		}
  	}
  }
  return(thresholds)
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

# The list of valid coordinate arguments, without 'thresholds'
roc.utils.valid.coords <- c("specificity", "sensitivity", "accuracy",
	"tn", "tp", "fn", "fp",
	"npv", "ppv", "fdr",
	"fpr", "tpr", "tnr", "fnr", 
	"1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv",
	"precision", "recall",
	"youden", "closest.topleft")

# Arguments which can be returned by coords
# @param threshold: FALSE for smooth.roc where threshold isn't valid
roc.utils.match.coords.ret.args <- function(x, threshold = TRUE) {
	valid.ret.args <- roc.utils.valid.coords
	if (threshold) {
		valid.ret.args <- c("threshold", valid.ret.args)
	}

	if ("all" %in% x) {
		if (length(x) > 1) {
			stop("ret='all' can't be used with other 'ret' options.")
		}
		x <- valid.ret.args
	}
	x <- replace(x, x=="topleft", "closest.topleft")
	x <- replace(x, x=="t", "threshold")
	x <- replace(x, x=="npe", "1-npv")
	x <- replace(x, x=="ppe", "1-ppv")
	return(match.arg(x, valid.ret.args, several.ok=TRUE))
}

# Arguments which can be used as input for coords
# @param threshold: FALSE for smooth.roc where threshold isn't valid
roc.utils.match.coords.input.args <- function(x, threshold = TRUE) {
	valid.args <- roc.utils.valid.coords
	if (threshold) {
		valid.args <- c("threshold", valid.args)
	}
	x <- replace(x, x=="topleft", "closest.topleft")
	x <- replace(x, x=="t", "threshold")
	x <- replace(x, x=="npe", "1-npv")
	x <- replace(x, x=="ppe", "1-ppv")
	matched <- match.arg(x, valid.args, several.ok=FALSE)
	# We only handle monotone coords
	if (! coord.is.monotone[matched]) {
		stop(sprintf("Coordinate '%s' is not monotone and cannot be used as input.", matched))
	}
	return(matched)
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

# Load package namespace 'pkg'. 
# Input: package name
# Returns: TRUE upon success (or if the package was already loaded)
# Stops if the package can't be loaded
load.suggested.package <- function(pkg) {
	if (requireNamespace(pkg)) {
		return(TRUE)
	}
	else if (interactive()) {
		if (getRversion() < "3.5.0") { # utils::askYesNo not available
			message(sprintf("Package %s not available, do you want to install it now?", pkg))
			auto.install <- utils::menu(c("Yes", "No")) == 1
		}
		else {
			auto.install <- utils::askYesNo(sprintf("Package %s not available, do you want to install it now?", pkg))
		}
		if (isTRUE(auto.install)) {
			utils::install.packages(pkg)
			if (requireNamespace(pkg)) {
				return(TRUE)
			}
			else {
				stop(sprintf("Installation of %s failed!", pkg))
			}
		}
	}
	stop(sprintf("Package '%s' not available.", pkg))
}


# Calculate coordinates
# @param roc: the roc curve, used to guess if data is in percent and number of cases and controls.
# @param thr, se, sp
# @param best.weights: see coords 
# @return data.frame
roc.utils.calc.coords <- function(roc, thr, se, sp, best.weights) {
	ncases <- ifelse(methods::is(roc, "smooth.roc"), length(attr(roc, "roc")$cases), length(roc$cases))
	ncontrols <- ifelse(methods::is(roc, "smooth.roc"), length(attr(roc, "roc")$controls), length(roc$controls))
	substr.percent <- ifelse(roc$percent, 100, 1)
	
	tp <- se * ncases / substr.percent
	fn <- ncases - tp
	tn <- sp * ncontrols / substr.percent
	fp <- ncontrols - tn
	npv <- substr.percent * tn / (tn + fn)
	ppv <- substr.percent * tp / (tp + fp)
	#res <- matrix(NA, nrow = length(ret), ncol = length(se))
	#if ("tp" %in% ret) {}
	accuracy <- substr.percent * (tp + tn) / (ncases + ncontrols)
	precision <- ppv
	recall <- tpr <- se
	fpr <- substr.percent - sp
	tnr <- sp
	fnr <- substr.percent * fn / (tp + fn)
	fdr <- substr.percent * fp / (tp + fp)
	youden <- roc.utils.optim.crit(se, sp, substr.percent, best.weights, "youden")
	closest.topleft <- - roc.utils.optim.crit(se, sp, substr.percent, best.weights, "closest.topleft") / substr.percent
	
	return(cbind(
		threshold=thr,
		sensitivity=se, 
		specificity=sp, 
		accuracy=accuracy, 
		tn=tn, 
		tp=tp,
		fn=fn,
		fp=fp,
		npv=npv,
		ppv=ppv,
		tpr=tpr,
		tnr=tnr,
		fpr=fpr,
		fnr=fnr,
		fdr=fdr,
		"1-specificity"=substr.percent-sp,
		"1-sensitivity"=substr.percent-se,
		"1-accuracy"=substr.percent-accuracy,
		"1-npv"=substr.percent-npv,
		"1-ppv"=substr.percent-ppv,
		precision=precision,
		recall=recall,
		youden=youden,
		closest.topleft=closest.topleft
	))
}

# Match arbitrary user-supplied thresholds to the threshold of the ROC curve.
# We need to be careful to assign x to the right thresholds around exact data point
# values. This means this function cannot look at the ROC thresholds themselves
# but must instead use the predictor values to assess the thresholds exactly.
# Returns the indices of the thresholds x.
# @param roc: the roc curve
# @param x: the threshold to determine indices
# @return integer vector of indices along roc$thresholds/roc$se/roc$sp.
roc.utils.thr.idx <- function(roc, x) {
	cut_points <- sort(unique(roc$predictor))
	thr_idx <- rep(NA_integer_, length(x))
	if (roc$direction == "<") {
		cut_points <- c(cut_points, Inf)
		j <- 1
		o <- order(x)
		for (i in seq_along(x)) {
			t <- x[o[i]]
			while (cut_points[j] < t) {
				j <- j + 1
			}
			thr_idx[o[i]] <- j
		}
	}
	else {
		cut_points <- c(rev(cut_points), -Inf)
		j <- 1
		o <- order(x, decreasing = TRUE)
		for (i in seq_along(x)) {
			t <- x[o[i]]
			while (cut_points[j] > t) {
				j <- j + 1
			}
			thr_idx[o[i]] <- j
		}
	}
	return(thr_idx)
}


# Get optimal criteria Youden or Closest Topleft
# @param se, sp: the roc curve
# @param max: the maximum value, 1 or 100, based on percent. Namely ifelse(percent, 100, 1)
# @param weights: see coords(best.weights=)
# @param method: youden or closest.topleft coords(best.method=)
# @return numeric vector along roc$thresholds/roc$se/roc$sp.
roc.utils.optim.crit <- function(se, sp, max, weights, method) {
	if (is.numeric(weights) && length(weights) == 2) {
		r <- (1 - weights[2]) / (weights[1] * weights[2]) # r should be 1 by default
	}
	else {
		stop("'best.weights' must be a numeric vector of length 2")
	}
	if (weights[2] <= 0 || weights[2] >= 1) {
		stop("prevalence ('best.weights[2]') must be in the interval ]0,1[.")
	}
	
	if (method == "youden") {
		optim.crit <- se + r * sp
	}
	else if (method == "closest.topleft" || method == "topleft") {
		optim.crit <- - ((max - se)^2 + r * (max - sp)^2)
	}
	return(optim.crit)
}

coord.is.monotone <- c(
	"threshold"=TRUE,
	"sensitivity"=TRUE,
	"specificity"=TRUE,
	"accuracy"=FALSE,
	"tn"=TRUE, 
	"tp"=TRUE,
	"fn"=TRUE,
	"fp"=TRUE,       
	"npv"=FALSE,
	"ppv"=FALSE,
	"tpr"=TRUE,
	"tnr"=TRUE,
	"fpr"=TRUE,
	"fnr"=TRUE,
	"fdr"=FALSE,
	"1-specificity"=TRUE,
	"1-sensitivity"=TRUE,
	"1-accuracy"=FALSE,
	"1-npv"=FALSE,
	"1-ppv"=FALSE,
	"precision"=FALSE,
	"recall"=TRUE,
	"youden"=FALSE,
	"closest.topleft"=FALSE
)

coord.is.decreasing <- c(
	"threshold"=NA, # Depends on direction
	"sensitivity"=TRUE,
	"specificity"=FALSE,
	"accuracy"=NA,
	"tn"=FALSE, 
	"tp"=TRUE,
	"fn"=FALSE,
	"fp"=TRUE,       
	"npv"=NA,
	"ppv"=NA,
	"tpr"=TRUE,
	"tnr"=FALSE,
	"fpr"=TRUE,
	"fnr"=FALSE,
	"fdr"=NA,
	"1-specificity"=TRUE,
	"1-sensitivity"=FALSE,
	"1-accuracy"=NA,
	"1-npv"=NA,
	"1-ppv"=NA,
	"precision"=NA,
	"recall"=TRUE,
	"youden"=NA,
	"closest.topleft"=NA
)

# Get response and predictor(s) from a formula.
# This function takes care of all the logic to handle
# weights, subset, na.action etc. It handles formulas with 
# and without data. It rejects weights and certain na.actions.
# @param formula
# @param data
# @param data.missing 
# @param call the call from the parent function
# @param ... the ... from the parent function
# @return a list with 3 elements: response (vector), predictor.names (character),
#         predictors (data.frame).
roc.utils.extract.formula <- function(formula, data, data.missing, call, ...) {
	# Get predictors (easy)
	if (data.missing) {
		predictors <- attr(terms(formula), "term.labels")
	}
	else {
		predictors <- attr(terms(formula, data = data), "term.labels")
	}
	
	indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(call), nomatch=0)
	if (indx[1] == 0) {
		stop("A formula argument is required")
	}
	# Keep the standard arguments and run them in model.frame
	temp <- call[c(1,indx)]  
	temp[[1]] <- as.name("model.frame")
	# Only na.pass and na.fail should be used
	if (indx[5] != 0) {
		na.action.value = as.character(call[indx[5]])
		if (! na.action.value %in% c("na.pass", "na.fail")) {
			warning(paste0(sprintf("Value %s of na.action is not supported ", na.action.value),
						   "and will break pairing in roc.test and are.paired. ",
						   "Please use 'na.rm = TRUE' instead."))
		}
	}
	else {
		temp$na.action = "na.pass"
	}
	# Adjust call with data from caller
	if (data.missing) {
		temp$data <- NULL
	}
	
	# Run model.frame in the parent
	m <- eval.parent(temp, n = 2)
	
	if (!is.null(model.weights(m))) stop("weights are not supported")
	
	return(list(response.name = names(m)[1],
				response = model.response(m),
				predictor.names = predictors,
				predictors = m[predictors]))
}
