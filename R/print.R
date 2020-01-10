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

print.smooth.roc <- function(x, digits=max(3, getOption("digits") - 3), call=TRUE, ...) {
  # do we print the call?
  if (call)
    cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
  # Always print number of patients, controls, thresholds, levels?
  print.dataline(attr(x, "roc")) # take this from original roc

  # Smoothing
  cat("Smoothing: ")
  if (is.null(x$smoothing.args)) {
    cat("density with controls: ", as.character(x$call[match("density.controls", names(x$call))]), "; and cases: ", as.character(x$call[match("density.cases", names(x$call))]), "\n", sep="")
  }
  else if (x$smoothing.args$method == "density")
    cat("density (bandwidth: ", x$smoothing.args$bw, "; adjust: ", ifelse(is.null(x$smoothing.args$adjust), 1, x$smoothing.args$adjust), ")\n", sep="")
  else if (x$smoothing.args$method == "density") {
    cat("fitting ", x$fit.controls$densfun, " distribution for controls:\n", sep="")
    print(x$fit.controls$estimate)
    cat("fitting ", x$fit.cases$densfun, " distribution for cases:\n", sep="")
    print(x$fit.cases$estimate)
  }
  else
  	cat(x$smoothing.args$method, "\n")

  # AUC if exists
  if (!is.null(x$auc)) {
    print(x$auc, digits=digits, ...)
  }
  else
    cat("Area under the curve not computed.\n")

  # CI if exists, print it
  if(!is.null(x$ci)) {
    print(x$ci, digits=digits, ...)
  }

  invisible(x)
}

print.multiclass.roc <- function(x, digits=max(3, getOption("digits") - 3), call=TRUE, ...) {
  # do we print the call?
  if (call)
    cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
  # get predictor name
  if ("predictor" %in% names(x$call))
    predictor.name <- as.character(x$call[match("predictor", names(x$call))])
  else if (!is.null(x$call$formula)) 
    predictor.name <- attr(terms(as.formula(x$call$formula), data=x$data), "term.labels")
  # Get response
  if ("response" %in% names(x$call))
    response.name <- as.character(x$call[match("response", names(x$call))])
  else if (!is.null(x$call$formula)) {
    formula.attrs <- attributes(terms(as.formula(x$call$formula), data=x$data))
    response.name <- rownames(formula.attrs$factors)[formula.attrs$response]
  }
  cat("Data: ", predictor.name, " with ", length(x$levels), " levels of ", response.name, ": ", paste(x$levels, collapse=", "),  ".\n", sep="")

  # AUC if exists
  if (!is.null(x$auc)) {
    print(x$auc, digits=digits, ...)
  }
  else
    cat("Multi-class area under the curve not computed.\n")

  # CI if exists, print it
  if(!is.null(x$ci)) {
    print(x$ci, digits=digits, ...)
  }

  invisible(x)
}

print.mv.multiclass.roc <- function(x, digits=max(3, getOption("digits") - 3), call=TRUE, ...) {
	# do we print the call?
	if (call)
		cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
	# get predictor name
	if ("predictor" %in% names(x$call))
		predictor.name <- as.character(x$call[match("predictor", names(x$call))])
	else if (!is.null(x$call$formula)) 
		predictor.name <- attr(terms(as.formula(x$call$formula), data=x$data), "term.labels")
	# Get response
	if ("response" %in% names(x$call))
		response.name <- as.character(x$call[match("response", names(x$call))])
	else if (!is.null(x$call$formula)) {
		formula.attrs <- attributes(terms(as.formula(x$call$formula), data=x$data))
		response.name <- rownames(formula.attrs$factors)[formula.attrs$response]
	}
	cat("Data: multivariate predictor ", predictor.name, " with ", length(x$levels), " levels of ", response.name, ": ", paste(x$levels, collapse=", "),  ".\n", sep="")
	
	# AUC if exists
	if (!is.null(x$auc)) {
		print(x$auc, digits=digits, ...)
	}
	else
		cat("Multi-class area under the curve not computed.\n")
	
	# CI if exists, print it
	if(!is.null(x$ci)) {
		print(x$ci, digits=digits, ...)
	}
	
	invisible(x)
}

print.roc <- function(x, digits=max(3, getOption("digits") - 3), call=TRUE, ...) {
  # do we print the call?
  if (call)
    cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
  # Always print number of patients, controls, thresholds, levels?
  print.dataline(x)

  # AUC if exists
  if (!is.null(x$auc)) {
    print(x$auc, digits=digits, ...)
  }
  else
    cat("Area under the curve not computed.\n")

  # CI if exists, print it
  if(!is.null(x$ci)) {
    print(x$ci, digits=digits, ...)
  }

  invisible(x)
}

print.auc <- function(x, digits=max(3, getOption("digits") - 3), ...) {
  if (identical(attr(x, "partial.auc"), FALSE))
    cat("Area under the curve: ", signif(x, digits=digits), ifelse(attr(x, "percent"), "%", ""), "\n", sep="")
  else {
    cat(ifelse(identical(attr(x, "partial.auc.correct"), TRUE), "Corrected p", "P"), "artial area under the curve", sep="")
    cat(" (", attr(x, "partial.auc.focus"), " ", attr(x, "partial.auc")[1], ifelse(attr(x, "percent"), "%", ""), "-", attr(x, "partial.auc")[2], ifelse(attr(x, "percent"), "%", ""), ")", sep="")
    cat(": ", signif(x, digits=digits), ifelse(attr(x, "percent"), "%", ""), "\n", sep="")
  }
  invisible(x)
}

print.multiclass.auc <- function(x, digits=max(3, getOption("digits") - 3), ...) {
  if (identical(attr(x, "partial.auc"), FALSE))
    cat("Multi-class area under the curve: ", signif(x, digits=digits), ifelse(attr(x, "percent"), "%", ""), "\n", sep="")
  else {
    cat("Multi-class ", ifelse(identical(attr(x, "partial.auc.correct"), TRUE), "corrected ", ""), "partial area under the curve", sep="")
    cat(" (", attr(x, "partial.auc.focus"), " ", attr(x, "partial.auc")[1], ifelse(attr(x, "percent"), "%", ""), "-", attr(x, "partial.auc")[2], ifelse(attr(x, "percent"), "%", ""), ")", sep="")
    cat(": ", signif(x, digits=digits), ifelse(attr(x, "percent"), "%", ""), "\n", sep="")
  }
  invisible(x)
}

print.mv.multiclass.auc <- print.multiclass.auc

print.ci.auc <- function(x, digits=max(3, getOption("digits") - 3), ...) {
  signif.ci <- signif(x, digits=digits)
  cat(attr(x, "conf.level")*100, "% CI: ", sep="")
  cat(signif.ci[1], ifelse(attr(attr(x, "auc"), "percent"), "%", ""), "-", signif.ci[3], ifelse(attr(attr(x, "auc"), "percent"), "%", ""), sep="")
  if (attr(x, "method") == "delong")
    cat(" (DeLong)\n", sep="")
  else
    cat(" (", attr(x, "boot.n"), " ", ifelse(attr(x, "boot.stratified"), "stratified", "non-stratified"), " bootstrap replicates)\n", sep="")
  invisible(x)
}

print.ci.thresholds <- function(x, digits=max(3, getOption("digits") - 3), ...) {
  cat(attr(x, "conf.level")*100, "% CI", sep="")
  cat(" (", attr(x, "boot.n"), " ", ifelse(attr(x, "boot.stratified"), "stratified", "non-stratified"), " bootstrap replicates):\n", sep="")
  signif.sp <- signif(x$sp, digits=digits)
  signif.se <- signif(x$se, digits=digits)
  print(data.frame(thresholds=attr(x, "thresholds"), sp.low=signif.sp[,1], sp.median=signif.sp[,2], sp.high=signif.sp[,3], se.low=signif.se[,1], se.median=signif.se[,2], se.high=signif.se[,3]), row.names=FALSE)
  invisible(x)
}

print.ci.sp <- function(x, digits=max(3, getOption("digits") - 3), ...) {
  cat(attr(x, "conf.level")*100, "% CI", sep="")
  cat(" (", attr(x, "boot.n"), " ", ifelse(attr(x, "boot.stratified"), "stratified", "non-stratified"), " bootstrap replicates):\n", sep="")
  signif.sp <- signif(x, digits=digits)
  print(data.frame(se=attr(x, "sensitivities"), sp.low=signif.sp[,1], sp.median=signif.sp[,2], sp.high=signif.sp[,3]), row.names=FALSE)
  invisible(x)
}

print.ci.se <- function(x, digits=max(3, getOption("digits") - 3), ...) {
  cat(attr(x, "conf.level")*100, "% CI", sep="")
  cat(" (", attr(x, "boot.n"), " ", ifelse(attr(x, "boot.stratified"), "stratified", "non-stratified"), " bootstrap replicates):\n", sep="")
  signif.se <- signif(x, digits=digits)
  print(data.frame(sp=attr(x, "specificities"), se.low=signif.se[,1], se.median=signif.se[,2], se.high=signif.se[,3]), row.names=FALSE)
  invisible(x)
}

print.ci.coords <- function(x, digits=max(3, getOption("digits") - 3), ...) {
  cat(attr(x, "conf.level")*100, "% CI", sep="")
  cat(" (", attr(x, "boot.n"), " ", ifelse(attr(x, "boot.stratified"), "stratified", "non-stratified"), " bootstrap replicates):\n", sep="")

  table <- do.call(cbind, x)
  table <- signif(table, digits = digits)
  table <- cbind(x = attr(x, "x"), as.data.frame(table))
  
  colnames.grid <- expand.grid(c("low", "median", "high"), attr(x, "ret"))
  colnames.vec <- paste(colnames.grid$Var2, colnames.grid$Var1, sep=".")
  colnames(table) <- c(attr(x, "input"), colnames.vec)
  rownames(table) <- attr(x, "x")
  
  print(table, row.names=length(attr(x, "ret")) > 1)
  invisible(x)
}

print.dataline <- function(x) {
  # Case / Controls call
  if ("cases" %in%  names(x$call) && "controls" %in% names(x$call)) {
    cat("Data: ", length(x$controls), " controls ", x$direction, " ", length(x$cases), " cases.\n", sep="")
  }
  else {
    # get predictor name
    if ("predictor" %in% names(x$call))
      predictor.name <- as.character(x$call[match("predictor", names(x$call))])
    else if (!is.null(x$call$formula)) 
      predictor.name <- attr(terms(as.formula(x$call$formula), data=x$data), "term.labels")
    else
      return()
    # Get response
    if ("response" %in% names(x$call))
      response.name <- as.character(x$call[match("response", names(x$call))])
    else if (!is.null(x$call$formula)) {
      formula.attrs <- attributes(terms(as.formula(x$call$formula), data=x$data))
      response.name <- rownames(formula.attrs$factors)[formula.attrs$response]
    }
    else if ("x" %in% names(x$call))
      response.name <- as.character(x$call[match("x", names(x$call))])
    else
      return()
    cat("Data: ", predictor.name, " in ", length(x$controls), " controls (", response.name, " ", x$levels[1], ") ", x$direction, " ", length(x$cases), " cases (", response.name, " ", x$levels[2], ").\n", sep="")
  }
}
