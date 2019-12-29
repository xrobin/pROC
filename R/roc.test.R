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

roc.test <- function(...) {
	UseMethod("roc.test")
}

roc.test.formula <- function (formula, data, ...) {
	data.missing <- missing(data)
	call <- match.call()
	roc.data <- roc.utils.extract.formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = call)
	if (length(roc.data$predictor.name) != 2) {
		stop("Invalid formula: exactly 2 predictors are required in a formula of type response~predictor1+predictor2.")
	}
	response <- roc.data$response
	predictors <- roc.data$predictors
	
	testres <- roc.test.default(response, predictors, ...)
	testres$call <- call
	# data.names for pretty print()ing
	if (data.missing) {
		testres$data.names <- sprintf("%s and %s by %s (%s, %s)", roc.data$predictor.names[1], roc.data$predictor.names[2], roc.data$response.name, testres$roc1$levels[1], testres$roc1$levels[2])
	}
	else {
		testres$data.names <- sprintf("%s and %s in %s by %s (%s, %s)", roc.data$predictor.names[1], roc.data$predictor.names[2], deparse(substitute(data)), roc.data$response.name, testres$roc1$levels[1], testres$roc1$levels[2])
	}
	
	return(testres)
}

roc.test.default <- function(response, predictor1, predictor2=NULL, na.rm=TRUE, method=NULL, ...) {
	if (is.matrix(predictor1) | is.data.frame(predictor1)) {
		if (!is.null(predictor2))
			stop("Predictor2 must not be specified if predictor1 is a matrix or a data.frame.")
		if (dim(predictor1)[2] == 2 & length(response) == dim(predictor1)[1]) {
			roc1 <- roc(response, predictor1[,1], ...)
			roc2 <- roc(response, predictor1[,2], ...)
			if (!is.null(names(predictor1)))
				data.names <- sprintf("%s and %s in %s by %s (%s, %s)", names(predictor1)[1], names(predictor1)[2], deparse(substitute(predictor1)), deparse(substitute(response)), roc1$levels[1], roc1$levels[2])
			else if (!is.null(colnames(predictor1)))
				data.names <- sprintf("%s and %s in %s by %s (%s, %s)", colnames(predictor1)[1], colnames(predictor1)[2], deparse(substitute(predictor1)), deparse(substitute(response)), roc1$levels[1], roc1$levels[2])
			else
				data.names <- sprintf("%s by %s (%s, %s)", deparse(substitute(predictor1)), deparse(substitute(response)), roc1$levels[1], roc1$levels[2])
		}
		else {
			stop("Wrong dimension for predictor1 as a matrix or a data.frame.")
		}
	}
	else {
		if (missing(predictor2))
			stop("Missing argument predictor2 with predictor1 as a vector.")
		# Need to remove NAs
		if (na.rm) {
			nas <- is.na(response) | is.na(predictor1) | is.na(predictor2)
			response <- response[!nas]
			predictor1 <- predictor1[!nas]
			predictor2 <- predictor2[!nas]
		}
		roc1 <- roc(response, predictor1, ...)
		roc2 <- roc(response, predictor2, ...)
		call <- match.call()
		data.names <- sprintf("%s and %s by %s (%s, %s)", deparse(call$predictor1), deparse(call$predictor2), deparse(call$response), roc1$levels[1], roc1$levels[2])
	}
	test <- roc.test.roc(roc1, roc2, method=method, ...)
	test$data.names <- data.names
	return(test)
}

roc.test.auc <- function(roc1, roc2, ...) {
	# First save the names
	data.names <- paste(deparse(substitute(roc1)), "and", deparse(substitute(roc2)))
	# Change roc1 from an auc to a roc object but keep the auc specifications
	auc1 <- roc1
	attr(auc1, "roc") <- NULL
	roc1 <- attr(roc1, "roc")
	roc1$auc <- auc1
	# Pass to roc.test.roc
	testres <- roc.test.roc(roc1, roc2, ...)
	testres$call <- match.call()
	testres$data.names <- data.names
	return(testres)
}

roc.test.smooth.roc <- function(roc1, roc2, ...) {
	testres <- roc.test.roc(roc1, roc2, ...)
	testres$call <- match.call()
	testres$data.names <- paste(deparse(substitute(roc1)), "and", deparse(substitute(roc2)))
	return(testres)
}

roc.test.roc <- function(roc1, roc2,
						 method=c("delong", "bootstrap", "venkatraman", "sensitivity", "specificity"),
						 sensitivity=NULL, specificity=NULL,
						 alternative = c("two.sided", "less", "greater"),
						 paired=NULL,
						 reuse.auc=TRUE,
						 boot.n=2000, boot.stratified=TRUE,
						 ties.method="first",
						 progress=getOption("pROCProgress")$name,
						 parallel=FALSE,
						 ...) {
	alternative <- match.arg(alternative)
	data.names <- paste(deparse(substitute(roc1)), "and", deparse(substitute(roc2)))
	# If roc2 is an auc, take the roc but keep the auc specifications
	if (methods::is(roc2, "auc")) {
		auc2 <- roc2
		attr(auc2, "roc") <- NULL
		roc2 <- attr(roc2, "roc")
		roc2$auc <- auc2
	}
	
	if (roc.utils.is.perfect.curve(roc1) && roc.utils.is.perfect.curve(roc2)) {
		warning("roc.test() of two ROC curves with AUC == 1 has always p.value = 1 and can be misleading.")
	}
	
	# store which objects are smoothed, and how
	smoothing.args <- list()
	if (methods::is(roc1, "smooth.roc")) {
		smoothing.args$roc1 <- roc1$smoothing.args
		smoothing.args$roc1$smooth <- TRUE
		roc1 <- attr(roc1, "roc")
	}
	else {
		smoothing.args$roc1 <- list(smooth=FALSE)
	}
	if (methods::is(roc2, "smooth.roc")) {
		smoothing.args$roc2 <- roc2$smoothing.args
		smoothing.args$roc2$smooth <- TRUE
		roc2 <- attr(roc2, "roc")
	}
	else {
		smoothing.args$roc2 <- list(smooth=FALSE)
	}
	
	# Check if we do a paired or unpaired roc.test
	if (is.null(paired)) {
		# then determine whether the rocs are paired or not
		rocs.are.paired <- are.paired(roc1, roc2, return.paired.rocs=TRUE, reuse.auc=TRUE, reuse.ci=FALSE, reuse.smooth=TRUE)
		if (rocs.are.paired) {
			paired <- TRUE
			roc1 <- attr(rocs.are.paired, "roc1")
			roc2 <- attr(rocs.are.paired, "roc2")
		}
		else {
			paired <- FALSE
			roc1 <- roc1
			roc2 <- roc2
		}
	}
	else if (paired) {
		# make sure the rocs are really paired
		rocs.are.paired <- rocs.are.paired <- are.paired(roc1, roc2, return.paired.rocs=TRUE, reuse.auc=TRUE, reuse.ci=FALSE, reuse.smooth=TRUE)
		if (! rocs.are.paired) 
			stop("The paired ROC test cannot be applied to unpaired curves.")
		roc1 <- attr(rocs.are.paired, "roc1")
		roc2 <- attr(rocs.are.paired, "roc2")
	}
	else { # assume unpaired
		rocs.are.paired <- are.paired(roc1, roc2, return.paired.rocs=FALSE)
		if (rocs.are.paired) 
			warning("The ROC curves seem to be paired. Consider performing a paired roc.test.")
		roc1 <- roc1
		roc2 <- roc2
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
				stop("DeLong's test is not supported for partial AUC. Use method=\"bootstrap\" instead.")
			}
			if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth) {
				stop("DeLong's test is not supported for smoothed ROCs. Use method=\"bootstrap\" instead.")
			}
			if (roc1$direction != roc2$direction)
				warning("DeLong's test should not be applied to ROC curves with a different direction.")
		}
		else if (method == "venkatraman") {
			if (has.partial.auc(roc1))
				stop("Partial AUC is not supported for Venkatraman's test.")
			if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth)
				stop("Venkatraman's test is not supported for smoothed ROCs")
			if (roc1$direction != roc2$direction)
				warning("Venkatraman's test should not be applied to ROC curves with different directions.")
			if (alternative != "two.sided") {
				stop("Only two-sided tests are available for Venkatraman.")
			}
		}
	}
	
	# Prepare the return value htest
	if (smoothing.args$roc1$smooth)
		estimate <- do.call("smooth.roc", c(list(roc=roc1), smoothing.args$roc1))$auc
	else
		estimate <- roc1$auc
	if (smoothing.args$roc2$smooth)
		estimate <- c(estimate, do.call("smooth.roc", c(list(roc=roc2), smoothing.args$roc2))$auc)
	else
		estimate <- c(estimate, roc2$auc)
	if (identical(attr(roc1$auc, "partial.auc"), FALSE)) {
		nest <- paste(ifelse(smoothing.args$roc1$smooth, "Smoothed ", ""), "AUC of roc1", sep="")
	}
	else {
		nest <- paste(ifelse (attr(roc1$auc, "partial.auc.correct"), "Corrected ", ""),
					  ifelse (smoothing.args$roc1$smooth, "Smoothed ", ""),
					  "pAUC (", attr(roc1$auc, "partial.auc")[1], "-", attr(roc1$auc, "partial.auc")[2], " ", attr(roc1$auc, "partial.auc.focus"),
					  ") of roc1", sep="")
	}
	if (identical(attr(roc2$auc, "partial.auc"), FALSE)) {
		nest <- c(nest, paste(ifelse(smoothing.args$roc2$smooth, "Smoothed ", ""), "AUC of roc2", sep=""))
	}
	else {
		nest <- c(nest, paste(ifelse (attr(roc2$auc, "partial.auc.correct"), "Corrected ", ""),
							  ifelse (smoothing.args$roc2$smooth, "Smoothed ", ""),
							  "pAUC (", attr(roc2$auc, "partial.auc")[1], "-", attr(roc2$auc, "partial.auc")[2], " ", attr(roc2$auc, "partial.auc.focus"),
							  ") of roc2", sep=""))
	}
	nest <- sub("Corrected Smoothed", "Corrected smoothed", nest) # no upper on smoothed if corrected.
	names(estimate) <- nest
	null.value <- 0
	names(null.value) <- "difference in AUC"
	htest <- list(
		alternative = alternative,
		data.names = data.names,
		estimate = estimate,
		null.value = null.value
	)
	class(htest) <- "htest"
	
	if (method == "delong") {
		if (paired) {
			stat <- delong.paired.test(roc1, roc2)
			names(stat) <- "Z"
			htest$statistic <- stat
			htest$method <- "DeLong's test for two correlated ROC curves"
			
			if (alternative == "two.sided")
				pval <- 2*pnorm(-abs(stat))
			else if (alternative == "greater")
				pval <- pnorm(-stat)
			else
				pval <- pnorm(stat)
			htest$p.value <- pval
		}
		else {
			stats <- delong.unpaired.test(roc1, roc2)
			stat <- stats[1]
			df <- stats[2]
			htest$statistic <- c("D"=stat)
			htest$parameter <- c("df"=df)
			htest$method <- "DeLong's test for two ROC curves"
			
			if (alternative == "two.sided")
				pval <- 2*pt(-abs(stat), df=df)
			else if (alternative == "greater")
				pval <- pt(-stat, df=df)
			else
				pval <- pt(stat, df=df)
			htest$p.value <- pval
		}
	}
	else if (method == "venkatraman") {
		if(class(progress) != "list")
			progress <- roc.utils.get.progress.bar(progress, title="Venkatraman ROC test", label="Permutations in progress...", ...)
		if (paired) {
			stats <- venkatraman.paired.test(roc1, roc2, boot.n, ties.method, progress, parallel)
			htest$method <- "Venkatraman's test for two paired ROC curves"
		}
		else {
			stats <- venkatraman.unpaired.test(roc1, roc2, boot.n, ties.method, progress, parallel)
			htest$method <- "Venkatraman's test for two unpaired ROC curves"
		}
		stat <- stats[[1]]
		names(stat) <- "E"
		htest$statistic <- stat
		parameter <- c(boot.n)
		names(parameter) <- "boot.n"
		htest$parameter <- parameter
		pval <- sum(stats[[2]]>=stats[[1]])/boot.n
		htest$p.value <- pval
		names(null.value) <- "difference in ROC operating points"
		htest$estimate <- NULL # AUC not relevant in venkatraman
	}
	else { # method == "bootstrap" or "sensitivity" or "specificity"
		# Check if called with density.cases or density.controls
		if (is.null(smoothing.args) || is.numeric(smoothing.args$density.cases) || is.numeric(smoothing.args$density.controls))
			stop("Cannot compute the statistic on ROC curves smoothed with numeric density.controls and density.cases.")
		
		if(class(progress) != "list")
			progress <- roc.utils.get.progress.bar(progress, title="Bootstrap ROC test", label="Bootstrap in progress...", ...)
		
		if (method == "specificity") {
			if (! is.numeric(specificity) || length(specificity) != 1) {
				stop("Argument 'specificity' must be numeric of length 1 for a specificity test.")
			}
			stat <- bootstrap.test(roc1, roc2, "sp", specificity, paired, boot.n, boot.stratified, smoothing.args, progress, parallel)
			if (paired)
				htest$method <- "Specificity test for two correlated ROC curves"
			else
				htest$method <- "Specificity test for two ROC curves"
		}
		else if (method == "sensitivity") {
			if (! is.numeric(sensitivity) || length(sensitivity) != 1) {
				stop("Argument 'sensitivity' must be numeric of length 1 for a sensitivity test.")
			}
			stat <- bootstrap.test(roc1, roc2, "se", sensitivity, paired, boot.n, boot.stratified, smoothing.args, progress, parallel)
			if (paired)
				htest$method <- "Sensitivity test for two correlated ROC curves"
			else
				htest$method <- "Sensitivity test for two ROC curves"
		}
		else {
			stat <- bootstrap.test(roc1, roc2, "boot", NULL, paired, boot.n, boot.stratified, smoothing.args, progress, parallel)
			if (paired)
				htest$method <- "Bootstrap test for two correlated ROC curves"
			else
				htest$method <- "Bootstrap test for two ROC curves"
		}
		stat <- as.vector(stat) # remove auc attributes
		names(stat) <- "D"
		htest$statistic <- stat
		parameter <- c(boot.n, boot.stratified)
		names(parameter) <- c("boot.n", "boot.stratified")
		htest$parameter <- parameter
		
		if (alternative == "two.sided")
			pval <- 2*pnorm(-abs(stat))
		else if (alternative == "greater")
			pval <- pnorm(-stat)
		else
			pval <- pnorm(stat)
		htest$p.value <- pval
	}
	
	htest$roc1 <- roc1
	htest$roc2 <- roc2
	# Remove name from p value
	htest$p.value <- unname(htest$p.value)
	# Restore smoothing if necessary
	if (smoothing.args$roc1$smooth)
		htest$roc1 <- do.call("smooth.roc", c(list(roc=roc1), smoothing.args$roc1))
	if (smoothing.args$roc2$smooth)
		htest$roc2 <- do.call("smooth.roc", c(list(roc=roc2), smoothing.args$roc2))
	return(htest)
}
