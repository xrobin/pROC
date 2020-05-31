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

roc <- function(...) {
  UseMethod("roc")
}

roc.formula <- function (formula, data, ...) {
	data.missing <- missing(data)
	roc.data <- roc.utils.extract.formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = match.call())
	response <- roc.data$response
	predictors <- roc.data$predictors
	
	if (length(response) == 0) {
		stop("Error in the formula: a response is required in a formula of type response~predictor.")
	}
	
	if (ncol(predictors) == 1) {
		roc <- roc.default(response, predictors[, 1], ...)
		roc$call <- match.call()
		if (!is.null(roc$smooth))
			attr(roc, "roc")$call <- roc$call
		return(roc)
	}
	else if (ncol(predictors) > 1) {
		roclist <- lapply(roc.data$predictor.names, function(predictor, formula, m.data, call, ...) {
			# Get one ROC
			roc <- roc.default(response, m.data[[predictor]], ...)
			# Update the call to reflect the parents
			formula[3] <- call(predictor) # replace the predictor in formula
			call$formula <- formula # Replace modified formula
			roc$call <- call
			return(roc)
		}, formula = formula, m.data = predictors, call = match.call(), ...)
		# Set the list names
		names(roclist) <- roc.data$predictor.names
		return(roclist)
	}
	else {
		stop("Invalid formula:at least 1 predictor is required in a formula of type response~predictor.")
	}
}

roc.data.frame <- function(data, response, predictor, 
                           ret = c("roc", "coords", "all_coords"),
                           ...) {
  ret <- match.arg(ret)
  
  if (is.character(substitute(response))) {
  	response_name <- response
  }
  else {
  	response_name <- deparse(substitute(response))
  }
  
  if (is.character(substitute(predictor))) {
  	predictor_name <- predictor
  }
  else {
  	predictor_name <- deparse(substitute(predictor))
  }
  
  if (any(! c(response_name, predictor_name) %in% colnames(data))) {
  	# Some column is not in data. This could be a genuine error or the user not aware or NSE and wants to use roc_ instead
  	warning("This method uses non-standard evaluation (NSE). Did you want to use the `roc_` function instead?")
  }
  
  r <- roc_(data, response_name, predictor_name, ret = ret, ...)
  
  if (ret == "roc") {
    r$call <- match.call()
  }
  return(r)
}

roc_ <- function(data, response, predictor,
                 ret = c("roc", "coords", "all_coords"),
                 ...) {
  ret <- match.arg(ret)
  
  # Ensure the data contains the columns we need
  # In case of an error we want to show the name of the data. If the function
  # was called from roc.data.frame we want to deparse in that environment instead
  if (sys.nframe() > 1 && deparse(sys.calls()[[sys.nframe()-1]][[1]]) == "roc.data.frame") {
  	data_name <- deparse(substitute(data, parent.frame(n = 1)))
  }
  else {
  	data_name <- deparse(substitute(data))
  }
  if (! response %in% colnames(data)) {
  	stop(sprintf("Column %s not present in data %s", response, data_name))
  }
  if (! predictor %in% colnames(data)) {
  	stop(sprintf("Column '%s' not present in data %s", predictor, data_name))
  }
  
  r <- roc(data[[response]], data[[predictor]], ...)
  
  if (ret == "roc") {
    r$call <- match.call()
    return(r)
  }
  else if (ret == "coords") {
    co <- coords(r, x = "all", transpose = FALSE)
    rownames(co) <- NULL
    return(co)
  }
  else if (ret == "all_coords") {
    co <- coords(r, x = "all", ret="all", transpose = FALSE)
    rownames(co) <- NULL
    return(co)
  }
}

roc.default <- function(response, predictor,
                        controls, cases,
                        density.controls, density.cases,
                        # data interpretation
                        levels=base::levels(as.factor(response)), # precise the levels of the responses as c("control group", "positive group"). Can be used to ignore some response levels.
                        
                        percent=FALSE, # Must sensitivities, specificities and AUC be reported in percent? Note that if TRUE, and you want a partial area, you must pass it in percent also (partial.area=c(100, 80))
                        na.rm=TRUE,
                        direction=c("auto", "<", ">"), # direction of the comparison. Auto: automatically define in which group the median is higher and take the good direction to have an AUC >= 0.5
                        algorithm=6,
						quiet = FALSE,

                        # what computation must be done
                        smooth=FALSE, # call smooth.roc on the current object
                        auc=TRUE, # call auc.roc on the current object
                        ci=FALSE, # call ci.roc on the current object
                        plot=FALSE, # call plot.roc on the current object

                        # disambiguate method/n for ci and smooth
                        smooth.method="binormal",
						smooth.n=512,
                        ci.method=NULL,
                        # capture density for smooth.roc here (do not pass to graphical functions)
                        density=NULL,
                        
                        # further arguments passed to plot, auc, ci, smooth.
                        ... 
                       
                        ) {
  # Check arguments
  direction <- match.arg(direction)
  # Response / Predictor
  if (!missing(response) && !is.null(response) && !missing(predictor) && !is.null(predictor)) {
  	# Forbid case/controls
  	if ((!missing(cases) && !is.null(cases)) || (!missing(controls) && !is.null(controls))) {
  		stop("'cases/controls' argument incompatible with 'response/predictor'.")
  	}
  	# Forbid density
  	if ((!missing(density.cases) && !is.null(density.cases)) || (!missing(density.controls) && !is.null(density.controls))) {
  		stop("'density.*' arguments incompatible with 'response/predictor'.")
  	}
  	
    original.predictor <- predictor # store a copy of the original predictor (before converting ordered to numeric and removing NA)
    original.response <- response # store a copy of the original predictor (before converting ordered to numeric)
    
    # Validate levels
    if (missing(levels)) {
    	if (length(levels) > 2) {
    		warning("'response' has more than two levels. Consider setting 'levels' explicitly or using 'multiclass.roc' instead")
    		levels <- levels[1:2]
    	}
    	else if (length(levels) < 2) {
    		stop("'response' must have two levels")
    	}
    	ifelse(quiet, invisible, message)(sprintf("Setting levels: control = %s, case = %s", levels[1], levels[2]))
    }
    else if (length(levels) != 2) {
    	stop("'levels' argument must have length 2")
    }

    # ensure predictor is numeric or ordered
    if (!is.numeric(predictor)) {
      if (is.ordered(predictor)) {
      	predictor <- tryCatch(
      		{
	      		as.numeric(as.character(predictor))
	      	},
	      	warning = function(warn) {
	      		warning("Ordered predictor converted to numeric vector. Threshold values will not correspond to values in predictor.")
	      		return(as.numeric(predictor))
	      	}
      	)
      }
      else {
        stop("Predictor must be numeric or ordered.")
      }
    }
    if (is.matrix(predictor)) {
    	warning("Deprecated use a matrix as predictor. Unexpected results may be produced, please pass a numeric vector.")
    }
    if (is.matrix(response)) {
    	warning("Deprecated use a matrix as response. Unexpected results may be produced, please pass a vector or factor.")
    }
    # also make sure response and predictor are vectors of the same length
    if (length(predictor) != length(response)) {
      stop("Response and predictor must be vectors of the same length.")
    }
    # remove NAs if requested
    if (na.rm) {
      nas <- is.na(response) | is.na(predictor)
      if (any(nas)) {
        na.action <- grep(TRUE, nas)
        class(na.action) <- "omit"
        response <- response[!nas]
        attr(response, "na.action") <- na.action
        predictor <- predictor[!nas]
        attr(predictor, "na.action") <- na.action
      }
    }
    else if(any(is.na(c(predictor[response==levels[1]], predictor[response==levels[2]], response)))) # Unable to compute anything if there is any NA in the response or in the predictor data we want to consider !
      return(NA)
    splitted <- split(predictor, response)
    controls <- splitted[[as.character(levels[1])]]
    if (length(controls) == 0)
      stop("No control observation.")
    cases <- splitted[[as.character(levels[2])]]
    if (length(cases) == 0)
      stop("No case observation.")

    # Remove patients not in levels
    patients.in.levels <- response %in% levels
    if (!all(patients.in.levels)) {
      response <- response[patients.in.levels]
      predictor <- predictor[patients.in.levels]
    }
    
    # Check infinities
    if (any(which <- is.infinite(predictor))) {
    	warning("Infinite values(s) in predictor, cannot build a valid ROC curve. NaN returned instead.")
    	return(NaN)
    }
  }

  # Cases / Controls
  else if (!missing(cases) && !is.null(cases) && !missing(controls) && !is.null(controls)) {
  	# Forbid density
  	if ((!missing(density.cases) && !is.null(density.cases)) || (!missing(density.controls) && !is.null(density.controls))) {
  		stop("'density.*' arguments incompatible with 'response/predictor'.")
  	}
  	# remove nas
  	if (na.rm) {
  		if (any(is.na(controls)))
  			controls <- na.omit(controls)
  		if (any(is.na(cases)))
  			cases <- na.omit(cases)
  	}
  	else if (any(is.na(c(controls, cases)))) # Unable to compute anything if there is any NA in the data we want to consider !
  		return(NA)
  	# are there empty cats?
  	if (length(controls) == 0)
  		stop("No control observation.")
  	if (length(cases) == 0)
  		stop("No case observation.")
  	
    # check data consistency
    if (is.ordered(cases)) {
    	if (is.ordered(controls)) {
    		if (identical(attr(cases, "levels"), attr(controls, "levels"))) {
    			# merge
    			original.predictor <- ordered(c(as.character(cases), as.character(controls)), levels = attr(controls, "levels"))
    			# Predictor, control and cases must be numeric from now on
    			predictor <- as.numeric(original.predictor)
    			controls <- as.numeric(controls)
    			cases <- as.numeric(cases)
    		}
    		else {
    			stop("Levels of cases and controls differ.")
    		}
    	}
    	else {
    		stop("Cases are of ordered type but controls are not.")
    	}
    }
    else if (is.numeric(cases)) {
    	if (is.numeric(controls)) {
    		# build response/predictor
    		predictor <- c(controls, cases)
    		original.predictor <- predictor
    	}
    	else {
    		stop("Cases are of numeric type but controls are not.")
    	}
    }
    else { 
    	stop("Cases and controls must be numeric or ordered.")
    }
  	
  	# Check infinities
  	if (any(which <- is.infinite(predictor))) {
  		warning("Infinite values(s) in predictor, cannot build a valid ROC curve. NaN returned instead.")
  		return(NaN)
  	}
    
    response <- c(rep(0, length(controls)), rep(1, length(cases)))
    original.response <- response
    levels <- c(0, 1)
  }
  else if (!missing(density.cases) && !is.null(density.cases) && !missing(density.controls) && !is.null(density.controls)) {
    if (!is.numeric(density.cases) || !is.numeric(density.controls))
      stop("'density.cases' and 'density.controls' must be numeric values of density (over the y axis).")
    if (direction == "auto")
      dir <- "<"
    else
      dir <- direction
    smooth.roc <- smooth.roc.density(density.controls=density.controls, density.cases=density.cases, percent=percent, direction=dir)
    class(smooth.roc) <- "smooth.roc"
    smooth.roc <- sort(smooth.roc) # sort se and sp
    # anchor SE/SP at 0/100
    smooth.roc$specificities <- c(0, as.vector(smooth.roc$specificities), ifelse(percent, 100, 1))
    smooth.roc$sensitivities <- c(ifelse(percent, 100, 1), as.vector(smooth.roc$sensitivities), 0)
    smooth.roc$percent <- percent # keep some basic roc specifications
    smooth.roc$direction <- direction
    smooth.roc$call <- match.call()
    if (auc) {
      smooth.roc$auc <- auc(smooth.roc, ...)
      if (direction == "auto" && smooth.roc$auc < roc.utils.min.partial.auc.auc(smooth.roc$auc)) {
        smooth.roc <- roc.default(density.controls=density.controls, density.cases=density.cases, levels=levels,
                                  percent=percent, direction=">", auc=auc, ci=ci, plot=plot, ...)
        smooth.roc$call <- match.call()
        return(smooth.roc)
      }
    }
    if (ci)
      warning("CI can not be computed with densities.")
    if (plot)
      plot.roc(smooth.roc, ...)
    return(smooth.roc)
  }
  else {
    stop("No valid data provided.")
  }

  if (direction == "auto" && median(controls) <= median(cases)) {
  	direction <- "<"
  	ifelse(quiet, invisible, message)("Setting direction: controls < cases")
  }
  else if (direction == "auto" && median(controls) > median(cases)) {
  	direction <- ">"
  	ifelse(quiet, invisible, message)("Setting direction: controls > cases")
  }
  
  # smooth with densities, but only density was provided, not density.controls/cases
  if (smooth) {
    if (missing(density.controls))
      density.controls <- density
    if (missing(density.cases))
      density.cases <- density
  }
  
  # Choose algorithm
  if (isTRUE(algorithm == 6)) {
    if (is.numeric(predictor)) {
      algorithm <- 2
    }
    else {
      algorithm <- 3
    }
  }
  else if (isTRUE(algorithm == 0)) {
  	load.suggested.package("microbenchmark")
    cat("Starting benchmark of algorithms 2 and 3, 10 iterations...\n")
    thresholds <- roc.utils.thresholds(c(controls, cases), direction)
    benchmark <- microbenchmark::microbenchmark(
#      "1" = roc.utils.perfs.all.safe(thresholds=thresholds, controls=controls, cases=cases, direction=direction),
      "2" = roc.utils.perfs.all.fast(thresholds=thresholds, controls=controls, cases=cases, direction=direction),
      "3" = rocUtilsPerfsAllC(thresholds=thresholds, controls=controls, cases=cases, direction=direction),
      times = 10
    )
    print(summary(benchmark))
    if (any(is.na(benchmark))) {
      warning("Microbenchmark returned NA. Using default algorithm 1.")
      algorithm <- 2
    }
    algorithm <- as.integer(names(which.min(tapply(benchmark$time, benchmark$expr, sum))))
    cat(sprintf("Selecting algorithm %s.\n", algorithm))
  }
  else if (isTRUE(algorithm == 5)) {
    thresholds <- length(roc.utils.thresholds(c(controls, cases), direction))
    if (thresholds > 55) { # critical number determined in inst/extra/algorithms.speed.test.R
      algorithm <- 2
    } else {
      algorithm <- 3
    }
  }
  
  if (isTRUE(algorithm == 2)) {
    fun.sesp <- roc.utils.perfs.all.fast
  }
  else if (isTRUE(algorithm  == 3)) {
    fun.sesp <- rocUtilsPerfsAllC
  }
  else if (isTRUE(algorithm ==  1)) {
    fun.sesp <- roc.utils.perfs.all.safe
  }
  else if (isTRUE(algorithm == 4)) {
    fun.sesp <- roc.utils.perfs.all.test
  }
  else {
    stop("Unknown algorithm (must be 0, 1, 2, 3, 4 or 5).")
  }

  roc <- roc.cc.nochecks(controls, cases,
             percent=percent,
             direction=direction,
             fun.sesp=fun.sesp,
             smooth = smooth, density.cases = density.cases,  density.controls = density.controls, smooth.method = smooth.method, smooth.n = smooth.n,
             auc, ...)
  
  roc$call <- match.call()
  if (smooth) {
    attr(roc, "roc")$call <- roc$call
    attr(roc, "roc")$original.predictor <- original.predictor
    attr(roc, "roc")$original.response <- original.response
    attr(roc, "roc")$predictor <- predictor
    attr(roc, "roc")$response <- response
    attr(roc, "roc")$levels <- levels
  }
  roc$original.predictor <- original.predictor
  roc$original.response <- original.response
  roc$predictor <- predictor
  roc$response <- response
  roc$levels <- levels
  
  if (auc) {
  	attr(roc$auc, "roc") <- roc
  }
  
  # compute CI
  if (ci)
    roc$ci <- ci(roc, method=ci.method, ...)
  # plot
  if (plot)
    plot.roc(roc, ...)
  
  # return roc
  return(roc)
}

#' Creates a ROC object from response, predictor, ... without argument checking. Not to be exposed to the end user
roc.rp.nochecks <- function(response, predictor, levels, ...) {
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(levels[1])]]
  if (length(controls) == 0)
    stop("No control observation.")
  cases <- splitted[[as.character(levels[2])]]
  if (length(cases) == 0)
    stop("No case observation.")
  roc.cc.nochecks(controls, cases, ...)
}

#' Creates a ROC object from controls, cases, ... without argument checking. Not to be exposed to the end user
roc.cc.nochecks <- function(controls, cases, percent, direction, fun.sesp, smooth, smooth.method, smooth.n, auc, ...) {
  # create the roc object
  roc <- list()
  class(roc) <- "roc"
  roc$percent <- percent

  # compute SE / SP
  thresholds <- roc.utils.thresholds(c(controls, cases), direction)
  perfs <- fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=direction)

  se <- perfs$se
  sp <- perfs$sp

  if (percent) {
    se <- se*100
    sp <- sp*100
  }

  # store the computations in the roc object
  roc$sensitivities <- se
  roc$specificities <- sp
  roc$thresholds <- thresholds
  roc <- sort(roc)
  roc$direction <- direction
  roc$cases <- cases
  roc$controls <- controls
  roc$fun.sesp <- fun.sesp
  
  if (smooth) {
    roc <- smooth.roc(roc, method=smooth.method, n=smooth.n, ...)
  }
  # compute AUC
  if (auc)
    roc$auc <- auc.roc(roc, ...)
  
  return(roc)
}
