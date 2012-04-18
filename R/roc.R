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

roc <- function(...) {
  UseMethod("roc")
}

roc.formula <- function (formula, data, ...){
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  term.labels <- attr(attr(m, "terms"), "term.labels")
  response <- model.extract(m, "response")

  if (length(term.labels) != 1) {
    stop("Invalid formula: exactly 1 predictor is required in a formula of type response~predictor.")
  }
  if (length(response) == 0) {
    stop("Error in the formula: a response is required in a formula of type response~predictor.")
  }

  roc <- roc.default(response, m[[term.labels]], ...)
  roc$call <- match.call()
  return(roc)
}

roc.default <- function(response, predictor,
                        controls, cases,
                        density.controls, density.cases,
                        # data interpretation
                        levels=base::levels(as.factor(response)), # precise the levels of the responses as c("control group", "positive group"). Can be used to ignore some response levels.
                        
                        percent=FALSE, # Must sensitivities, specificities and AUC be reported in percent? Note that if TRUE, and you want a partial area, you must pass it in percent also (partial.area=c(100, 80))
                        na.rm=TRUE,
                        direction=c("auto", "<", ">"), # direction of the comparison. Auto: automatically define in which group the median is higher and take the good direction to have an AUC >= 0.5

                        # what computation must be done
                        smooth=FALSE, # call smooth.roc on the current object
                        auc=TRUE, # call auc.roc on the current object
                        ci=FALSE, # call ci.roc on the current object
                        plot=FALSE, # call plot.roc on the current object

                        # disambiguate method for ci and smooth
                        smooth.method="binormal",
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
    original.predictor <- predictor # store a copy of the original predictor (before converting ordered to numeric)
    # ensure predictor is numeric
    if (!is.numeric(predictor)) {
      if (is.ordered(predictor))
        predictor <- as.numeric(predictor)
      else
        stop("Predictor must be numeric or ordered.")
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
    patients.in.levels <- response == levels[1] | response == levels[2]
    response <- response[patients.in.levels]
    predictor <- predictor[patients.in.levels]
  }

  # Cases / Controls
  else if (!missing(cases) && !is.null(cases) && !missing(controls) && !is.null(controls)) {
    # check data consistency
    if (length(controls) == 0)
      stop("No control observation.")
    if (length(cases) == 0)
      stop("No case observation.")
    if (!is.numeric(cases))
      stop("Cases must be numeric or ordered.")
    if (!is.numeric(controls))
      stop("Controls must be numeric or ordered.")
    # build response/predictor
    response <- c(rep(0, length(controls)), rep(1, length(cases)))
    predictor <- c(controls, cases)
    original.predictor <- c(controls, cases)
    # remove nas
    if (na.rm) {
      if (any(is.na(controls)))
        controls <- na.omit(controls)
      if (any(is.na(cases)))
        controls <- na.omit(cases)
    }
    else if (any(is.na(c(controls, cases)))) # Unable to compute anything if there is any NA in the data we want to consider !
      return(NA)
    # are there empty cats?
    if (length(controls) == 0)
      stop("No control observation.")
    if (length(cases) == 0)
      stop("No case observation.")
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
      if (direction == "auto" && smooth.roc$auc < 0.5) {
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

  # create the roc object
  roc <- list()
  class(roc) <- "roc"
  roc$levels <- levels
  roc$percent <- percent
  roc$call <- match.call()

  if (direction == "auto" && median(controls) <= median(cases))
    direction <- "<"
  else if (direction == "auto" && median(controls) > median(cases))
    direction <- ">"

  # compute SE / SP
  thresholds <- roc.utils.thresholds(c(controls, cases))
  perfs <- roc.utils.perfs.all(thresholds=thresholds, predictor=predictor, response=response, ncontrols=length(controls), ncases=length(cases), direction=direction, levels=levels)
  se <- perfs$se
  sp <- perfs$sp

  if (length(thresholds) != length(se)) {
    stop("New SE/SP computation is wrong (inconsistent number of se/sp)")
  }

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
  roc$original.predictor <- original.predictor
  roc$predictor <- predictor
  roc$response <- response

  # smooth?
  if (smooth) {
    if (missing(density.controls))
      density.controls <- density
    if (missing(density.cases))
      density.cases <- density
    roc <- smooth.roc(roc, density=density, density.controls=density.controls, density.cases=density.cases, method=smooth.method, ...)
    roc$call <- match.call()
  }

  # compute AUC
  if (auc)
    roc$auc <- auc(roc, ...)
  # compute CI
  if (ci)
    roc$ci <- ci(roc, method=ci.method, ...)
  # plot
  if (plot)
    plot.roc(roc, ...)

  # return roc
  return(roc)
}
