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

ci <- function(...) {
  UseMethod("ci")
}

ci.formula <- function(formula, data, ...) {
	data.missing <- missing(data)
	roc.data <- roc.utils.extract.formula(formula, data, ..., 
										  data.missing = data.missing,
										  call = match.call())
	if (length(roc.data$predictor.name) > 1) {
		stop("Only one predictor supported in 'ci'.")
	}
	response <- roc.data$response
	predictor <- roc.data$predictors[, 1]
	ci.roc(roc(response, predictor, ...), ...)
}

ci.default <- function(response, predictor, ...) {
	roc <- roc.default(response, predictor, ci = FALSE, ...)
	if (methods::is(roc, "smooth.roc")) {
		return(ci.roc(smooth.roc = roc, ...))
	}
	else {
		return(ci.roc(roc = roc, ...))
	}
}

ci.smooth.roc <- function(smooth.roc, of = c("auc", "sp", "se", "coords"), ...) {
  of <- match.arg(of)
  
  if (of == "auc")
    ci <- ci.auc.smooth.roc(smooth.roc, ...)
  else if (of == "sp")
    ci <- ci.sp.smooth.roc(smooth.roc, ...)
  else if (of == "se")
    ci <- ci.se.smooth.roc(smooth.roc, ...)
  else if (of == "coords")
  	ci <- ci.coords.smooth.roc(smooth.roc, ...)
  else
  	stop(sprintf("Unknown 'of' for CI: %s", of))

  return(ci)
}

ci.roc <- function(roc, of = c("auc", "thresholds", "sp", "se", "coords"), ...) {
  of <- match.arg(of)
  
  if (of == "auc")
    ci <- ci.auc.roc(roc, ...)
  else if (of == "thresholds")
    ci <- ci.thresholds.roc(roc, ...)
  else if (of == "sp")
    ci <- ci.sp.roc(roc, ...)
  else if (of == "se")
    ci <- ci.se.roc(roc, ...)
  else if (of == "coords")
  	ci <- ci.coords.roc(roc, ...)
  else
  	stop(sprintf("Unknown 'of' for CI: %s", of))

  return(ci)
}

ci.multiclass.roc <- function(multiclass.roc, of = "auc", ...) {
	stop("CI of a multiclass ROC curve not implemented")
}

ci.multiclass.auc <- function(multiclass.auc, of = "auc", ...) {
	stop("CI of a multiclass AUC not implemented")
}

