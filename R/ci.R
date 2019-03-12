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
	# Get the data. Use standard code from survival::coxph as suggested by Terry Therneau
	Call <- match.call()
	indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(Call), nomatch=0)
	if (indx[1] == 0) {
		stop("A formula argument is required")
	}
	# Keep the standard arguments and run them in model.frame
	temp <- Call[c(1,indx)]  
	temp[[1]] <- as.name('model.frame')
	m <- eval(temp, parent.frame())
	
	if (!is.null(model.weights(m))) stop("weights are not supported")
	
	response <- model.response(m)
	predictor <- m[[attr(terms(formula), "term.labels")]]
	ci.roc(roc(response, predictor, ...), ...)
}

ci.default <- function(response, predictor, ...) {
  ci.roc(roc.default(response, predictor, ...), ...)
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

