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

lines.roc <- function(x, ...) {
  UseMethod("lines.roc")
}

lines.roc.formula <- function(x, data, subset, na.action, ...) {
	data.missing <- missing(data)
	call <- match.call()
	names(call)[2] <- "formula" # forced to be x by definition of lines
	roc.data <- roc.utils.extract.formula(formula=x, data, subset, na.action, ..., 
										  data.missing = data.missing,
										  call = call)
	if (length(roc.data$predictor.name) > 1) {
		stop("Only one predictor supported in 'lines.roc'.")
	}
	response <- roc.data$response
	predictor <- roc.data$predictors[, 1]
	roc <- roc(response, predictor, ...)
	lines.roc.roc(roc, ...)
	roc$call <- match.call()
	invisible(roc)
}

lines.roc.default <- function(x, predictor, ...) {
  roc <- roc(x, predictor, ...)
  lines.roc.roc(roc, ...)
  roc$call <- match.call()
  invisible(roc)
}

lines.roc.smooth.roc <- lines.smooth.roc <- function(x, ...) {
  lines.roc.roc(x, ...) # force usage of lines.roc.roc
}

lines.roc.roc <- function(x, lwd=2, ...) {
  suppressWarnings(lines(x$specificities, x$sensitivities, lwd=lwd, ...))
  invisible(x)
}
