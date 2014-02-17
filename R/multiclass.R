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

multiclass.roc <- function(...)
  UseMethod("multiclass.roc")

multiclass.roc.formula <- function(formula, data, ...) {
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  term.labels <- attr(attr(m, "terms"), "term.labels")
  response <- model.extract(m, "response")

  if (length(response) == 0) {
    stop("Error in the formula: a response is required in a formula of type response~predictor.")
  }

  multiclass.roc <- multiclass.roc.default(response, m[[term.labels]], ...)
  multiclass.roc$call <- match.call()
  return(multiclass.roc)
}

multiclass.roc.default <- function(response, predictor,
                                   levels=base::levels(as.factor(response)),
                                   percent=FALSE, # Must sensitivities, specificities and AUC be reported in percent? Note that if TRUE, and you want a partial area, you must pass it in percent also (partial.area=c(100, 80))
                                   # what computation must be done
                                   #auc=TRUE, # call auc.roc on the current object
                                   #ci=FALSE, # call ci.roc on the current object
                                   ...) {
  multiclass.roc <- list(
                         response = response,
                         predictor = predictor,
                         percent = percent,
                         call=match.call())
  class(multiclass.roc) <- "multiclass.roc"

  if ("ordered" %in% class(response) && any(names(table(response))[table(response) == 0] %in% levels)) {
    missing.levels <- names(table(response))[table(response) == 0]
    missing.levels.requested <- missing.levels[missing.levels %in% levels]
    warning(paste("No observation for response level(s):", paste(missing.levels.requested, collapse=", ")))
    levels <- levels[!(levels %in% missing.levels.requested)]
  }
  multiclass.roc$levels <- levels
  
  rocs <- combn(levels, 2, function(X, response, predictor, percent, ...) {
    roc(response, predictor, levels=X, percent=percent, auc=FALSE, ci=FALSE, ...)
  }, simplify=FALSE, response=response, predictor=predictor, percent=percent, ...)

  multiclass.roc$rocs <- rocs

  # Makes no sense to turn auc off, so remove this option
  #if (auc)
    multiclass.roc$auc <- auc.multiclass.roc(multiclass.roc, ...)
  # CI is not implemented yet.
  #if (ci)
  #  multiclass.roc$ci <- ci.multiclass.roc(multiclass.roc, ...)

  return(multiclass.roc)
}
