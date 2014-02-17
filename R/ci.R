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
  ci.roc(roc.formula(formula, data, ...), ...)
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

  return(ci)
}
