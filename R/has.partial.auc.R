# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2011-2014 Xavier Robin, Alexandre Hainard, Natacha Turck,
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

has.partial.auc <- function(roc) {
  UseMethod("has.partial.auc")
}

has.partial.auc.auc <- function(roc) {
  if (is.null(roc)) {
    return(NULL)
  }

  is.numeric(attr(roc, "partial.auc")) && length(attr(roc, "partial.auc") == 2)
}

has.partial.auc.smooth.roc <- function(roc) {
  return(has.partial.auc.roc(roc))
}

has.partial.auc.roc <- function(roc) {
  return(has.partial.auc.auc(roc$auc))
}
