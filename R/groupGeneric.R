# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison.
# Copyright (C) 2014 Xavier Robin
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

# A single function for all classes: R only dispatches to the method when
# both operands have identical methods (e.g. an auc and a ci.auc)
Ops.auc <- Ops.ci.se <- Ops.ci.sp <- Ops.ci.auc <- function(e1, e2) {
  e1 <- remove.auc.ci.attributes(e1)
  # e2 is missing for unary operators
  if (!missing(e2)) {
    e2 <- remove.auc.ci.attributes(e2)
  }
  NextMethod()
}

Math.auc <- function(x, ...) {
  attributes(x) <- NULL
  NextMethod()
}


Math.ci.se <- Math.ci.sp <- Math.ci.auc <- function(x, ...) {
  x <- remove.ci.attributes(x)
  NextMethod()
}

remove.auc.ci.attributes <- function(x) {
  if (methods::is(x, "auc")) {
    attributes(x) <- NULL
  } else if (methods::is(x, "ci.se") || methods::is(x, "ci.sp") || methods::is(x, "ci.auc")) {
    x <- remove.ci.attributes(x)
  }
  return(x)
}

remove.ci.attributes <- function(ci) {
  attr(ci, "conf.level") <- NULL
  attr(ci, "boot.n") <- NULL
  attr(ci, "boot.stratified") <- NULL
  attr(ci, "specificities") <- NULL
  attr(ci, "sensitivities") <- NULL
  attr(ci, "roc") <- NULL
  attr(ci, "method") <- NULL
  attr(ci, "auc") <- NULL
  class(ci) <- class(ci)[-(1:2)]
  return(ci)
}
