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

.register_ggplot_add_ggroc_layer <- function(...) {
  registerS3method(
    "ggplot_add",
    "ggroc_layer",
    ggplot_add.ggroc_layer,
    envir = asNamespace("ggplot2")
  )
}

.onLoad <- function(lib, pkg) {
  # ggplot2 is in Suggests; register ggplot_add only when that namespace exists.
  if (isNamespaceLoaded("ggplot2")) {
    .register_ggplot_add_ggroc_layer()
  }
  setHook(packageEvent("ggplot2", "onLoad"), .register_ggplot_add_ggroc_layer)
}

.onAttach <- function(lib, pkg) {
  packageStartupMessage("Type 'citation(\"pROC\")' for a citation.")
}
