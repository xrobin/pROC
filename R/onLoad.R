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

.onLoad <- function(lib, pkg) {
  # Generate progressbar option with smart default values
  if (is.null(getOption("pROCProgress"))) {
    if (interactive()) {
      # Check the presence of tcltk
      tcltk.present <- length(find.package("tcltk", quiet=TRUE)) > 0
      if (Sys.info()[['sysname']] == "Darwin") {
        options("pROCProgress" = list(name = "text", width = NA, char = "=", style = 3)) # Tcltk looks broken in some macs (just hangs forever)
      }
      else if (!is.null(getOption("STERM")) && getOption("STERM") == "iESS")
        options("pROCProgress" = list(name = "text", width = NA, char = "=", style = 1))
      else if (.Platform$OS.type == "windows")
        options("pROCProgress" = list(name = "win", width = 300))
      else if (tcltk.present && Sys.getenv("DISPLAY") != "")
        options("pROCProgress" = list(name = "tk", width = 300))
      else
        options("pROCProgress" = list(name = "text", width = NA, char = "=", style = 3))
    }
    else {
      options("pROCProgress" = list(name = "none"))
    }
  }
}

.onAttach <- function(lib, pkg) {
  packageStartupMessage("Type 'citation(\"pROC\")' for a citation.")
}
