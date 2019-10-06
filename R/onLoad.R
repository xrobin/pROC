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
      if (!is.null(getOption("STERM")) && getOption("STERM") == "iESS")
        options("pROCProgress" = list(name = "text", width = NA, char = "=", style = 1))
      else if (.Platform$OS.type == "windows")
        options("pROCProgress" = list(name = "win", width = 300))
      else
        options("pROCProgress" = list(name = "text", width = NA, char = "=", style = 3))
    }
    else {
      options("pROCProgress" = list(name = "none"))
    }
  }
}

.parseRcppVersion <- function(rcpp.version) {
  # Parses Rcpp version integer into a string.
  # Eg "65538" -> "1.0.2"
  rcpp.version <- as.integer(rcpp.version)
  major <- rcpp.version %/% 65536
  rcpp.version <- rcpp.version - major * 65536
  minor <- rcpp.version %/% 256
  rcpp.version <- rcpp.version - minor * 256
  rev <- rcpp.version
  return(sprintf("%s.%s.%s", major, minor, rev))
}

.checkRcppVersion <- function() {
  # Check runtime version of Rcpp is the same than we had at compile time
  runtime_version <- package_version(utils::packageVersion("Rcpp"))
  build_version <- package_version(.parseRcppVersion(RcppVersion()))
  if (runtime_version != build_version) {
    warning(sprintf("It seems pROC was compiled with Rcpp version %s, but %s is available now. Please re-install pROC to avoid problems: install.packages(\"pROC\").",
                    build_version,runtime_version))
  }
}

.onAttach <- function(lib, pkg) {
  packageStartupMessage("Type 'citation(\"pROC\")' for a citation.")
}
