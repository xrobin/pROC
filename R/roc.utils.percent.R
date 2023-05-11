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

# Helper functions to safely convert ROC objects from percent=TRUE to percent=FALSE
# and inversely. These are internal and experimental. They shouldn't be exposed
# to the end user.

# Returns a ROC curve with percent=FALSE
roc_utils_unpercent <- function(x) {
  UseMethod("roc_utils_unpercent")
}

roc_utils_unpercent.roc <- function(x) {
	if (x$percent) {
		if (! is.null(x$auc)) {
			x$auc <- roc_utils_unpercent(x$auc)
		}
		x$sensitivities <- x$sensitivities / 100
		x$specificities <- x$specificities / 100
		x$percent <- FALSE
		if (!is.null(x$call)) {
		  x$call$percent <- FALSE
		}
		if (!is.null(x$ci)) {
			x$ci <- roc_utils_unpercent(x$ci)
		}
	}
	
	return(x)
}

roc_utils_unpercent.auc <- function(x) {
	if (attr(x, "percent")) {
		newx <- x / 100
		attributes(newx) <- attributes(x)
		x <- newx
		attr(x, "percent") <- FALSE
		if (is.numeric(attr(x, "partial.auc"))) {
			attr(x, "partial.auc") <- attr(x, "partial.auc") / 100
		}
		if (! is.null(attr(x, "roc"))) {
			attr(x, "roc") <- roc_utils_unpercent(attr(x, "roc"))
		}
	}
	return(x)
}

roc_utils_unpercent.ci.auc <- function(x) {
	if (attr(attr(x, "auc"), "percent")) {
		x[] <- x / 100
		attr(x, "auc") <- roc_utils_unpercent(attr(x, "auc"))
	}
	return(x)
}

roc_utils_unpercent.ci.thresholds <- function(x) {
	if (attr(x, "roc")$percent) {
		x$sensitivity[] <- x$sensitivity / 100
		x$specificity[] <- x$specificity / 100
		attr(x, "roc") <- roc_utils_unpercent(attr(x, "roc"))
	}
	return(x)
}

roc_utils_unpercent.ci.sp <- function(x) {
	if (attr(x, "roc")$percent) {
		x[] <- x / 100
		attr(x, "sensitivities") <- attr(x, "sensitivities") / 100
		rownames(x) <- attr(x, "sensitivities")
		attr(x, "roc") <- roc_utils_unpercent(attr(x, "roc"))
	}
	return(x)
}

roc_utils_unpercent.ci.se <- function(x) {
	if (attr(x, "roc")$percent) {
		x[] <- x / 100
		attr(x, "specificities") <- attr(x, "specificities") / 100
		rownames(x) <- attr(x, "specificities")
		attr(x, "roc") <- roc_utils_unpercent(attr(x, "roc"))
	}
	return(x)
}

roc_utils_unpercent.ci.coords <- function(x) {
	stop("Cannot convert ci.coords object to percent = FALSE")
}

# Returns a ROC curve with percent=TRUE
roc_utils_topercent <- function(x) {
  UseMethod("roc_utils_topercent")
}

roc_utils_topercent.roc <- function(x) {
	if (! x$percent) {
		if (! is.null(x$auc)) {
			x$auc <- roc_utils_topercent(x$auc)
		}
		x$sensitivities <- x$sensitivities * 100
		x$specificities <- x$specificities * 100
		x$percent <- TRUE
		if (!is.null(x$call)) {
		  x$call$percent <- TRUE
		}
		if (!is.null(x$ci)) {
			x$ci <- roc_utils_topercent(x$ci)
		}
	}
 
  return(x)
}

roc_utils_topercent.auc <- function(x) {
	if (! attr(x, "percent")) {
		newx <- x * 100
		attributes(newx) <- attributes(x)
		x <- newx
		attr(x, "percent") <- TRUE
		if (is.numeric(attr(x, "partial.auc"))) {
			attr(x, "partial.auc") <- attr(x, "partial.auc") * 100
		}
		if (! is.null(attr(x, "roc"))) {
			attr(x, "roc") <- roc_utils_topercent(attr(x, "roc"))
		}
	}
	return(x)
}

roc_utils_topercent.ci.auc <- function(x) {
	if (! attr(attr(x, "auc"), "percent")) {
		x[] <- x * 100
		attr(x, "auc") <- roc_utils_topercent(attr(x, "auc"))
	}
	return(x)
}

roc_utils_topercent.ci.thresholds <- function(x) {
	if (! attr(x, "roc")$percent) {
		x$sensitivity[] <- x$sensitivity * 100
		x$specificity[] <- x$specificity * 100
		attr(x, "roc") <- roc_utils_topercent(attr(x, "roc"))
	}
	return(x)
}

roc_utils_topercent.ci.sp <- function(x) {
	if (! attr(x, "roc")$percent) {
		x[] <- x * 100
		attr(x, "sensitivities") <- attr(x, "sensitivities") * 100
		rownames(x) <- paste(attr(x, "sensitivities"), "%", sep="")
		attr(x, "roc") <- roc_utils_topercent(attr(x, "roc"))
	}
	return(x)
}

roc_utils_topercent.ci.se <- function(x) {
	if (! attr(x, "roc")$percent) {
		x[] <- x * 100
		attr(x, "specificities") <- attr(x, "specificities") * 100
		rownames(x) <- paste(attr(x, "specificities"), "%", sep="")
		attr(x, "roc") <- roc_utils_topercent(attr(x, "roc"))
	}
	return(x)
}

roc_utils_topercent.ci.coords <- function(x) {
	stop("Cannot convert ci.coords object to percent = TRUE")
}