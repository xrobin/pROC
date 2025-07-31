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

plot.ci.thresholds <- function(x, length = .01 * ifelse(attr(x, "roc")$percent, 100, 1), col = par("fg"), ...) {
  bounds <- cbind(x$sp, x$se)
  apply(bounds, 1, function(x, ...) {
    suppressWarnings(segments(x[2], x[4], x[2], x[6], col = col, ...))
    suppressWarnings(segments(x[2] - length, x[4], x[2] + length, x[4], col = col, ...))
    suppressWarnings(segments(x[2] - length, x[6], x[2] + length, x[6], col = col, ...))
    suppressWarnings(segments(x[1], x[5], x[3], x[5], col = col, ...))
    suppressWarnings(segments(x[1], x[5] + length, x[1], x[5] - length, col = col, ...))
    suppressWarnings(segments(x[3], x[5] + length, x[3], x[5] - length, col = col, ...))
  }, ...)
  invisible(x)
}

plot.ci.sp <- function(x, type = c("bars", "shape"), length = .01 * ifelse(attr(x, "roc")$percent, 100, 1), col = ifelse(type == "bars", par("fg"), "gainsboro"), no.roc = FALSE, ...) {
  type <- match.arg(type)
  if (type == "bars") {
    sapply(1:dim(x)[1], function(n, ...) {
      se <- attr(x, "sensitivities")[n]
      suppressWarnings(segments(x[n, 1], se, x[n, 3], se, col = col, ...))
      suppressWarnings(segments(x[n, 1], se - length, x[n, 1], se + length, col = col, ...))
      suppressWarnings(segments(x[n, 3], se - length, x[n, 3], se + length, col = col, ...))
    }, ...)
  } else {
    if (length(x[, 1]) < 15) {
      warning("Low definition shape.")
    }
    suppressWarnings(polygon(c(1 * ifelse(attr(x, "roc")$percent, 100, 1), x[, 1], 0, rev(x[, 3]), 1 * ifelse(attr(x, "roc")$percent, 100, 1)), c(0, attr(x, "sensitivities"), 1 * ifelse(attr(x, "roc")$percent, 100, 1), rev(attr(x, "sensitivities")), 0), col = col, ...))
    if (!no.roc) {
      plot(attr(x, "roc"), add = TRUE)
    }
  }
  invisible(x)
}


plot.ci.se <- function(x, type = c("bars", "shape"), length = .01 * ifelse(attr(x, "roc")$percent, 100, 1), col = ifelse(type == "bars", par("fg"), "gainsboro"), no.roc = FALSE, ...) {
  type <- match.arg(type)
  if (type == "bars") {
    sapply(1:dim(x)[1], function(n, ...) {
      sp <- attr(x, "specificities")[n]
      suppressWarnings(segments(sp, x[n, 1], sp, x[n, 3], col = col, ...))
      suppressWarnings(segments(sp - length, x[n, 1], sp + length, x[n, 1], col = col, ...))
      suppressWarnings(segments(sp - length, x[n, 3], sp + length, x[n, 3], col = col, ...))
    }, ...)
  } else {
    if (length(x[, 1]) < 15) {
      warning("Low definition shape.")
    }
    suppressWarnings(polygon(c(0, attr(x, "specificities"), 1 * ifelse(attr(x, "roc")$percent, 100, 1), rev(attr(x, "specificities")), 0), c(1 * ifelse(attr(x, "roc")$percent, 100, 1), x[, 1], 0, rev(x[, 3]), 1 * ifelse(attr(x, "roc")$percent, 100, 1)), col = col, ...))
    if (!no.roc) {
      plot(attr(x, "roc"), add = TRUE)
    }
  }
  invisible(x)
}

plot.ci.coords <- function(x, type = c("bars", "shape"), length = NULL, col = ifelse(type == "bars", par("fg"), "gainsboro"), ...) {
  type <- match.arg(type)
  if (length(x) > 1) {
    warning(sprintf("'ci.coords' object contains multiple coordinates, only %s will be plotted", names(x)[1]))
  }
  if (is.null(length)) {
    x_range <- range(attr(x, "x"))
    length <- (x_range[2] - x_range[1]) / length(attr(x, "x")) / 5
  }
  if (type == "bars") {
    x_val <- attr(x, "x")
    suppressWarnings(segments(x_val, x[[1]][, 1], x_val, x[[1]][, 3], col = col, ...))
    suppressWarnings(segments(x_val - length, x[[1]][, 1], x_val + length, x[[1]][, 1], col = col, ...))
    suppressWarnings(segments(x_val - length, x[[1]][, 3], x_val + length, x[[1]][, 3], col = col, ...))
  } else {
    if (length(x[[1]][, 1]) < 15) {
      warning("Low definition shape.")
    }
    suppressWarnings(polygon(c(attr(x, "x"), rev(attr(x, "x"))), c(x[[1]][, 1], rev(x[[1]][, 3])), col = col, ...))
  }
  invisible(x)
}
