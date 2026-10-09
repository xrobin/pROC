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

plot.roc <- function(x, ...) {
  UseMethod("plot.roc")
}

plot.roc.formula <- function(x, data, subset, na.action, ...) {
  data.missing <- missing(data)
  call <- match.call()
  names(call)[2] <- "formula" # forced to be x by definition of plot
  roc.data <- roc_utils_extract_formula(
    formula = x, data, subset, na.action, ...,
    data.missing = data.missing,
    call = call
  )
  if (length(roc.data$predictor.name) > 1) {
    stop("Only one predictor supported in 'plot.roc'.")
  }
  response <- roc.data$response
  predictor <- roc.data$predictors[, 1]

  roc <- roc(response, predictor, plot = TRUE, ...)
  if (inherits(roc, c("roc", "smooth.roc"))) {
    roc$call <- match.call()
  }
  invisible(roc)
}

plot.roc.default <- function(x, predictor, ...) {
  roc <- roc(x, predictor, plot = TRUE, ...)
  if (inherits(roc, c("roc", "smooth.roc"))) {
    roc$call <- match.call()
  }
  invisible(roc)
}

plot.roc.smooth.roc <- plot.smooth.roc <- function(x, ...) {
  invisible(plot.roc.roc(x, ...)) # force usage of plot.roc.roc: only print.thres not working
}

plot.roc.roc <- function(x,
                         add = FALSE,
                         reuse.auc = TRUE,
                         axes = TRUE,
                         legacy.axes = FALSE,
                         xlim = if (x$percent) {
                           c(100, 0)
                         } else {
                           c(1, 0)
                         },
                         ylim = if (x$percent) {
                           c(0, 100)
                         } else {
                           c(0, 1)
                         },
                         xlab = ifelse(x$percent, ifelse(legacy.axes, "100 - Specificity (%)", "Specificity (%)"), ifelse(legacy.axes, "1 - Specificity", "Specificity")),
                         ylab = ifelse(x$percent, "Sensitivity (%)", "Sensitivity"),
                         asp = 1,
                         mar = c(4, 4, 2, 2) + .1,
                         mgp = c(2.5, 1, 0),
                         # col, lty and lwd for the ROC line only
                         col = par("col"),
                         lty = par("lty"),
                         lwd = 2,
                         type = "l",
                         # Identity line
                         identity = !add,
                         identity.col = "darkgrey",
                         identity.lty = 1,
                         identity.lwd = 1,
                         # Print the thresholds on the plot
                         print.thres = FALSE,
                         print.thres.pch = 20,
                         print.thres.adj = c(-.05, 1.25),
                         print.thres.col = "black",
                         print.thres.pattern = NULL,
                         print.thres.cex = par("cex"),
                         print.thres.pattern.cex = print.thres.cex,
                         print.thres.best.method = NULL,
                         print.thres.best.weights = c(1, 0.5),
                         # Print the AUC on the plot
                         print.auc = FALSE,
                         print.auc.pattern = NULL,
                         print.auc.x = ifelse(x$percent, 50, .5),
                         print.auc.y = ifelse(x$percent, 50, .5),
                         print.auc.adj = c(0, 1),
                         print.auc.col = col,
                         print.auc.cex = par("cex"),
                         # Grid
                         grid = FALSE,
                         grid.v = {
                           if (is.logical(grid) && grid[1] == TRUE) {
                             seq(0, 1, 0.1) * ifelse(x$percent, 100, 1)
                           } else if (is.numeric(grid)) {
                             seq(0, ifelse(x$percent, 100, 1), grid[1])
                           } else {
                             NULL
                           }
                         },
                         grid.h = {
                           if (length(grid) == 1) {
                             grid.v
                           } else if (is.logical(grid) && grid[2] == TRUE) {
                             seq(0, 1, 0.1) * ifelse(x$percent, 100, 1)
                           } else if (is.numeric(grid)) {
                             seq(0, ifelse(x$percent, 100, 1), grid[2])
                           } else {
                             NULL
                           }
                         },
                         # for grid.lty, grid.lwd and grid.col, a length 2 value specifies both values for vertical (1) and horizontal (2) grid
                         grid.lty = 3,
                         grid.lwd = 1,
                         grid.col = "#DDDDDD",
                         # Polygon for the auc
                         auc.polygon = FALSE,
                         auc.polygon.col = "gainsboro", # Other arguments can be passed to polygon() using "..." (for these two we cannot)
                         auc.polygon.lty = par("lty"),
                         auc.polygon.density = NULL,
                         auc.polygon.angle = 45,
                         auc.polygon.border = NULL,
                         # Should we show the maximum possible area as another polygon?
                         max.auc.polygon = FALSE,
                         max.auc.polygon.col = "#EEEEEE", # Other arguments can be passed to polygon() using "..." (for these two we cannot)
                         max.auc.polygon.lty = par("lty"),
                         max.auc.polygon.density = NULL,
                         max.auc.polygon.angle = 45,
                         max.auc.polygon.border = NULL,
                         # Confidence interval
                         ci = !is.null(x$ci),
                         ci.type = c("bars", "shape", "no"),
                         ci.col = ifelse(ci.type == "bars", par("fg"), "gainsboro"),
                         # Hooks to draw add-ons underneath (panel.first) or
                         # on top of (panel.last) everything plot.roc draws
                         panel.first = NULL,
                         panel.last = NULL,
                         ...) {
  percent <- x$percent

  if (max.auc.polygon | auc.polygon | print.auc) { # we need the auc here
    if (is.null(x$auc) | !reuse.auc) {
      x$auc <- auc(x, ...)
    }
  }

  # get and sort the sensitivities and specificities
  se <- sort(x$sensitivities, decreasing = TRUE)
  sp <- sort(x$specificities, decreasing = FALSE)
  if (!add) {
    opar <- par(mar = mar, mgp = mgp)
    on.exit(par(opar))
    # type="n" to plot background lines and polygon shapes first. We will add the line later. axes=FALSE, we'll add them later according to legacy.axis
    suppressWarnings(plot(x$specificities, x$sensitivities, xlab = xlab, ylab = ylab, type = "n", axes = FALSE, xlim = xlim, ylim = ylim, lwd = lwd, asp = asp, ...))

    # As we had axes=FALSE we need to add them again unless axes=FALSE
    if (axes) {
      box()
      # axis behave differently when at and labels are passed (no decimals on 1 and 0),
      # so handle each case separately and consistently across axes
      if (legacy.axes) {
        lab.at <- axTicks(side = 1)
        lab.labels <- format(ifelse(x$percent, 100, 1) - lab.at)
        suppressWarnings(axis(side = 1, at = lab.at, labels = lab.labels, ...))
        lab.at <- axTicks(side = 2)
        suppressWarnings(axis(side = 2, at = lab.at, labels = format(lab.at), ...))
      } else {
        suppressWarnings(axis(side = 1, ...))
        suppressWarnings(axis(side = 2, ...))
      }
    }
  }

  # Plot the grid
  # make sure grid.lty, grid.lwd and grid.col are at least of length 2
  grid.lty <- rep(grid.lty, length.out = 2)
  grid.lwd <- rep(grid.lwd, length.out = 2)
  grid.col <- rep(grid.col, length.out = 2)
  if (!is.null(grid.v)) {
    suppressWarnings(abline(v = grid.v, lty = grid.lty[1], col = grid.col[1], lwd = grid.lwd[1], ...))
  }
  if (!is.null(grid.h)) {
    suppressWarnings(abline(h = grid.h, lty = grid.lty[2], col = grid.col[2], lwd = grid.lwd[2], ...))
  }

  # panel.first: drawn here, after the grid and before the max-AUC polygon,
  # so add-on polygons (polygon_auc, polygon_max_auc) land exactly where
  # auc.polygon/max.auc.polygon draw them below. Only forced on a new plot:
  # on add=TRUE it would redraw a previous curve's background a second time.
  if (!add) {
    panel.first
  }

  # Plot the polygon displaying the maximal area
  if (max.auc.polygon) {
    roc_utils_draw_max_auc_polygon(x$auc,
      col = max.auc.polygon.col, lty = max.auc.polygon.lty,
      density = max.auc.polygon.density, angle = max.auc.polygon.angle,
      border = max.auc.polygon.border, ...
    )
  }
  # Plot the ci shape
  if (ci && !methods::is(x$ci, "ci.auc")) {
    ci.type <- match.arg(ci.type)
    if (ci.type == "shape") {
      plot(x$ci, type = "shape", col = ci.col, no.roc = TRUE, ...)
    }
  }
  # Plot the polygon displaying the actual area
  if (auc.polygon) {
    roc_utils_draw_auc_polygon(x, x$auc,
      col = auc.polygon.col, lty = auc.polygon.lty,
      density = auc.polygon.density, angle = auc.polygon.angle,
      border = auc.polygon.border, ...
    )
  }
  # Identity line
  if (identity) suppressWarnings(abline(ifelse(percent, 100, 1), -1, col = identity.col, lwd = identity.lwd, lty = identity.lty, ...))
  # Actually plot the ROC curve
  suppressWarnings(lines(sp, se, type = type, lwd = lwd, col = col, lty = lty, ...))
  # Plot the ci bars
  if (ci && !methods::is(x$ci, "ci.auc")) {
    if (ci.type == "bars") {
      plot(x$ci, type = "bars", col = ci.col, ...)
    }
  }
  if (is.null(print.thres.pattern)) {
    print.thres.pattern <- if (!methods::is(x, "smooth.roc") && roc_utils_is_ordered_roc(x)) {
      ifelse(x$percent, "%s (%.1f%%, %.1f%%)", "%s (%.3f, %.3f)")
    } else {
      ifelse(x$percent, "%.1f (%.1f%%, %.1f%%)", "%.3f (%.3f, %.3f)")
    }
  }

  # Print the thresholds on the curve if print.thres is TRUE
  if (isTRUE(print.thres)) {
    print.thres <- "best"
  }
  special.thres <- coords_special_x(print.thres, roc = x, keywords = c("no", "all", "local maximas", "best"))
  if (!is.na(special.thres)) {
    print.thres <- special.thres
  }
  if (identical(print.thres, "no") || identical(print.thres, FALSE)) {
    # do nothing
  } else if (methods::is(x, "smooth.roc")) {
    if (is.numeric(print.thres)) {
      stop("Numeric 'print.thres' unsupported on a smoothed ROC plot.")
    } else if (print.thres == "all" || print.thres == "local maximas") {
      stop("'all' and 'local maximas' 'print.thres' unsupported on a smoothed ROC plot.")
    } else if (print.thres == "best") {
      co <- coords(x, print.thres, best.method = print.thres.best.method, best.weights = print.thres.best.weights, transpose = FALSE)
      suppressWarnings(points(co$specificity, co$sensitivity, pch = print.thres.pch, cex = print.thres.cex, col = print.thres.col, ...))
      suppressWarnings(text(co$specificity, co$sensitivity, sprintf(print.thres.pattern, NA, co$specificity, co$sensitivity), adj = print.thres.adj, cex = print.thres.pattern.cex, col = print.thres.col, ...))
    }
  } else if (is.numeric(print.thres) || is.character(print.thres) || is.ordered(print.thres)) {
    co <- coords(x, print.thres, best.method = print.thres.best.method, best.weights = print.thres.best.weights, transpose = FALSE)
    suppressWarnings(points(co$specificity, co$sensitivity, pch = print.thres.pch, cex = print.thres.cex, col = print.thres.col, ...))
    thres.lab <- if (is.ordered(co$threshold) || is.character(co$threshold)) as.character(co$threshold) else co$threshold
    suppressWarnings(text(co$specificity, co$sensitivity, sprintf(print.thres.pattern, thres.lab, co$specificity, co$sensitivity), adj = print.thres.adj, cex = print.thres.pattern.cex, col = print.thres.col, ...))
  }

  # Print the AUC on the plot
  if (print.auc) {
    ci.auc.obj <- if (ci && methods::is(x$ci, "ci.auc")) x$ci else NULL
    roc_utils_draw_auc_text(x$auc,
      ci = ci.auc.obj, pattern = print.auc.pattern,
      xy = c(print.auc.x, print.auc.y), adj = print.auc.adj,
      col = print.auc.col, cex = print.auc.cex, ...
    )
  }

  # panel.last: drawn on top of everything plot.roc draws, for symmetry with
  # plot.default. Unlike panel.first, forced on every call (add=TRUE too):
  # each call's own panel.last is independent per-curve annotation, not tied
  # to the one underlying plot frame.
  panel.last

  invisible(x)
}
