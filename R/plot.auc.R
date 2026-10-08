# Internal drawing helpers shared by plot.roc.roc (R/plot.roc.R) and the
# exported add-on functions below, so both produce identical output by
# construction: plot.roc(auc.polygon=TRUE, max.auc.polygon=TRUE, print.auc=TRUE)
# and polygon_max_auc()/polygon_auc()/text() draw exactly the same shapes.

roc_utils_draw_max_auc_polygon <- function(auc, col = "#EEEEEE", lty = par("lty"),
                                            density = NULL, angle = 45, border = NULL, ...) {
  percent <- attr(auc, "percent")
  partial.auc <- attr(auc, "partial.auc")
  partial.auc.focus <- attr(auc, "partial.auc.focus")
  if (identical(partial.auc, FALSE)) {
    map.y <- c(0, 1, 1, 0) * ifelse(percent, 100, 1)
    map.x <- c(1, 1, 0, 0) * ifelse(percent, 100, 1)
  } else if (partial.auc.focus == "sensitivity") {
    map.y <- c(partial.auc[2], partial.auc[2], partial.auc[1], partial.auc[1])
    map.x <- c(0, 1, 1, 0) * ifelse(percent, 100, 1)
  } else {
    map.y <- c(0, 1, 1, 0) * ifelse(percent, 100, 1)
    map.x <- c(partial.auc[2], partial.auc[2], partial.auc[1], partial.auc[1])
  }
  suppressWarnings(polygon(map.x, map.y, col = col, lty = lty, border = border, density = density, angle = angle, ...))
  invisible(NULL)
}

roc_utils_draw_auc_polygon <- function(roc, auc, col = "gainsboro", lty = par("lty"),
                                        density = NULL, angle = 45, border = NULL, ...) {
  percent <- attr(auc, "percent")
  partial.auc <- attr(auc, "partial.auc")
  partial.auc.focus <- attr(auc, "partial.auc.focus")
  se <- sort(roc$sensitivities, decreasing = TRUE)
  sp <- sort(roc$specificities, decreasing = FALSE)

  if (identical(partial.auc, FALSE)) {
    suppressWarnings(polygon(c(sp, 0), c(se, 0), col = col, lty = lty, border = border, density = density, angle = angle, ...))
    return(invisible(NULL))
  }

  if (partial.auc.focus == "sensitivity") {
    x.all <- rev(se)
    y.all <- rev(sp)
  } else {
    x.all <- sp
    y.all <- se
  }
  # find the SEs and SPs in the interval
  in.window <- roc_utils_partial_auc_window(x.all, partial.auc)
  x.int <- x.all[in.window]
  y.int <- y.all[in.window]
  # if a limit is not exactly present in the curve, interpolate it
  upper <- roc_utils_interpolate_partial_auc_boundary(x.all, y.all, partial.auc[1])
  if (!is.null(upper)) {
    x.int <- c(x.int, upper$x)
    y.int <- c(y.int, upper$y)
  }
  lower <- roc_utils_interpolate_partial_auc_boundary(x.all, y.all, partial.auc[2])
  if (!is.null(lower)) {
    x.int <- c(lower$x, x.int)
    y.int <- c(lower$y, y.int)
  }
  # anchor to baseline
  x.int <- c(partial.auc[2], x.int, partial.auc[1])
  y.int <- c(0, y.int, 0)
  if (partial.auc.focus == "sensitivity") {
    # for SE, invert x and y again
    suppressWarnings(polygon(y.int, x.int, col = col, lty = lty, border = border, density = density, angle = angle, ...))
  } else {
    suppressWarnings(polygon(x.int, y.int, col = col, lty = lty, border = border, density = density, angle = angle, ...))
  }
  invisible(NULL)
}

roc_utils_draw_auc_text <- function(auc, ci = NULL, pattern = NULL, xy = NULL,
                                     adj = c(0, 1), col = par("col"), cex = par("cex"), ...) {
  if (is.null(xy)) {
    percent <- attr(auc, "percent")
    xy <- c(ifelse(percent, 50, .5), ifelse(percent, 50, .5))
  }
  label <- ggroc_auc_label(auc, ci, pattern)
  suppressWarnings(text(xy[1], xy[2], label, adj = adj, cex = cex, col = col, ...))
  invisible(label)
}

# polygon_auc(): the AUC (or partial AUC) area, matching
# plot.roc(auc.polygon = TRUE).

polygon_auc <- function(x, ...) {
  UseMethod("polygon_auc")
}

polygon_auc.auc <- function(x, col = "gainsboro", lty = par("lty"), density = NULL,
                             angle = 45, border = NULL, ...) {
  roc_utils_stop_if_no_device("polygon_auc")
  roc <- attr(x, "roc")
  roc_utils_draw_auc_polygon(roc, x, col = col, lty = lty, density = density, angle = angle, border = border, ...)
  invisible(x)
}

polygon_auc.roc <- function(x, ...) {
  roc_utils_stop_if_no_auc(x)
  polygon_auc(x$auc, ...)
  invisible(x)
}

polygon_auc.smooth.roc <- polygon_auc.roc

polygon_auc.multiclass.auc <- function(x, ...) {
  stop("polygon_auc() does not support 'multiclass.auc' objects.")
}

polygon_auc.mv.multiclass.auc <- polygon_auc.multiclass.auc

# polygon_max_auc(): the 100% (or partial) max-area rectangle, matching
# plot.roc(max.auc.polygon = TRUE).

polygon_max_auc <- function(x, ...) {
  UseMethod("polygon_max_auc")
}

polygon_max_auc.auc <- function(x, col = "#EEEEEE", ...) {
  roc_utils_stop_if_no_device("polygon_max_auc")
  roc_utils_draw_max_auc_polygon(x, col = col, ...)
  invisible(x)
}

polygon_max_auc.roc <- function(x, ...) {
  roc_utils_stop_if_no_auc(x)
  polygon_max_auc(x$auc, ...)
  invisible(x)
}

polygon_max_auc.smooth.roc <- polygon_max_auc.roc

polygon_max_auc.multiclass.auc <- function(x, ...) {
  stop("polygon_max_auc() does not support 'multiclass.auc' objects.")
}

polygon_max_auc.mv.multiclass.auc <- polygon_max_auc.multiclass.auc

# text() methods: the AUC (+ CI) label, matching plot.roc(print.auc = TRUE).

text.auc <- function(x, ci = NULL, pattern = NULL, xy = NULL, adj = c(0, 1),
                      col = par("col"), cex = par("cex"), ...) {
  roc_utils_stop_if_no_device("text.auc")
  invisible(roc_utils_draw_auc_text(x, ci = ci, pattern = pattern, xy = xy, adj = adj, col = col, cex = cex, ...))
}

text.ci.auc <- function(x, ...) {
  text.auc(attr(x, "auc"), ci = x, ...)
}

text.roc <- function(x, ...) {
  roc_utils_stop_if_no_auc(x)
  text.auc(x$auc, ci = x$ci, ...)
}

text.smooth.roc <- text.roc

text.multiclass.auc <- function(x, ...) {
  stop("text() does not support 'multiclass.auc' objects.")
}

text.mv.multiclass.auc <- text.multiclass.auc
