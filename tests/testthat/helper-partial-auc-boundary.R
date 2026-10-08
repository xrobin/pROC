# Helpers for test-partial-auc-boundary.R.
#
# auc(), coords(x = "all") and the partial-AUC polygon (plot.auc.R) each
# need the curve value at the two partial.auc bounds, and historically each
# interpolated it with its own code. These helpers provide an independent
# reference ("oracle") implementation plus accessors that extract, for one
# roc object and one partial.auc window, what each of the three consumers
# produces.

# trace()'s tracer always runs in the global environment: go through an
# explicit, dedicated global binding instead of a local `<<-`.
partial_auc_capture_polygon <- function(expr) {
  assign(".test_pab_captured_polygon", NULL, envir = globalenv())
  suppressMessages(trace(graphics::polygon,
    tracer = quote(assign(".test_pab_captured_polygon", list(x = x, y = y), envir = globalenv())),
    print = FALSE
  ))
  on.exit({
    suppressMessages(try(untrace(graphics::polygon), silent = TRUE))
    rm(".test_pab_captured_polygon", envir = globalenv())
  })
  force(expr)
  get(".test_pab_captured_polygon", envir = globalenv())
}

partial_auc_shoelace <- function(x, y) {
  abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
}

partial_auc_trapezoid <- function(x, y) {
  if (length(x) < 2) {
    return(0)
  }
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}

# Orient (x, y) the way pROC walks the curve for the given focus: x
# increasing, ties in x kept in pROC's order (y decreasing).
partial_auc_orient <- function(sp, se, focus) {
  if (focus == "sensitivity") {
    x <- se
    y <- sp
  } else {
    x <- sp
    y <- se
  }
  o <- order(x, -y)
  list(x = x[o], y = y[o])
}

# Independent reference: the vertices (and area) of the curve restricted to
# [lo, hi]. Points already in the window are kept as-is (with exact
# equality, as the package does); a missing bound is linearly interpolated
# between the two points straddling it.
partial_auc_oracle <- function(x, y, lo, hi) {
  interp <- function(b) {
    i <- match(TRUE, x > b)
    y[i - 1] + (b - x[i - 1]) / (x[i] - x[i - 1]) * (y[i] - y[i - 1])
  }
  inside <- x >= lo & x <= hi
  xs <- x[inside]
  ys <- y[inside]
  if (!(lo %in% xs)) {
    xs <- c(lo, xs)
    ys <- c(interp(lo), ys)
  }
  if (!(hi %in% xs)) {
    xs <- c(xs, hi)
    ys <- c(ys, interp(hi))
  }
  list(x = xs, y = ys, area = partial_auc_trapezoid(xs, ys))
}

# What the three consumers return for `roc` and the window `partial.auc`
# (given in the roc's own scale, i.e. percent if roc$percent).
# The result is expressed in the roc's scale, x = focus axis, y = other axis.
partial_auc_consumers <- function(roc, partial.auc, focus) {
  scale <- ifelse(roc$percent, 100, 1)
  auc.value <- as.numeric(auc(roc, partial.auc = partial.auc, partial.auc.focus = focus, partial.auc.correct = FALSE))

  roc.partial <- roc
  roc.partial$auc <- auc(roc, partial.auc = partial.auc, partial.auc.focus = focus, partial.auc.correct = FALSE)
  coords.res <- suppressWarnings(coords(roc.partial, "all", ret = c("specificity", "sensitivity"), transpose = FALSE))
  if (is.null(coords.res)) {
    coords.vertices <- NULL
  } else {
    sp <- coords.res[, "specificity"]
    se <- coords.res[, "sensitivity"]
    coords.vertices <- partial_auc_orient(sp, se, focus)
  }

  polygon.res <- partial_auc_capture_polygon({
    plot(roc.partial)
    polygon_auc(roc.partial)
  })
  # polygon() is called with (x, y) swapped for sensitivity focus and
  # anchored to the baseline: swap back, the area is what we compare.
  if (focus == "sensitivity") {
    polygon.x <- polygon.res$y
    polygon.y <- polygon.res$x
  } else {
    polygon.x <- polygon.res$x
    polygon.y <- polygon.res$y
  }
  list(
    auc = auc.value,
    coords = coords.vertices,
    polygon = list(x = polygon.x, y = polygon.y),
    # the polygon is drawn in plot units: area in scale^2, auc() in scale^2 / scale
    polygon.area = partial_auc_shoelace(polygon.x, polygon.y) / scale
  )
}

# A deterministic collection of curves exercising the edge cases listed in
# the issue: ties in the predictor, ordered (few-level) predictors,
# degenerate curves, percent scale and reversed direction.
partial_auc_curves <- function() {
  data(aSAH, envir = environment())
  ties <- roc(c(0, 0, 0, 1, 1, 1, 0, 1), c(1, 2, 3, 3, 4, 5, 5, 6), quiet = TRUE)
  list(
    s100b = roc(aSAH$outcome, aSAH$s100b, quiet = TRUE),
    s100b.percent = roc(aSAH$outcome, aSAH$s100b, percent = TRUE, quiet = TRUE),
    s100b.direction = roc(aSAH$outcome, aSAH$s100b, direction = ">", quiet = TRUE),
    ndka = roc(aSAH$outcome, aSAH$ndka, quiet = TRUE),
    wfns.ordered = roc(aSAH$outcome, aSAH$wfns, quiet = TRUE),
    wfns.numeric = roc(aSAH$outcome, as.numeric(aSAH$wfns), quiet = TRUE),
    ties = ties,
    ties.percent = roc(c(0, 0, 0, 1, 1, 1, 0, 1), c(1, 2, 3, 3, 4, 5, 5, 6), percent = TRUE, quiet = TRUE),
    two.points = roc(c(0, 1), c(1, 2), quiet = TRUE),
    constant = roc(c(0, 0, 1, 1), c(1, 1, 1, 1), quiet = TRUE)
  )
}

# A set of partial.auc windows for a curve and focus, in the roc's scale:
# exactly on points, strictly between points, both bounds inside the same
# segment, one bound on a point, and the full range.
partial_auc_windows <- function(roc, focus) {
  scale <- ifelse(roc$percent, 100, 1)
  x <- partial_auc_orient(roc$specificities, roc$sensitivities, focus)$x
  pts <- sort(unique(x))
  n <- length(pts)
  windows <- list(
    full = c(max(pts), min(pts)),
    between = c(0.95, 0.9) * scale,
    between.low = c(0.45, 0.05) * scale
  )
  if (n >= 3) {
    windows$on.points <- c(pts[n - 1], pts[2])
    windows$hi.on.point <- c(pts[n - 1], mean(pts[1:2]))
    windows$lo.on.point <- c(mean(pts[(n - 1):n]), pts[2])
  }
  if (n >= 2) {
    mid <- floor(n / 2)
    seg <- c(pts[mid], pts[mid + 1])
    windows$same.segment <- c(seg[1] + 0.75 * diff(seg), seg[1] + 0.25 * diff(seg))
  }
  windows
}

# coords() results restricted to a partial AUC window for the "best",
# "local maximas" (and, for smooth curves, "all") special values. Only the
# specificity, sensitivity and threshold columns are recorded.
partial_auc_coords_snapshot <- function() {
  curves <- partial_auc_curves()
  coords_matrix <- function(roc, x, ...) {
    res <- tryCatch(
      suppressWarnings(coords(roc, x, ret = c("specificity", "sensitivity", "threshold"), transpose = FALSE, ...)),
      error = function(e) NULL
    )
    if (is.null(res)) NULL else unname(as.matrix(res))
  }
  out <- list()
  for (curve.name in c("s100b", "s100b.percent", "ties", "wfns.ordered")) {
    for (focus in c("specificity", "sensitivity")) {
      roc <- curves[[curve.name]]
      windows <- partial_auc_windows(roc, focus)
      for (window.name in intersect(names(windows), c("between", "on.points", "same.segment", "hi.on.point"))) {
        roc.partial <- roc
        roc.partial$auc <- auc(roc, partial.auc = windows[[window.name]], partial.auc.focus = focus, partial.auc.correct = FALSE)
        id <- paste(curve.name, focus, window.name, sep = "/")
        out[paste(id, "local maximas", sep = "/")] <- list(coords_matrix(roc.partial, "local maximas"))
        out[paste(id, "best youden", sep = "/")] <- list(coords_matrix(roc.partial, "best", best.method = "youden"))
        out[paste(id, "best closest.topleft", sep = "/")] <- list(coords_matrix(roc.partial, "best", best.method = "closest.topleft"))
      }
    }
  }
  smooth.roc <- smooth(curves$s100b, n = 20)
  for (focus in c("specificity", "sensitivity")) {
    for (bounds in list(c(0.99, 0.9), c(0.93, 0.91))) {
      smooth.partial <- smooth.roc
      smooth.partial$auc <- auc(smooth.roc, partial.auc = bounds, partial.auc.focus = focus, partial.auc.correct = FALSE)
      id <- paste("smooth", focus, paste(bounds, collapse = "-"), sep = "/")
      out[paste(id, "all", sep = "/")] <- list(coords_matrix(smooth.partial, "all"))
      out[paste(id, "best youden", sep = "/")] <- list(coords_matrix(smooth.partial, "best", best.method = "youden"))
      out[paste(id, "best closest.topleft", sep = "/")] <- list(coords_matrix(smooth.partial, "best", best.method = "closest.topleft"))
    }
  }
  out
}
