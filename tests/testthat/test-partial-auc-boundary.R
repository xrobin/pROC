context("partial AUC boundary interpolation")

# auc(), coords(x = "all") and the partial-AUC polygon must agree on where a
# partial-AUC region starts and ends, for any curve and any window. These
# tests check them against an independent reference implementation
# (partial_auc_oracle) and against values recorded before the shared
# boundary-interpolation primitive was introduced (issue #145).
# See helper-partial-auc-boundary.R.

curves <- partial_auc_curves()
foci <- c("specificity", "sensitivity")

# Run `code(curve.name, focus, window.name, roc, bounds, oracle)` on every
# curve / focus / window combination with a non-empty window.
for_each_window <- function(code) {
  for (curve.name in names(curves)) {
    roc <- curves[[curve.name]]
    for (focus in foci) {
      windows <- partial_auc_windows(roc, focus)
      for (window.name in names(windows)) {
        bounds <- windows[[window.name]]
        if (bounds[1] <= bounds[2]) next
        oriented <- partial_auc_orient(roc$specificities, roc$sensitivities, focus)
        oracle <- partial_auc_oracle(oriented$x, oriented$y, bounds[2], bounds[1])
        code(curve.name, focus, window.name, roc, bounds, oracle)
      }
    }
  }
}

test_that("the window collection covers every case from the issue", {
  ids <- character()
  for_each_window(function(curve.name, focus, window.name, ...) {
    ids <<- c(ids, paste(curve.name, focus, window.name, sep = "/"))
  })
  expect_identical(ids, names(partial_auc_boundary_expected_auc))
  for (kind in c("full", "between", "on.points", "hi.on.point", "lo.on.point", "same.segment")) {
    expect_true(any(grepl(paste0("/", kind, "$"), ids)), info = kind)
  }
  expect_true(any(grepl("percent", ids)))
  expect_true(any(grepl("direction", ids)))
  expect_true(any(grepl("ties", ids)))
  expect_true(any(grepl("two.points", ids)))
})

test_that("auc(), coords() and polygon_auc() all match the reference implementation", {
  pdf(NULL)
  on.exit(dev.off())
  for_each_window(function(curve.name, focus, window.name, roc, bounds, oracle) {
    info <- paste(curve.name, focus, window.name)
    scale <- ifelse(roc$percent, 100, 1)
    res <- partial_auc_consumers(roc, bounds, focus)
    expect_equal(res$auc, oracle$area / scale, info = info)
    expect_equal(res$polygon.area, oracle$area / scale, info = info)
    expect_false(is.null(res$coords), info = info)
    expect_equal(res$coords$x, oracle$x, info = info)
    expect_equal(res$coords$y, oracle$y, info = info)
    # The polygon is anchored to the baseline: vertices on the focus axis
    # span exactly the window.
    expect_equal(range(res$polygon$x), range(bounds), info = info)
  })
})

test_that("auc(), coords() and polygon_auc() agree with each other on random curves", {
  pdf(NULL)
  on.exit(dev.off())
  set.seed(145)
  checked <- 0
  while (checked < 40) {
    n.obs <- sample(c(8, 20, 60), 1)
    response <- rbinom(n.obs, 1, 0.5)
    if (length(unique(response)) < 2) next
    # round() creates ties in the predictor, hence vertical/horizontal runs
    predictor <- round(rnorm(n.obs) + response, sample(0:2, 1))
    roc <- roc(response, predictor, quiet = TRUE)
    focus <- sample(foci, 1)
    x <- partial_auc_orient(roc$specificities, roc$sensitivities, focus)$x
    pts <- unique(x)
    if (length(pts) < 3) next
    bounds <- switch(sample(c("random", "on.points", "hi.on.point", "lo.on.point"), 1),
      random = runif(2),
      on.points = sample(pts, 2),
      hi.on.point = c(sample(pts, 1), runif(1)),
      lo.on.point = c(runif(1), sample(pts, 1))
    )
    bounds <- sort(bounds, decreasing = TRUE)
    if (bounds[1] - bounds[2] < 1e-6 || bounds[1] > max(pts) || bounds[2] < min(pts)) next
    oriented <- partial_auc_orient(roc$specificities, roc$sensitivities, focus)
    oracle <- partial_auc_oracle(oriented$x, oriented$y, bounds[2], bounds[1])
    res <- partial_auc_consumers(roc, bounds, focus)
    info <- paste("random case", checked, focus, paste(format(bounds, digits = 17), collapse = " / "))
    expect_equal(res$auc, oracle$area, info = info)
    expect_equal(res$polygon.area, oracle$area, info = info)
    expect_equal(res$coords$x, oracle$x, info = info)
    expect_equal(res$coords$y, oracle$y, info = info)
    checked <- checked + 1
  }
})

test_that("a bound exactly on an existing point does not add a vertex, in all three", {
  pdf(NULL)
  on.exit(dev.off())
  roc <- curves$s100b
  pts <- sort(unique(roc$specificities))
  bounds <- c(pts[20], pts[8])
  res <- partial_auc_consumers(roc, bounds, "specificity")
  in.window <- roc$specificities >= bounds[2] & roc$specificities <= bounds[1]
  expect_equal(length(res$coords$x), sum(in.window))
  # polygon: window vertices + the two baseline anchors
  expect_equal(length(res$polygon$x), sum(in.window) + 2)
})

test_that("a bound strictly between two points adds exactly one interpolated vertex", {
  pdf(NULL)
  on.exit(dev.off())
  roc <- curves$s100b
  pts <- sort(unique(roc$specificities))
  bounds <- c(mean(pts[20:21]), pts[8])
  res <- partial_auc_consumers(roc, bounds, "specificity")
  in.window <- roc$specificities >= bounds[2] & roc$specificities <= bounds[1]
  expect_equal(length(res$coords$x), sum(in.window) + 1)
  expect_equal(max(res$coords$x), bounds[1])
  expect_equal(length(res$polygon$x), sum(in.window) + 1 + 2)
})

test_that("a floating-point near-miss of an existing point is interpolated, not matched (exact equality)", {
  pdf(NULL)
  on.exit(dev.off())
  roc <- curves$s100b
  p <- sort(unique(roc$specificities))[20]
  # exactly on the point: no interpolation
  on.point <- partial_auc_consumers(roc, c(p, 0.5), "specificity")
  # one ulp away: a different, interpolated, boundary, but the same region
  eps <- .Machine$double.eps
  near <- partial_auc_consumers(roc, c(p + eps, 0.5), "specificity")
  expect_equal(near$auc, on.point$auc, tolerance = 1e-10)
  expect_equal(max(near$coords$x), p + eps)
  expect_equal(max(on.point$coords$x), p)
  expect_gt(length(near$coords$x), length(on.point$coords$x))
})

test_that("auc(), coords() and polygon_auc() agree on smooth.roc curves", {
  pdf(NULL)
  on.exit(dev.off())
  smooth.roc <- smooth(curves$s100b, n = 20)
  for (focus in foci) {
    for (bounds in list(c(0.99, 0.9), c(0.93, 0.91))) {
      info <- paste(focus, paste(bounds, collapse = "-"))
      smooth.roc$auc <- auc(smooth.roc, partial.auc = bounds, partial.auc.focus = focus, partial.auc.correct = FALSE)
      coords.res <- suppressWarnings(coords(smooth.roc, "all", ret = c("specificity", "sensitivity"), transpose = FALSE))
      polygon.res <- partial_auc_capture_polygon({
        plot(smooth.roc)
        polygon_auc(smooth.roc)
      })
      expect_equal(partial_auc_shoelace(polygon.res$x, polygon.res$y), as.numeric(smooth.roc$auc), info = info)
      focus.axis <- coords.res[, focus]
      expect_equal(range(focus.axis), range(bounds), info = info)
      oriented <- partial_auc_orient(coords.res[, "specificity"], coords.res[, "sensitivity"], focus)
      expect_equal(partial_auc_trapezoid(oriented$x, oriented$y), as.numeric(smooth.roc$auc), info = info)
    }
  }
})

test_that("percent and fraction scales describe the same region", {
  pdf(NULL)
  on.exit(dev.off())
  for (focus in foci) {
    fraction <- partial_auc_consumers(curves$ties, c(0.9, 0.3), focus)
    percent <- partial_auc_consumers(curves$ties.percent, c(90, 30), focus)
    expect_equal(percent$auc, fraction$auc * 100, info = focus)
    expect_equal(percent$coords$x, fraction$coords$x * 100, info = focus)
    expect_equal(percent$coords$y, fraction$coords$y * 100, info = focus)
  }
})

test_that("auc() with partial.auc.correct is the McClish transform of the same region", {
  roc <- curves$s100b
  bounds <- c(0.95, 0.9)
  raw <- as.numeric(auc(roc, partial.auc = bounds, partial.auc.correct = FALSE))
  corrected <- as.numeric(auc(roc, partial.auc = bounds, partial.auc.correct = TRUE))
  min.area <- diff(rev(bounds)) * (1 - mean(bounds)) # area under the chance line
  max.area <- diff(rev(bounds))
  expect_equal(corrected, 0.5 * (1 + (raw - min.area) / (max.area - min.area)))
})

test_that("recorded auc() values are unchanged by the refactor", {
  expected <- partial_auc_boundary_expected_auc
  actual <- numeric()
  for_each_window(function(curve.name, focus, window.name, roc, bounds, oracle) {
    actual[paste(curve.name, focus, window.name, sep = "/")] <<- as.numeric(
      auc(roc, partial.auc = bounds, partial.auc.focus = focus, partial.auc.correct = FALSE)
    )
  })
  expect_identical(names(actual), names(expected))
  # Default expect_equal() tolerance, no looser: the refactor must not change
  # the load-bearing auc() output.
  expect_equal(actual, expected)
})

test_that("recorded coords() and polygon vertices are unchanged by the refactor", {
  pdf(NULL)
  on.exit(dev.off())
  for (id in names(partial_auc_boundary_expected_vertices)) {
    parts <- strsplit(id, "/")[[1]]
    roc <- curves[[parts[1]]]
    focus <- parts[2]
    bounds <- partial_auc_windows(roc, focus)[[parts[3]]]
    res <- partial_auc_consumers(roc, bounds, focus)
    expected <- partial_auc_boundary_expected_vertices[[id]]
    expect_equal(lapply(res$coords, unname), expected$coords, info = id)
    expect_equal(lapply(res$polygon, unname), expected$polygon, info = id)
  }
})
