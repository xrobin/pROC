context("geom_roc_threshold")

test_that("geom_roc_threshold adds a point and a label in front of the ROC", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) + geom_roc_threshold(r.s100b, thresholds = "best")
  expect_equal(layer_geom_classes(p), c("GeomLine", "GeomPoint", "GeomText"))
})

test_that("geom_roc_threshold labels ordered levels, not integer codes", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  wfns.txt <- aSAH$wfns
  levels(wfns.txt) <- c("very low", "low", "medium", "high", "very high")
  r <- roc(aSAH$outcome, wfns.txt, quiet = TRUE)
  p <- ggroc(r) + geom_roc_threshold(r, thresholds = "best")
  text_i <- match("GeomText", layer_geom_classes(p))
  labels <- p@layers[[text_i]]$data$label
  expect_match(labels, "high", fixed = TRUE)
  expect_false(grepl("^[0-9]", labels))
})

test_that("geom_roc_threshold numeric labels still use midpoint format", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) + geom_roc_threshold(r.s100b, thresholds = "best")
  text_i <- match("GeomText", layer_geom_classes(p))
  labels <- p@layers[[text_i]]$data$label
  expect_match(labels, "^[0-9.]+ \\(")
})

test_that("geom_roc_threshold screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    print(ggroc(r.s100b) + geom_roc_threshold(r.s100b, thresholds = "best"))
  }
  expect_ggroc_doppelganger("geom_roc_threshold.screenshot", test_screenshot)
})

test_that("geom_roc_threshold gives clear errors for thresholds a smoothed curve does not have", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  s <- smooth(r.s100b)
  expect_error(geom_roc_threshold(s, thresholds = "local maximas"), "'local maximas' 'thresholds' unsupported on a smoothed ROC curve")
  expect_error(geom_roc_threshold(s, thresholds = "local"), "'local maximas' 'thresholds' unsupported on a smoothed ROC curve")
  expect_error(geom_roc_threshold(s, thresholds = 0.5), "Numeric 'thresholds' unsupported on a smoothed ROC curve")
  expect_error(geom_roc_threshold(s, thresholds = c(0.5, 1)), "Numeric 'thresholds' unsupported on a smoothed ROC curve")
  expect_s3_class(geom_roc_threshold(s, thresholds = "best"), "ggroc_layer")
})

test_that("geom_roc_threshold percent labels keep 3 decimals for the threshold", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  r.percent <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE, quiet = TRUE)
  p <- ggroc(r.percent) + geom_roc_threshold(r.percent, thresholds = "best")
  text_i <- match("GeomText", layer_geom_classes(p))
  expect_identical(p@layers[[text_i]]$data$label, "0.205 (80.6%, 63.4%)")
})

test_that("geom_roc_threshold('all') marks only real thresholds on a partial AUC curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  r <- roc(aSAH$outcome, aSAH$s100b, partial.auc = c(1, 0.85), quiet = TRUE)
  co <- coords(r, "all")
  p <- ggroc(r) + geom_roc_threshold(r, thresholds = "all")
  labels <- p@layers[[match("GeomText", layer_geom_classes(p))]]$data$label
  expect_length(labels, sum(!is.na(co$threshold)))
  expect_false(any(startsWith(labels, "NA")))
  expect_equal(nrow(p@layers[[match("GeomPoint", layer_geom_classes(p))]]$data), sum(!is.na(co$threshold)))
})

test_that("geom_roc_threshold uses best.method and best.weights as plot.roc", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  point.data <- function(layer) {
    p <- ggroc(r.wfns) + layer
    p@layers[[match("GeomPoint", layer_geom_classes(p))]]$data
  }
  layer <- geom_roc_threshold(r.wfns, thresholds = "best", best.method = "closest.topleft")
  # not passed to the geoms ("Ignoring unknown parameters")
  expect_no_warning(d <- point.data(layer))
  expect_identical(as.character(d$threshold), as.character(coords(r.wfns, "best", best.method = "closest.topleft")$threshold))
  layer <- geom_roc_threshold(r.wfns, thresholds = "best", best.weights = c(5, 0.2))
  expect_identical(as.character(point.data(layer)$threshold), as.character(coords(r.wfns, "best", best.weights = c(5, 0.2))$threshold))
  # the default is still Youden's point
  expect_identical(as.character(point.data(geom_roc_threshold(r.wfns))$threshold), as.character(coords(r.wfns, "best")$threshold))
})
