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
  msg <- "Only 'best' 'thresholds' can be marked on a smoothed ROC curve: it has no thresholds."
  for (thresholds in list("all", "local maximas", "local", 0.5, c(0.5, 1))) {
    expect_error(geom_roc_threshold(s, thresholds = thresholds), msg, fixed = TRUE)
  }
  expect_s3_class(geom_roc_threshold(s, thresholds = "best"), "ggroc_layer")
  expect_s3_class(geom_roc_threshold(s, thresholds = "b"), "ggroc_layer")
})

test_that("geom_roc_threshold threshold labels don't depend on percent", {
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

test_that("geom_roc_threshold gives an empty layer when no threshold is in the partial AUC window", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  r <- roc(controls = c(1, 3), cases = c(2, 4), partial.auc = c(0.8, 0.6), quiet = TRUE)
  expect_warning(layer <- geom_roc_threshold(r), "No coordinates found")
  p <- ggroc(r) + layer
  expect_equal(nrow(p@layers[[match("GeomPoint", layer_geom_classes(p))]]$data), 0)
  expect_length(p@layers[[match("GeomText", layer_geom_classes(p))]]$data$label, 0)
  pdf(NULL)
  on.exit(dev.off())
  expect_error(print(p), NA)
})
