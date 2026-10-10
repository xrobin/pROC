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
