context("geom_roc_threshold")

test_that("geom_roc_threshold adds a point and a label in front of the ROC", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) + geom_roc_threshold(r.s100b, thresholds = "best")
  expect_equal(layer_geom_classes(p), c("GeomLine", "GeomPoint", "GeomText"))
})

test_that("geom_roc_threshold screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    print(ggroc(r.s100b) + geom_roc_threshold(r.s100b, thresholds = "best"))
  }
  expect_ggroc_doppelganger("geom_roc_threshold.screenshot", test_screenshot)
})
