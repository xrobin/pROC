context("geom_roc_identity")

test_that("geom_roc_identity is inserted behind the ROC curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) + geom_roc_identity(r.s100b)
  expect_equal(layer_geom_classes(p), c("GeomSegment", "GeomLine"))
})

test_that("geom_roc_identity sits in front of the AUC polygon", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) +
    geom_polygon_auc(r.s100b$auc) +
    geom_roc_identity(r.s100b)
  expect_equal(
    layer_geom_classes(p),
    c("GeomPolygon", "GeomSegment", "GeomLine")
  )
})

test_that("geom_roc_identity colour can be set with color=", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) + geom_roc_identity(r.s100b, color = "lightblue")
  expect_warning(ggplot2::ggplot_build(p), NA)
  expect_equal(p@layers[[1]]$aes_params$colour, "lightblue")
})

test_that("geom_roc_identity screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    print(ggroc(r.s100b) + geom_roc_identity(r.s100b))
  }
  expect_ggroc_doppelganger("geom_roc_identity.screenshot", test_screenshot)
})

test_that("geom_roc_identity works with percent and legacy.axes", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    print(
      ggroc(r.s100b.percent, legacy.axes = TRUE) +
        geom_roc_identity(r.s100b.percent)
    )
  }
  expect_ggroc_doppelganger("geom_roc_identity.percent.legacy.screenshot", test_screenshot)
})
