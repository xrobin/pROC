context("geom_polygon_auc")

layer_geom_class <- function(layer) {
  class(layer$geom)[[1]]
}

layer_geom_classes <- function(plot) {
  unname(vapply(plot@layers, layer_geom_class, character(1)))
}

test_that("geom_polygon_auc is inserted behind the ROC curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) + geom_polygon_auc(r.s100b$auc)
  expect_equal(
    layer_geom_classes(p),
    c("GeomPolygon", "GeomLine")
  )
})

test_that("partial AUC polygons stay behind the ROC curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  auc_sp <- auc(r.s100b, partial.auc = c(0.8, 0.9), partial.auc.focus = "sp")
  auc_se <- auc(r.s100b, partial.auc = c(0.8, 0.9), partial.auc.focus = "se")
  p <- ggroc(r.s100b) + geom_polygon_auc(auc_se) + geom_polygon_auc(auc_sp)
  expect_equal(
    layer_geom_classes(p),
    c("GeomPolygon", "GeomPolygon", "GeomLine")
  )
})

test_that("geom_polygon_auc works", {
  test_geom_polygon_auc_screenshot <- function() {
    print(ggroc(r.s100b) + geom_polygon_auc(r.s100b$auc))
  }
  expect_ggroc_doppelganger("geom_polygon_auc.screenshot", test_geom_polygon_auc_screenshot)
})

test_that("geom_polygon_auc works with percent and legacy.axes = TRUE", {
  test_geom_polygon_auc_percent_legacy_screenshot <- function() {
    print(ggroc(r.s100b.percent, legacy.axes = TRUE) + geom_polygon_auc(r.s100b.percent$auc, legacy.axes = TRUE))
  }
  expect_ggroc_doppelganger("geom_polygon_auc.percent.legacy.screenshot", test_geom_polygon_auc_percent_legacy_screenshot)
})


test_that("geom_polygon_auc works with percent and legacy.axes = TRUE", {
  test_geom_polygon_auc_partial_screenshot <- function() {
    auc_sp <- auc(r.s100b, partial.auc = c(0.8, 0.9), partial.auc.focus = "sp")
    auc_se <- auc(r.s100b, partial.auc = c(0.8, 0.9), partial.auc.focus = "se")
    print(ggroc(r.s100b) +
      geom_polygon_auc(auc_se) +
      geom_polygon_auc(auc_sp))
  }
  expect_ggroc_doppelganger("geom_polygon_auc.partial.screenshot", test_geom_polygon_auc_partial_screenshot)
})
