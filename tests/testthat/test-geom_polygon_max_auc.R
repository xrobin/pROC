context("geom_polygon_max_auc")

test_that("geom_polygon_max_auc is inserted behind the ROC curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) + geom_polygon_max_auc(r.s100b$auc)
  expect_equal(layer_geom_classes(p), c("GeomPolygon", "GeomLine"))
})

test_that("geom_polygon_max_auc sits behind geom_polygon_auc", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) +
    geom_polygon_max_auc(r.s100b$auc) +
    geom_polygon_auc(r.s100b$auc)
  expect_equal(
    layer_geom_classes(p),
    c("GeomPolygon", "GeomPolygon", "GeomLine")
  )
})

test_that("geom_polygon_max_auc screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    print(ggroc(r.s100b) + geom_polygon_max_auc(r.s100b$auc) + geom_polygon_auc(r.s100b$auc))
  }
  expect_ggroc_doppelganger("geom_polygon_max_auc.screenshot", test_screenshot)
})

test_that("geom_polygon_max_auc works with partial AUC", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    auc_sp <- auc(r.s100b, partial.auc = c(0.9, 1))
    print(ggroc(r.s100b) + geom_polygon_max_auc(auc_sp) + geom_polygon_auc(auc_sp))
  }
  expect_ggroc_doppelganger("geom_polygon_max_auc.partial.screenshot", test_screenshot)
})
