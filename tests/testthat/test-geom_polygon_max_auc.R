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

test_that("geom_polygon_max_auc fill defaults to the theme's ink and paper", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  fill_of <- function(p) unique(ggplot2::ggplot_build(p)$data[[1]]$fill)
  p <- ggroc(r.s100b) + geom_polygon_max_auc(r.s100b)
  expect_equal(fill_of(p), "#CCCCCCFF")
  expect_equal(fill_of(p + ggplot2::theme_grey(ink = "white", paper = "black")), "#333333FF")
  expect_equal(fill_of(ggroc(r.s100b) + geom_polygon_max_auc(r.s100b, fill = "red")), "red")
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

test_that("geom_polygon_max_auc gives a clear error on a curve without AUC", {
  r <- roc(aSAH$outcome, aSAH$s100b, auc = FALSE, quiet = TRUE)
  expect_error(geom_polygon_max_auc(r), "'data' has no 'auc'")
})
