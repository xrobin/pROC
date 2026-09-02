context("geom_polygon_auc")

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

test_that("geom_polygon_auc stays behind the ROC line when other layers are on top", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) +
    ggplot2::geom_point() +
    geom_polygon_auc(r.s100b$auc)
  expect_equal(
    layer_geom_classes(p),
    c("GeomPolygon", "GeomLine", "GeomPoint")
  )
})

test_that("geom_polygon_auc follows the parent plot's axes", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b, legacy.axes = TRUE) + geom_polygon_auc(r.s100b$auc)
  expect_match(
    rlang::as_label(p@layers[[1]]$mapping$x),
    "1-specificity",
    fixed = TRUE
  )
})

test_that("geom_polygon_auc works on a ggroc.list plot", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(list(r.s100b, smooth(r.s100b)), colour = "red", linetype = 2, linewidth = 3) +
    geom_polygon_auc(r.s100b$auc)
  expect_error(ggplot2::ggplot_build(p), NA)
})

test_that("geom_polygon_auc closes a smooth ROC below the identity line", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  sm <- smooth(r.s100b)
  p <- ggroc(sm) + geom_polygon_auc(sm)
  d <- p@layers[[1]]$data
  expect_equal(as.numeric(d$specificity[nrow(d)]), 0)
  expect_equal(as.numeric(d$sensitivity[nrow(d)]), 0)
})

test_that("geom_polygon_auc works", {
  test_geom_polygon_auc_screenshot <- function() {
    print(ggroc(r.s100b) + geom_polygon_auc(r.s100b$auc))
  }
  expect_ggroc_doppelganger("geom_polygon_auc.screenshot", test_geom_polygon_auc_screenshot)
})

test_that("geom_polygon_auc works with percent and legacy.axes = TRUE", {
  test_geom_polygon_auc_percent_legacy_screenshot <- function() {
    print(ggroc(r.s100b.percent, legacy.axes = TRUE) + geom_polygon_auc(r.s100b.percent$auc))
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
