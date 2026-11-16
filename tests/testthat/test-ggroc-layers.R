context("ggroc layers")

test_that("background layers sit immediately before the ROC line", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  line <- ggplot2::geom_line()
  poly1 <- ggplot2::geom_polygon()
  poly2 <- ggplot2::geom_polygon()
  point <- ggplot2::geom_point()
  result <- pROC:::ggroc_insert_last_before_roc_line(list(poly1, line, point, poly2))
  expect_equal(
    unname(vapply(result, layer_geom_class, character(1))),
    c("GeomPolygon", "GeomPolygon", "GeomLine", "GeomPoint")
  )
})

test_that("ggroc_legacy_axes_from_plot reads the parent x mapping", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  expect_false(pROC:::ggroc_legacy_axes_from_plot(ggroc(r.s100b)))
  expect_true(pROC:::ggroc_legacy_axes_from_plot(ggroc(r.s100b, legacy.axes = TRUE)))
})

test_that("layers stop with a clear error on ggplot2 < 4.0.0", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  skip_if_not_installed("testthat", minimum_version = "3.1.7") # local_mocked_bindings()
  local_mocked_bindings(ggroc_ggplot2_version = function() "3.5.2", .package = "pROC")
  expect_error(ggroc(r.s100b) + geom_roc_identity(r.s100b), "require ggplot2 >= 4.0.0, but ggplot2 3.5.2 is loaded")
  expect_error(ggroc(r.s100b) + geom_polygon_auc(r.s100b), "require ggplot2 >= 4.0.0")
})

test_that("layers given a list draw each curve in its own facet", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)
  panels_of <- function(layer) {
    p <- ggroc(rocs) + ggplot2::facet_grid(. ~ name) + layer
    b <- ggplot2::ggplot_build(p)
    expect_equal(as.character(b$layout$layout$name), names(rocs))
    roc_line <- which(layer_geom_classes(p) == "GeomLine")
    lapply(b$data[-roc_line], function(d) unique(as.integer(d$PANEL)))
  }
  # The panels hold whatever the bootstrap resamples, so no seed is needed
  ci_ses <- lapply(rocs, ci.se, specificities = seq(0, 1, 0.05), boot.n = 10, progress = "none")
  ci_ths <- lapply(rocs, ci.thresholds, thresholds = "best", boot.n = 10, progress = "none")
  expect_equal(panels_of(geom_polygon_max_auc(rocs)), list(1L, 2L, 3L))
  expect_equal(panels_of(geom_polygon_ci(ci_ses)), list(1L, 2L, 3L))
  expect_equal(panels_of(geom_errorbar_ci(ci_ses)), list(1L, 2L, 3L))
  # One layer per bar direction or per point and label
  expect_equal(panels_of(geom_errorbar_ci(ci_ths)), list(1L, 1L, 2L, 2L, 3L, 3L))
  expect_equal(panels_of(geom_roc_threshold(rocs)), list(1L, 1L, 2L, 2L, 3L, 3L))
})

test_that("layers given a list work in facets", {
  test_ggroc_layers_list_facet_screenshot <- function() {
    rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)
    print(ggroc(rocs) + ggplot2::facet_grid(. ~ name) +
      geom_polygon_max_auc(rocs) +
      geom_polygon_auc(rocs) +
      geom_roc_threshold(rocs))
  }
  expect_ggroc_doppelganger("ggroc.layers.list.facet.screenshot", test_ggroc_layers_list_facet_screenshot)
})
