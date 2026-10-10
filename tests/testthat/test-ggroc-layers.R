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

test_that("ggroc layers can be added with any ggplot2 version", {
  # No minimum version: the layers must not rely on ggplot2 4 S7 properties
  skip_if_not_installed("ggplot2")
  p <- ggroc(r.s100b)
  p <- p + geom_polygon_auc(r.s100b) + geom_roc_identity(r.s100b) + geom_text_auc(r.s100b)
  expect_equal(
    unname(vapply(p$layers, layer_geom_class, character(1))),
    c("GeomPolygon", "GeomSegment", "GeomLine", "GeomText")
  )
  expect_silent(ggplot2::ggplot_build(p))
  expect_true(pROC:::ggroc_legacy_axes_from_plot(ggroc(r.s100b, legacy.axes = TRUE)))
})
