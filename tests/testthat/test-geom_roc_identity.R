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

test_that("geom_roc_identity colour defaults to the theme's ink and paper and can be set with color=", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  colour_of <- function(p) unique(ggplot2::ggplot_build(p)$data[[1]]$colour)
  p <- ggroc(r.s100b) + geom_roc_identity(r.s100b)
  expect_equal(colour_of(p), "#A8A8A8FF")
  expect_equal(colour_of(p + ggplot2::theme_grey(ink = "white", paper = "black")), "#575757FF")
  p <- ggroc(r.s100b) + geom_roc_identity(r.s100b, color = "lightblue")
  expect_warning(ggplot2::ggplot_build(p), NA)
  expect_equal(p@layers[[1]]$aes_params$colour, "lightblue")
})

test_that("geom_roc_identity adds a mapping to its own", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  p <- ggroc(r.s100b) + geom_roc_identity(r.s100b, mapping = ggplot2::aes(linetype = "chance", colour = "chance"))
  d <- ggplot2::layer_data(p, 1)
  expect_equal(c(d$linetype, d$colour), c("solid", "#F8766D"))
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
