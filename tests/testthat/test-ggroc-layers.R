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

test_that("the layers of a curve take its colour unless one is given", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)
  # The colours hold whatever the bootstrap resamples, so no seed is needed
  ci_ses <- lapply(rocs, ci.se, specificities = seq(0, 1, 0.05), boot.n = 10, progress = "none")
  # Each layer, for one curve (i) or the list (NULL), and the aesthetic of its colour
  layers <- list(
    list(function(i, ...) geom_polygon_auc(if (is.null(i)) rocs else rocs[[i]], ...), "fill"),
    list(function(i, ...) geom_polygon_ci(if (is.null(i)) ci_ses else ci_ses[[i]], ...), "fill"),
    list(function(i, ...) geom_errorbar_ci(if (is.null(i)) ci_ses else ci_ses[[i]], ...), "colour")
  )
  # The colours of the layer, per curve, besides the ROC line
  colours_of <- function(p, aesthetic) {
    b <- ggplot2::ggplot_build(p)
    roc_line <- which(layer_geom_classes(p) == "GeomLine")
    unique(unlist(lapply(b$data[-roc_line], function(d) unique(d[[aesthetic]]))))
  }
  for (layer in layers) {
    make_layer <- layer[[1]]
    aesthetic <- layer[[2]]
    as_layer_colour <- function(colour) if (aesthetic == "fill") ggplot2::alpha(colour, 0.25) else colour
    # The colour mapped to the curves, with the scale of the plot
    p <- ggroc(rocs) + ggplot2::facet_grid(. ~ name) + make_layer(NULL) +
      ggplot2::scale_colour_brewer(palette = "Dark2")
    line <- ggplot2::layer_data(p, which(layer_geom_classes(p) == "GeomLine"))
    expect_equal(colours_of(p, aesthetic), as_layer_colour(unname(as.vector(tapply(line$colour, line$PANEL, unique)))))
    # It adds nothing to the legend
    expect_true(all(vapply(p@layers, function(l) isFALSE(l$show.legend) || inherits(l$geom, "GeomLine"), logical(1))))
    # The fixed colour of the ROC line
    expect_equal(colours_of(ggroc(r.s100b, colour = "red") + make_layer(1), aesthetic), as_layer_colour("red"))
    expect_equal(colours_of(ggroc(rocs, aes = "linetype", colour = "red") + make_layer(NULL), aesthetic), as_layer_colour("red"))
    # Neither: the theme's ink
    p <- ggroc(r.s100b) + make_layer(1)
    expect_equal(colours_of(p, aesthetic), as_layer_colour("black"))
    if (aesthetic == "fill") {
      expect_equal(colours_of(p + ggplot2::theme_grey(ink = "white", paper = "black"), aesthetic), as_layer_colour("white"))
    }
    # An explicit colour wins
    explicit <- stats::setNames(list("grey30"), aesthetic)
    expect_equal(colours_of(ggroc(rocs) + do.call(make_layer, c(list(NULL), explicit)), aesthetic), "grey30")
    # A mapped colour wins too: the default hue scale of 'fill' or
    # 'colour', without the alpha of the default fill
    mapped <- do.call(ggplot2::aes, stats::setNames(list(quote(name)), aesthetic))
    p <- ggroc(rocs) + ggplot2::facet_grid(. ~ name) + make_layer(NULL, mapping = mapped)
    line <- ggplot2::layer_data(p, which(layer_geom_classes(p) == "GeomLine"))
    expect_equal(colours_of(p, aesthetic), unname(as.vector(tapply(line$colour, line$PANEL, unique))))
  }
  # Thresholds stay in the theme's ink, as print.thres.col in plot.roc
  expect_equal(colours_of(ggroc(rocs) + geom_roc_threshold(rocs), "colour"), "black")
})

test_that("the layers of a list add a mapping to their own", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)
  # The mapped values hold whatever the bootstrap resamples, so no seed is needed
  ci_ses <- lapply(rocs, ci.se, specificities = seq(0, 1, 0.05), boot.n = 10, progress = "none")
  ci_ths <- lapply(rocs, ci.thresholds, thresholds = "best", boot.n = 10, progress = "none")
  mapping <- ggplot2::aes(linetype = name, colour = name)
  layers <- list(
    geom_polygon_max_auc(rocs, mapping = mapping),
    geom_polygon_auc(rocs, mapping = mapping),
    geom_polygon_ci(ci_ses, mapping = mapping),
    geom_errorbar_ci(ci_ses, mapping = mapping),
    geom_errorbar_ci(ci_ths, mapping = mapping),
    geom_roc_threshold(rocs, mapping = mapping),
    geom_text_auc(rocs, mapping = mapping)
  )
  for (layer in layers) {
    p <- ggroc(rocs, aes = c("linetype", "colour")) + layer
    b <- ggplot2::ggplot_build(p)
    roc_line <- which(layer_geom_classes(p) == "GeomLine")
    line <- b$data[[roc_line]]
    line_style <- unique(paste(line$colour, line$linetype))
    # Each curve's layer has the colour and linetype of the curve
    styles <- unique(unlist(lapply(b$data[-roc_line], function(d) paste(d$colour, d$linetype))))
    expect_equal(styles, line_style)
  }
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
