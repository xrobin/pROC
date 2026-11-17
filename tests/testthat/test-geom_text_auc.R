context("geom_text_auc")

test_that("geom_text_auc screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    print(ggroc(r.s100b) + geom_text_auc(r.s100b))
  }
  expect_ggroc_doppelganger("geom_text_auc.screenshot", test_screenshot)
})

test_that("geom_text_auc does not use the ci.auc stored in the curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  roc_ci <- roc(aSAH$outcome, aSAH$s100b, ci = TRUE, quiet = TRUE)
  p <- ggroc(roc_ci) + geom_text_auc(roc_ci)
  expect_identical(ggplot2::layer_data(p, 2)$label, ggroc_auc_label(roc_ci$auc))
})

test_that("geom_text_auc gives a clear error on a curve without AUC", {
  r <- roc(aSAH$outcome, aSAH$s100b, auc = FALSE, quiet = TRUE)
  expect_error(geom_text_auc(r), "'data' has no 'auc'")
})

test_that("geom_text_auc(roc, ci = ) uses the given ci.auc", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  ci1 <- ci.auc(r.s100b)
  expect_identical(ggplot2::layer_data(ggroc(r.s100b) + geom_text_auc(r.s100b, ci = ci1), 2)$label, ggroc_auc_label(r.s100b$auc, ci1))
  # it takes precedence over the stored CI
  roc_ci <- roc(aSAH$outcome, aSAH$s100b, ci = TRUE, quiet = TRUE)
  ci90 <- ci.auc(roc_ci, conf.level = 0.9)
  expect_identical(ggplot2::layer_data(ggroc(roc_ci) + geom_text_auc(roc_ci, ci = ci90), 2)$label, ggroc_auc_label(roc_ci$auc, ci90))
})

test_that("geom_text_auc takes a ci.auc as text.ci.auc does", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  ci1 <- ci.auc(r.s100b)
  expect_identical(
    ggplot2::layer_data(ggroc(r.s100b) + geom_text_auc(ci1), 2)$label,
    ggroc_auc_label(r.s100b$auc, ci1)
  )
})

test_that("geom_text_auc stacks the labels of a list without facets", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)
  d <- ggplot2::layer_data(ggroc(rocs) + geom_text_auc(rocs, vjust = 0), 2)
  expect_equal(as.integer(d$PANEL), c(1L, 1L, 1L))
  expect_equal(d$vjust, c(0, 1.5, 3))
  expect_equal(d$label, vapply(rocs, function(r) ggroc_auc_label(r$auc), "", USE.NAMES = FALSE))
})

test_that("geom_text_auc puts each label of a list at (x, y) in its facet", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)
  # Facets added after the layer
  d <- ggplot2::layer_data(ggroc(rocs) + geom_text_auc(rocs, x = 0.4, y = 0.2) + ggplot2::facet_grid(. ~ name), 2)
  expect_equal(as.integer(d$PANEL), 1:3)
  expect_equal(d$vjust, c(1, 1, 1))
  expect_equal(d$y, c(0.2, 0.2, 0.2))
})

test_that("geom_text_auc labels of a list follow the elements", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  ci1 <- ci.auc(r.s100b)
  # A ci.auc, an auc and an unnamed roc
  l <- list(s100b = ci1, wfns = r.wfns$auc, r.ndka)
  labels <- c(ggroc_auc_label(r.s100b$auc, ci1), ggroc_auc_label(r.wfns$auc), ggroc_auc_label(r.ndka$auc))
  p <- ggroc(list(s100b = r.s100b, wfns = r.wfns, r.ndka))
  expect_equal(ggplot2::layer_data(p + geom_text_auc(l), 2)$label, labels)
  # name = TRUE starts them with the names, also with a pattern
  expect_equal(ggplot2::layer_data(p + geom_text_auc(l, name = TRUE), 2)$label, paste(c("s100b", "wfns", "3"), labels))
  d <- ggplot2::layer_data(p + geom_text_auc(list(r.s100b, r.wfns), pattern = "%.2f", name = TRUE), 2)
  expect_equal(d$label, paste(1:2, sprintf("%.2f", c(r.s100b$auc, r.wfns$auc))))
})

test_that("geom_text_auc labels take the colour of their curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)
  colour_of <- function(p, i) unique(ggplot2::layer_data(p, i)$colour)
  # The colour mapped to the curves, with the scale of the plot
  for (p in list(ggroc(rocs) + geom_text_auc(rocs),
                 ggroc(rocs) + geom_text_auc(rocs) + ggplot2::scale_colour_brewer(palette = "Dark2"))) {
    line_colours <- tapply(ggplot2::layer_data(p, 1)$colour, ggplot2::layer_data(p, 1)$group, unique)
    expect_equal(ggplot2::layer_data(p, 2)$colour, unname(as.vector(line_colours)))
  }
  # The text adds no glyph to the legend
  p <- ggroc(rocs) + geom_text_auc(rocs)
  expect_false(p@layers[[2]]$show.legend)
  # The fixed colour of the ROC line
  expect_equal(colour_of(ggroc(r.s100b, colour = "red") + geom_text_auc(r.s100b), 2), "red")
  expect_equal(colour_of(ggroc(rocs, aes = "linetype", colour = "red") + geom_text_auc(rocs), 2), "red")
  # Neither: the ggplot2 default
  expect_equal(colour_of(ggroc(r.s100b) + geom_text_auc(r.s100b), 2), "black")
  expect_equal(colour_of(ggroc(rocs, aes = "linetype") + geom_text_auc(rocs), 2), "black")
  # An explicit colour wins
  expect_equal(colour_of(ggroc(rocs) + geom_text_auc(rocs, colour = "grey30"), 2), "grey30")
})

test_that("geom_text_auc works with a list, with and without facets", {
  rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)
  test_geom_text_auc_list_screenshot <- function() {
    print(ggroc(rocs) + geom_text_auc(rocs))
  }
  expect_ggroc_doppelganger("geom_text_auc.list.screenshot", test_geom_text_auc_list_screenshot)
  test_geom_text_auc_list_facet_screenshot <- function() {
    print(ggroc(rocs) + geom_text_auc(rocs) + ggplot2::facet_grid(. ~ name))
  }
  expect_ggroc_doppelganger("geom_text_auc.list.facet.screenshot", test_geom_text_auc_list_facet_screenshot)
})

test_that("geom_text_auc accepts a named vjust when labels are not stacked", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  d <- ggplot2::layer_data(ggroc(r.s100b) + geom_text_auc(r.s100b, vjust = "top"), 2)
  expect_equal(d$vjust, "top")
  rocs <- list(s100b = r.s100b, wfns = r.wfns)
  d <- ggplot2::layer_data(ggroc(rocs) + ggplot2::facet_grid(. ~ name) + geom_text_auc(rocs, vjust = "top"), 2)
  expect_equal(d$vjust, c("top", "top"))
  # Stacking needs numbers
  expect_warning(ggplot2::layer_data(ggroc(rocs) + geom_text_auc(rocs, vjust = "top"), 2), "numeric 'vjust'")
})

test_that("geom_text_auc refuses a 'ci' shared by the curves of a list", {
  rocs <- list(s100b = r.s100b, wfns = r.wfns)
  expect_error(geom_text_auc(rocs, ci = ci.auc(r.s100b)), "list of 'ci.auc'")
})

test_that("geom_text_auc adds nothing with an empty list, as the other layers", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  expect_equal(layer_geom_classes(ggroc(r.s100b) + geom_text_auc(list())), "GeomLine")
  expect_equal(layer_geom_classes(ggroc(r.s100b) + geom_polygon_auc(list())), "GeomLine")
})

test_that("geom_text_auc puts the labels of a list in facets added before the layer", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  rocs <- list(s100b = r.s100b, wfns = r.wfns, ndka = smooth(r.ndka))
  d <- ggplot2::layer_data(ggroc(rocs) + ggplot2::facet_wrap(~name) + geom_text_auc(rocs), 2)
  expect_equal(as.integer(d$PANEL), 1:3)
  expect_equal(d$vjust, c(1, 1, 1))
})
