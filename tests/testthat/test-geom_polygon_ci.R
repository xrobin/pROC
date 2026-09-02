context("geom_polygon_ci")

test_that("geom_polygon_ci se is inserted behind the ROC curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  set.seed(42)
  ci_se <- ci.se(r.s100b, specificities = seq(0, 1, 0.05), boot.n = 20)
  p <- ggroc(r.s100b) + geom_polygon_ci(ci_se)
  expect_equal(layer_geom_classes(p), c("GeomPolygon", "GeomLine"))
})

test_that("geom_polygon_ci se screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    set.seed(42)
    ci_se <- ci.se(r.s100b, specificities = seq(0, 1, 0.05), boot.n = 20)
    print(ggroc(r.s100b) + geom_polygon_ci(ci_se, fill = "#1c61b6", alpha = 0.5))
  }
  expect_ggroc_doppelganger("geom_polygon_ci.se.screenshot", test_screenshot)
})

test_that("geom_polygon_ci sp is inserted behind the ROC curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  set.seed(42)
  ci_sp <- ci.sp(r.s100b, sensitivities = seq(0, 1, 0.05), boot.n = 20)
  p <- ggroc(r.s100b) + geom_polygon_ci(ci_sp)
  expect_equal(layer_geom_classes(p), c("GeomPolygon", "GeomLine"))
})

test_that("geom_polygon_ci sp screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    set.seed(42)
    ci_sp <- ci.sp(r.s100b, sensitivities = seq(0, 1, 0.05), boot.n = 20)
    print(ggroc(r.s100b) + geom_polygon_ci(ci_sp, fill = "#1c61b6", alpha = 0.5))
  }
  expect_ggroc_doppelganger("geom_polygon_ci.sp.screenshot", test_screenshot)
})
