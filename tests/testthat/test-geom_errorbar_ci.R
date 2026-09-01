context("geom_errorbar_ci")

test_that("geom_errorbar_ci se bars are drawn in front of the ROC curve", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  set.seed(42)
  ci_se <- ci.se(r.s100b, specificities = seq(0, 1, 0.1), boot.n = 20)
  p <- ggroc(r.s100b) + geom_errorbar_ci(ci_se)
  expect_equal(layer_geom_classes(p), c("GeomLine", "GeomErrorbar"))
})

test_that("geom_errorbar_ci thresholds draws a cross in front of the ROC", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  set.seed(42)
  ci_th <- ci.thresholds(r.s100b, thresholds = "best", boot.n = 20)
  p <- ggroc(r.s100b) + geom_errorbar_ci(ci_th)
  expect_equal(layer_geom_classes(p), c("GeomLine", "GeomErrorbar", "GeomErrorbar"))
})

test_that("geom_errorbar_ci sp screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    set.seed(42)
    ci_sp <- ci.sp(r.s100b, sensitivities = seq(0, 1, 0.1), boot.n = 20)
    print(ggroc(r.s100b) + geom_errorbar_ci(ci_sp))
  }
  expect_ggroc_doppelganger("geom_errorbar_ci.sp.screenshot", test_screenshot)
})

test_that("geom_errorbar_ci thresholds screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    set.seed(42)
    ci_th <- ci.thresholds(r.s100b, thresholds = "best", boot.n = 20)
    print(
      ggroc(r.s100b) +
        geom_errorbar_ci(ci_th) +
        geom_roc_threshold(r.s100b, thresholds = "best")
    )
  }
  expect_ggroc_doppelganger("geom_errorbar_ci.thresholds.screenshot", test_screenshot)
})
