context("geom_text_auc")

test_that("geom_text_auc screenshot looks normal", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  test_screenshot <- function() {
    print(ggroc(r.s100b) + geom_text_auc(r.s100b))
  }
  expect_ggroc_doppelganger("geom_text_auc.screenshot", test_screenshot)
})

test_that("geom_text_auc includes ci.auc when present", {
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
  roc_ci <- roc(aSAH$outcome, aSAH$s100b, ci = TRUE, quiet = TRUE)
  p <- ggroc(roc_ci) + geom_text_auc(roc_ci)
  expect_true("GeomText" %in% layer_geom_classes(p))
})
