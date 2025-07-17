context("ggroc")



test_that("Ggroc screenshot looks normal", {
  skip_if(packageVersion("ggplot2") < "2.4")
  test_ggplot_screenshot <- function() {
    print(ggroc(r.s100b.percent, alpha = 0.5, colour = "red", linetype = 2, linewidth = 2))
  }
  expect_ggroc_doppelganger("ggroc.screenshot", test_ggplot_screenshot)
})

test_that("Ggroc works with legacy.axis and percent", {
  skip_if(packageVersion("ggplot2") < "2.4")

  # nothing
  test_ggplot_screenshot <- function() {
    print(ggroc(r.s100b))
  }
  expect_ggroc_doppelganger("ggroc.screenshot_base", test_ggplot_screenshot)

  # percent
  test_ggplot_screenshot <- function() {
    print(ggroc(r.s100b.percent))
  }
  expect_ggroc_doppelganger("ggroc.screenshot_percent", test_ggplot_screenshot)

  # legacy.axes
  test_ggplot_screenshot <- function() {
    print(ggroc(r.s100b, legacy.axes = TRUE))
  }
  expect_ggroc_doppelganger("ggroc.screenshot_legacy", test_ggplot_screenshot)

  # percent, legacy.axes
  test_ggplot_screenshot <- function() {
    print(ggroc(r.s100b.percent, legacy.axes = TRUE))
  }
  expect_ggroc_doppelganger("ggroc.screenshot_percent_legacy", test_ggplot_screenshot)
})

test_that("Ggroc list screenshot looks normal", {
  test_ggplot_list_screenshot <- function() {
    print(ggroc(list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka)))
  }
  expect_ggroc_doppelganger("ggroc.list.screenshot", test_ggplot_list_screenshot)
})

test_that("Ggroc list can take multiple aes", {
  test_ggplot_list_screenshot <- function() {
    print(ggroc(list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka), aes = c("c", "linetype")))
  }
  expect_ggroc_doppelganger("ggroc.list.multi.aes", test_ggplot_list_screenshot)
})

test_that("Ggroc list doesn't get merged with set colour", {
  test_ggplot_list_screenshot <- function() {
    print(ggroc(list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka), colour = "red"))
  }
  expect_ggroc_doppelganger("ggroc.list.colour", test_ggplot_list_screenshot)
})

test_that("Ggroc list extra aestetics screenshot looks normal", {
  test_ggplot_list_extra_aes_screenshot <- function() {
    print(ggroc(list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka), aes = "linetype", color = "red"))
  }
  expect_ggroc_doppelganger("ggroc.list.extra.aes.screenshot", test_ggplot_list_extra_aes_screenshot)
})

test_that("Ggroc list with group facet screenshot looks normal", {
  test_ggplot_list_group_facet_screenshot <- function() {
    library(ggplot2)
    g <- ggroc(list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka), aes = "group") + facet_grid(. ~ name)
    print(g)
  }
  expect_ggroc_doppelganger("ggroc.list.group.facet.screenshot", test_ggplot_list_group_facet_screenshot)
})

test_that("Ggroc aesthetics can be modified with scale_colour_manual", {
  test_ggplot_list_screenshot <- function() {
    print(ggroc(list(s100b = r.s100b, wfns = r.wfns, ndka = r.ndka), aes = c("c", "linetype")) +
      scale_colour_manual(values = c("purple", "yellow", "purple")))
  }
  expect_ggroc_doppelganger("ggroc.list.scale.colour.manual", test_ggplot_list_screenshot)
})


test_that("Ggroc screenshot looks normal with a single smooth.roc", {
  skip_if(packageVersion("ggplot2") < "2.4")
  test_ggplot_screenshot <- function() {
    print(ggroc(smooth(r.s100b), , alpha = 0.5, colour = "red", linetype = 2, linewidth = 2))
  }
  expect_ggroc_doppelganger("ggroc.smooth.screenshot", test_ggplot_screenshot)
})

test_that("Ggroc screenshot looks normal with a list of smooth.roc", {
  test_ggplot_screenshot <- function() {
    print(ggroc(list(s100b = smooth(r.s100b), wfns = smooth(r.wfns), ndka = smooth(r.ndka))))
  }
  expect_ggroc_doppelganger("ggroc.smooth.list.screenshot", test_ggplot_screenshot)
})
