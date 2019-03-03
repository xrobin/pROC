context("ggroc")



test_that("Ggroc screenshot looks normal", {
	test_ggplot_screenshot <- function() {
		print(ggroc(r.s100b.percent, , alpha = 0.5, colour = "red", linetype = 2, size = 2))
	}
	vdiffr::expect_doppelganger("ggroc.screenshot", test_ggplot_screenshot)
})

test_that("Ggroc list screenshot looks normal", {
	test_ggplot_list_screenshot <- function() {
		print(ggroc(list(s100b=r.s100b, wfns=r.wfns, ndka=r.ndka)))
	}
	vdiffr::expect_doppelganger("ggroc.list.screenshot", test_ggplot_list_screenshot)
})

test_that("Ggroc list can take multiple aes", {
	test_ggplot_list_screenshot <- function() {
		print(ggroc(list(s100b=r.s100b, wfns=r.wfns, ndka=r.ndka), aes=c("c", "linetype", "size")))
	}
	vdiffr::expect_doppelganger("ggroc.list.multi.aes", test_ggplot_list_screenshot)
})

test_that("Ggroc list extra aestetics screenshot looks normal", {
	test_ggplot_list_extra_aes_screenshot <- function() {
		print(ggroc(list(s100b=r.s100b, wfns=r.wfns, ndka=r.ndka), aes="linetype", color="red"))
	}
	vdiffr::expect_doppelganger("ggroc.list.extra.aes.screenshot", test_ggplot_list_extra_aes_screenshot)
})

test_that("Ggroc list with group facet screenshot looks normal", {
	test_ggplot_list_group_facet_screenshot <- function() {
		library(ggplot2)
		g <- ggroc(list(s100b=r.s100b, wfns=r.wfns, ndka=r.ndka), aes="group") + facet_grid(.~name)
		print(g)
	}
	vdiffr::expect_doppelganger("ggroc.list.group.facet.screenshot", test_ggplot_list_group_facet_screenshot)
})
