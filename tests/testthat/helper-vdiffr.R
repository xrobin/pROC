
# Skip expect_doppelganger if vdiffr is not installed
expect_doppelganger <- function(title, fig, ...) {
	testthat::skip_if_not_installed("vdiffr")
	vdiffr::expect_doppelganger(title, fig, ...)
}

# expect_doppelganger for ggroc
expect_ggroc_doppelganger <- function(title, fig, ...) {
	testthat::skip_if_not_installed("ggplot2")
	expect_doppelganger(title, fig, ...)
}
