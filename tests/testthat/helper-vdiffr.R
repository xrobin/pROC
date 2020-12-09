
# Skip expect_doppelganger if vdiffr is not installed
expect_doppelganger <- function(title, fig, path = NULL, ...) {
	testthat::skip_if_not_installed("vdiffr")
	vdiffr::expect_doppelganger(title, fig, path = path, ...)
}