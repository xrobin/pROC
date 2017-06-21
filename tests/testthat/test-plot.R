context("plot")

# Tests powered by vdiffr. 
# To update the reference with vdiffr:
# > library(vdiffr)
# > manage_cases()

test_that("plot draws correctly", {
	test_basic_plot <- function() plot(r)
	# S100b
	r <- r.s100b
	vdiffr::expect_doppelganger("basic-s100b", test_basic_plot)
	
	r <- r.ndka
	vdiffr::expect_doppelganger("basic-ndka", test_basic_plot)
	
	r <- r.wfns
	vdiffr::expect_doppelganger("basic-wfns", test_basic_plot)
})

test_that("legacy.axis works correctly", {
	r <- r.s100b
	test_legacy.axis_plot <- function() plot(r, legacy.axes=TRUE)
	vdiffr::expect_doppelganger("legacy.axes", test_legacy.axis_plot)
})