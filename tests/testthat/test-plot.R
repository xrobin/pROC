context("plot")

# Tests powered by vdiffr. 
# To update the reference with vdiffr:
# > library(vdiffr)
# > manage_cases()

test_that("plot draws correctly", {
	r <- roc(aSAH$outcome, aSAH$s100b)
	f <- function() plot(r)
	vdiffr::expect_doppelganger("basic", f)
})