library(pROC)
data(aSAH)


test_that("can get bar works", {
	my_progress <- pROC:::roc_utils_get_progress_bar("text")
	expect_equal(names(my_progress), c("init", "step", "term"))
})

test_that("progress = none works", {
	expect_silent(var(r.wfns, method="bootstrap", boot.n=2, progress="none"))
})
