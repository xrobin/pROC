library(pROC)
data(aSAH)


test_that("var with delong works", {
	expect_equal(var(r.wfns), 0.00146991470882363)
	expect_equal(var(r.ndka), 0.0031908105493913)
	expect_equal(var(r.s100b), 0.00266868245717244)
})


test_that("var works with auc", {
	expect_equal(var(auc(r.wfns)), 0.00146991470882363)
	expect_equal(var(auc(r.ndka)), 0.0031908105493913)
	expect_equal(var(auc(r.s100b)), 0.00266868245717244)
})
