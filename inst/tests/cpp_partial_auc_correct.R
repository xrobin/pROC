context("Partial auc correct (c++)")

test_that("Partial auc correct on null ROC works (c++)", {
	
	# Diagonal null ROC
	expect_that(pROC:::runit_partial_auc_correct(5e-05, 1, .99), equals(0.5))
	expect_that(pROC:::runit_partial_auc_correct(0.005, 1, .9), equals(0.5))
	expect_that(pROC:::runit_partial_auc_correct(0.02, 1, .8), equals(0.5))
	expect_that(pROC:::runit_partial_auc_correct(0.015, .9, .8), equals(0.5))
	expect_that(pROC:::runit_partial_auc_correct(0.24, .5, .1), equals(0.5))
	expect_that(pROC:::runit_partial_auc_correct(0.375, .5, .0), equals(0.5))
	expect_that(pROC:::runit_partial_auc_correct(0.08, .1, .0), equals(0.5))
})
		  
test_that("Partial auc correct on perfect ROC works (c++)", {
	
	# Diagonal null ROC
	expect_that(pROC:::runit_partial_auc_correct(.01, 1, .99), equals(1.0))
	expect_that(pROC:::runit_partial_auc_correct(0.1, 1, .9), equals(1.0))
	expect_that(pROC:::runit_partial_auc_correct(0.2, 1, .8), equals(1.0))
	expect_that(pROC:::runit_partial_auc_correct(0.1, .9, .8), equals(1.0))
	expect_that(pROC:::runit_partial_auc_correct(0.4, .5, .1), equals(1.0))
	expect_that(pROC:::runit_partial_auc_correct(0.5, .5, .0), equals(1.0))
	expect_that(pROC:::runit_partial_auc_correct(0.1, .1, .0), equals(1.0))
})

test_that("Partial auc correct on opposite ROC works (c++)", {
	
	# Diagonal null ROC
	expect_that(pROC:::runit_partial_auc_correct(0, 1, .99), equals(0))
	expect_that(pROC:::runit_partial_auc_correct(0, 1, .9), equals(0))
	expect_that(pROC:::runit_partial_auc_correct(0, 1, .8), equals(0))
	expect_that(pROC:::runit_partial_auc_correct(0, .9, .8), equals(0))
	expect_that(pROC:::runit_partial_auc_correct(0, .5, .1), equals(0))
	expect_that(pROC:::runit_partial_auc_correct(0, .5, .0), equals(0))
	expect_that(pROC:::runit_partial_auc_correct(0, .1, .0), equals(0))
})

test_that("Partial auc correct on WFNS ROC works (c++)", {
	
	# Diagonal null ROC
	expect_that(pROC:::runit_partial_auc_correct(0.000395121951219513, 1, .99), equals(0.517342811619071))
	expect_that(pROC:::runit_partial_auc_correct(0.0334417344173442, 1, .9), equals(0.649693339038653))
	expect_that(pROC:::runit_partial_auc_correct(0.0932791327913279, 1, .8), equals(0.703553146642578))
	expect_that(pROC:::runit_partial_auc_correct(0.0598373983739837, .9, .8), equals(0.763749402199904))
	expect_that(pROC:::runit_partial_auc_correct(0.38860909690178, .5, .1), equals(0.952537903757416))
	expect_that(pROC:::runit_partial_auc_correct(0.488134475939354, .5, .0), equals(0.952537903757416))
	expect_that(pROC:::runit_partial_auc_correct(0.0995253790375742, .1, .0), equals(0.952537903757416))
})

test_that("Partial auc correct on NDKA ROC works (c++)", {
	
	# Diagonal null ROC
	expect_that(pROC:::runit_partial_auc_correct(0.00024390243902439, 1, .99), equals(0.509743841157005))
	expect_that(pROC:::runit_partial_auc_correct(0.0107046070460705, 1, .9), equals(0.530024247610897))
	expect_that(pROC:::runit_partial_auc_correct(0.0384823848238482, 1, .8), equals(0.551339957844023))
	expect_that(pROC:::runit_partial_auc_correct(0.0277777777777778, .9, .8), equals(0.575163398692811))
	expect_that(pROC:::runit_partial_auc_correct(0.32320460704607, .5, .1), equals(0.680019196025294))
	expect_that(pROC:::runit_partial_auc_correct(0.416836043360434, .5, .0), equals(0.667344173441734))
	expect_that(pROC:::runit_partial_auc_correct(0.0936314363143631, .1, .0), equals(0.363143631436315))
})

test_that("Partial auc correct on S100B ROC works (c++)", {
	
	# Diagonal null ROC
	expect_that(pROC:::runit_partial_auc_correct(0.00292682926829269, 1, .99), equals(0.644564284838828))
	expect_that(pROC:::runit_partial_auc_correct(0.0327574525745257, 1, .9), equals(0.646091855655399))
	expect_that(pROC:::runit_partial_auc_correct(0.0805894308943089, 1, .8), equals(0.668303974706414))
	expect_that(pROC:::runit_partial_auc_correct(0.0478319783197832, .9, .8), equals(0.693129284234019))
	expect_that(pROC:::runit_partial_auc_correct(0.350567411924119, .5, .1), equals(0.794030883017164))
	expect_that(pROC:::runit_partial_auc_correct(0.448128387533875, .5, .0), equals(0.792513550135501))
	expect_that(pROC:::runit_partial_auc_correct(0.0975609756097561, .1, .0), equals(0.75609756097561))
})
