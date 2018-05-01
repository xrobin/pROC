library(pROC)
data(aSAH)

numacc.response <- c(2, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2)
numacc.predictor <- c(0.960602681556147, 0.0794407386056549, 0.144842404246611, 
			   0.931816485855784, 0.931816485855784, 0.97764041048215, 0.653549466997938699464, 
			   0.796401132206396, 0.427720540184519, 0.811278021288732, 0.0188323116581187, 
			   0.653549466997938588442, 0.653549466997938477419, 0.959111701445925, 0.931816485855784, 
			   0.663663279418747, 0.800100838413179, 0.780456095511079)
# Predictor has near-ties that can break numerical comparisons

test_that("AUC is consistent across algorithms with numerical near-ties", {
	r1 <- roc(numacc.response, numacc.predictor, algorithm=1)
	r2 <- roc(numacc.response, numacc.predictor, algorithm=2)
	r3 <- roc(numacc.response, numacc.predictor, algorithm=3)
	expect_equal(as.numeric(auc(r1)), as.numeric(auc(r2)))
	expect_equal(as.numeric(auc(r1)), as.numeric(auc(r3)))
})

test_that("AUC is consistent across algorithms with numerical near-ties and direction = >", {
	r1 <- roc(2-numacc.response, numacc.predictor, algorithm=1)
	r2 <- roc(2-numacc.response, numacc.predictor, algorithm=2)
	r3 <- roc(2-numacc.response, numacc.predictor, algorithm=3)
	expect_equal(as.numeric(auc(r1)), as.numeric(auc(r2)))
	expect_equal(as.numeric(auc(r1)), as.numeric(auc(r3)))
})

test_that("delong theta is consistent with auc", {
	r1 <- roc(numacc.response, numacc.predictor, algorithm=1)
	r2 <- roc(numacc.response, numacc.predictor, algorithm=2)
	r3 <- roc(numacc.response, numacc.predictor, algorithm=3)
	expect_equal(pROC:::delongPlacements(r1)$theta, as.numeric(auc(r1)))
	expect_equal(pROC:::delongPlacements(r2)$theta, as.numeric(auc(r2)))
	expect_equal(pROC:::delongPlacements(r3)$theta, as.numeric(auc(r3)))
})

test_that("delong theta is consistent with auc and direction = >", {
	r1 <- roc(2-numacc.response, numacc.predictor, algorithm=1)
	r2 <- roc(2-numacc.response, numacc.predictor, algorithm=2)
	r3 <- roc(2-numacc.response, numacc.predictor, algorithm=3)
	expect_equal(pROC:::delongPlacements(r1)$theta, as.numeric(auc(r1)))
	expect_equal(pROC:::delongPlacements(r2)$theta, as.numeric(auc(r2)))
	expect_equal(pROC:::delongPlacements(r3)$theta, as.numeric(auc(r3)))
})
