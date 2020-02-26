library(pROC)
data(aSAH)

context("smooth")

# Define some density functions

unif.density <- function(x, n, from, to, bw, kernel, ...) {
	smooth.x <- seq(from = from, to = to, length.out = n)
	smooth.y <- dunif(smooth.x, min = min(x), max = max(x))
	return(smooth.y)
}

norm.density <- function(x, n, from, to, bw, kernel, ...) {
	smooth.x <- seq(from = from, to = to, length.out = n)
	smooth.y <- dnorm(smooth.x, mean = mean(x), sd = sd(x))
	return(smooth.y)
}

lnorm.density <- function(x, n, from, to, bw, kernel, ...) {
	smooth.x <- seq(from = from, to = to, length.out = n)
	smooth.y <- dlnorm(smooth.x, meanlog = mean(x), sdlog = sd(x))
	return(smooth.y)
}

test_that("We fall back to the standard smooth", {
	tukey <- smooth(c(4, 1, 3, 6, 6, 4, 1, 6, 2, 4, 2))
	expect_is(tukey, "tukeysmooth")
	expect_equal(as.numeric(tukey), c(3, 3, 3, 3, 4, 4, 4, 4, 2, 2, 2))
})

test_that("smooth with a density function works", {
	smoothed <- smooth(r.ndka, method="density", density = unif.density, n = 10)
	expect_is(smoothed, "smooth.roc")
	expect_equal(smoothed$sensitivities, c(1, 1, 1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0, 0))
	expect_equal(smoothed$specificities, c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1))
	expect_equal(as.numeric(smoothed$auc), 0.9375)
})

test_that("smooth with two density functions works", {
	smoothed <- smooth(r.ndka, method="density", density.controls = norm.density, density.cases = lnorm.density, n = 10)
	expect_is(smoothed, "smooth.roc")
	expect_equal(smoothed$sensitivities, c(1, 1, 1, 0.635948942024884, 0.460070154191559, 0.344004532431686, 
										   0.25735248652959, 0.188201024566009, 0.130658598389315, 0.0813814489619488, 
										   0.0382893349015216, 0))
	expect_equal(smoothed$specificities, c(0, 0, 0.832138478872629, 0.99999996787709, 1, 1, 1, 1, 1, 1, 1, 1))
	expect_equal(as.numeric(smoothed$auc), 0.9694449)
})
