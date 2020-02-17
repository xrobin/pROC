library(pROC)
data(aSAH)

context("ci.formula")

test_that("bootstrap cov works with smooth and !reuse.auc", {
	skip_slow()
	if (R.version$minor >= "6.0") {
		RNGkind(sample.kind="Rounding")
	}
	
	for (pair in list(
			list(ci, list()),
			list(ci.se, list(boot.n = 10)),
			list(ci.sp, list(boot.n = 10)),
			list(ci.thresholds, list(boot.n = 10)),
			list(ci.coords, list(boot.n = 10, x = 0.5)),
			list(ci.auc, list()))) {
		fun <- pair[[1]]

		# First calculate ci with .default
		args.default <- c(
			list(response = aSAH$outcome,
				 predictor = aSAH$s100b),
			pair[[2]])
		set.seed(42) # For reproducible CI
		obs.default <- do.call(fun, args.default)

		# Then with .formula
		args.formula <- c(
			list(formula = outcome ~ s100b,
				 data = aSAH),
			pair[[2]])
		set.seed(42) # For reproducible CI
		obs.formula <- do.call(fun, args.formula)

		# Here we check both returned the same result
		# We ignore attributes, as we have different
		# roc objects, and unfortunately equivalent means
		# we only test near equality
		expect_equivalent(obs.default, obs.formula)
	}
})
