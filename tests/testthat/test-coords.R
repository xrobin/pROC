library(pROC)
data(aSAH)

context("coords")

test_that("coords with thresholds works", {
	return.rows <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "youden", "closest.topleft")
	obtained <- coords(r.s100b, "all", ret = return.rows, transpose=TRUE)
	expect_equal(obtained, expected.coords[return.rows,])
})

test_that("coords returns all thresholds by default", {
	obtained <- coords(r.s100b, transpose=TRUE)
	expect_equal(obtained, expected.coords[c("threshold", "specificity", "sensitivity"),])
	# but not if it's an empty numeric, as this might be indicative of user error
	expect_error(coords(r.s100b, numeric(0)), "length")
})


test_that("coords returns all thresholds by default with smooth.roc", {
	obtained <- coords(smooth(r.s100b))
	expect_equal(obtained, expected.coords.smooth[,c("specificity", "sensitivity")])
	# but not if it's an empty numeric, as this might be indicative of user error
	expect_error(coords(r.s100b, numeric(0)), "length")
})


test_that("coords returns all columns with ret = 'all' with smooth.roc", {
	obtained <- coords(smooth(r.s100b), ret = "all")
	expect_equal(obtained, expected.coords.smooth)
})


test_that("coords with transpose = FALSE works", {
	return.rows <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "youden", "closest.topleft")
	obtained <- coords(r.s100b, "all", ret = return.rows, transpose = FALSE)
	expect_equal(obtained, as.data.frame(t(expected.coords[return.rows,])))
	obtained <- coords(r.s100b, transpose = FALSE)
	expect_equal(obtained, as.data.frame(t(expected.coords[c("threshold", "specificity", "sensitivity"),])))
	
	# With drop=TRUE
	obtained <- coords(r.s100b, "all", ret = "se", transpose = FALSE, drop=TRUE)
	expect_is(obtained, "numeric")
	#  Not why drop.data.frame returns a list, skipping
	# obtained <- coords(r.s100b, "best", ret = "all", transpose = FALSE, drop=TRUE)
	
	# With drop=FALSE
	obtained <- coords(r.s100b, "all", ret = "se", transpose = FALSE, drop=FALSE)
	expect_is(obtained, "data.frame")
})


test_that("coords with ret='all' works", {
	obtained <- coords(r.s100b, "all", ret = "all", transpose=TRUE)
	expect_equal(dim(obtained), c(24, 51))
	expect_equal(obtained[rownames(expected.coords),], expected.coords)
})


test_that("coords with ret='all' doesn't accept additional options", {
	expect_error(coords(r.s100b, "all", ret = c("all", "thresholds")))
})


test_that("coords with percent works", {
	return.rows <- "all"
	percent.rows <- c("specificity", "sensitivity", "accuracy", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "youden", "closest.topleft", "fdr", "fpr", "tpr", "tnr", "fnr", "precision", "recall")
	obtained.percent <- coords(r.s100b.percent, "all", ret = return.rows, transpose=TRUE)
	# Adjust for percent
	obtained.percent[percent.rows,] <- obtained.percent[percent.rows,] / 100
	expect_equal(obtained.percent, expected.coords)
})


test_that("coords with local maximas thresholds works", {
	return.rows <- "all"
	obtained <- coords(r.s100b, "local maximas", ret = return.rows, transpose=TRUE)
	expected.thresholds <- c(-Inf, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.135, 0.155, 0.205, 0.245, 0.29, 0.325, 0.345, 0.395, 0.435, 0.475, 0.485, 0.51)
	expect_equal(as.vector(obtained["threshold",]), expected.thresholds)
	expect_equivalent(obtained, expected.coords[,expected.coords["threshold",] %in% expected.thresholds])
})


test_that("coords with best threshold works", {
	return.rows <- "all"
	obtained <- coords(r.s100b, "best", ret = return.rows, transpose=TRUE)
	expect_equal(obtained, expected.coords[,expected.coords["threshold",] == 0.205])
})


test_that("coords with arbitrary thresholds works", {
	return.rows <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "youden", "closest.topleft")
	obtained <- coords(r.s100b, c(0.205, 0.055), input = "threshold", ret = return.rows, transpose=TRUE)
	expect_equivalent(obtained, expected.coords[return.rows, c(18, 4)])
})

test_that("coords with arbitrary thresholds at exact data point works", {
	return.rows <- "all"
	expect_equal(sum(aSAH$s100b == 0.05),  3)
	expect_equal(sum(aSAH$s100b == 0.52),  1)
	obtained <- coords(r.s100b, c(0.05, 0.52), input = "threshold", ret = return.rows, transpose=TRUE)
	expect_equivalent(obtained[-1,], expected.coords[-1, c(3, 40)])
})

test_that("coords with arbitrary thresholds works with direction=>", {
	obtained <- coords(r.100b.reversed, c(0.05, 0.055, 0.205, 0.52), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"), transpose=TRUE)
	expect_equivalent(obtained, expected.coords.reverse)
})


test_that("coords with single arbitrary threshold works", {
	return.rows <- "all"
	obtained <- coords(r.s100b, c(0.205), input = "threshold", ret = return.rows, transpose=TRUE)
	expect_equal(obtained, expected.coords[, c(18), drop=T])
})


test_that("coords with arbitrary thresholds at exact data point works", {
	expect_equal(sum(aSAH$s100b == 0.05),  3)
	expect_equal(sum(aSAH$s100b == 0.52),  1)
	obtained <- coords(r.s100b, c(0.05), input = "threshold", ret = "all", transpose=TRUE)
	expect_equal(obtained[-1], expected.coords[-1, 3])
	obtained <- coords(r.s100b, c(0.52), input = "threshold", ret = "all", transpose=TRUE)
	expect_equal(obtained[-1], expected.coords[-1, 40])
})


test_that("coords with arbitrary thresholds works with direction=>", {
	obtained <- coords(r.100b.reversed, c(0.05), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"), transpose=TRUE)
	expect_equal(obtained, expected.coords.reverse[, 1])
	obtained <- coords(r.100b.reversed, c(0.055), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"), transpose=TRUE)
	expect_equal(obtained, expected.coords.reverse[, 2])
	obtained <- coords(r.100b.reversed, c(0.205), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"), transpose=TRUE)
	expect_equal(obtained, expected.coords.reverse[, 3])
	obtained <- coords(r.100b.reversed, c(0.52), input = "threshold", ret = c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv"), transpose=TRUE)
	expect_equal(obtained, expected.coords.reverse[, 4])
})


test_that("coords with sensitivity works", {
	obtained <- coords(r.s100b, seq(0, 1, .1), input = "sensitivity", ret = c("threshold", "specificity", "sensitivity"), transpose=TRUE)
	expect_equal(unname(obtained["threshold",]), c(Inf, rep(NA, 9), -Inf))
	expect_equal(unname(obtained["sensitivity",]), seq(0, 1, .1))
	expect_equal(unname(obtained["specificity",]), c(1, 1, 1, 0.972222222222222, 0.888888888888889, 0.833333333333333, 0.805555555555556, 0.56875, 0.447222222222222, 0.230555555555556, 0))
})


test_that("coords with sensitivity works with percent", {
	obtained <- coords(r.s100b.percent, seq(0, 100, 10), input = "sensitivity", ret = c("threshold", "specificity", "sensitivity"), transpose=TRUE)
	expect_equal(unname(obtained["threshold",]), c(Inf, rep(NA, 9), -Inf))
	expect_equal(unname(obtained["sensitivity",]), seq(0, 100, 10))
	expect_equal(unname(obtained["specificity",]), c(1, 1, 1, 0.972222222222222, 0.888888888888889, 0.833333333333333, 0.805555555555556, 0.56875, 0.447222222222222, 0.230555555555556, 0) * 100)
})


test_that("coords with specificity works", {
	obtained <- coords(r.s100b, seq(0, 1, .1), input = "specificity", ret = c("threshold", "specificity", "sensitivity"), transpose=TRUE)
	expect_equal(unname(obtained["threshold",]), c(-Inf, rep(NA, 9), 0.51))
	expect_equal(unname(obtained["specificity",]), seq(0, 1, .1))
	expect_equal(unname(obtained["sensitivity",]), c(1, 0.975609756097561, 0.921951219512195, 0.879674796747967, 0.823693379790941, 0.774390243902439, 0.675609756097561, 0.655284552845528, 0.634146341463415, 0.390243902439024, 0.292682926829268))
})


test_that("coords with specificity works with percent", {
	obtained <- coords(r.s100b.percent, seq(0, 100, 10), input = "specificity", ret = c("threshold", "specificity", "sensitivity"), transpose=TRUE)
	expect_equal(unname(obtained["threshold",]), c(-Inf, rep(NA, 9), 0.51))
	expect_equal(unname(obtained["specificity",]), seq(0, 100, 10))
	expect_equal(unname(obtained["sensitivity",]), c(1, 0.975609756097561, 0.921951219512195, 0.879674796747967, 0.823693379790941, 0.774390243902439, 0.675609756097561, 0.655284552845528, 0.634146341463415, 0.390243902439024, 0.292682926829268) * 100)
})


test_that("coords with specificity works with as.list", {
	obtained <- coords(r.s100b.percent, "best", ret = c("threshold", "specificity", "accuracy"), as.list = TRUE)
	expect_equal(obtained, list(
		threshold = 0.205,
		specificity = unname(expected.coords["specificity", 18]) * 100,
		accuracy = unname(expected.coords["accuracy", 18]) * 100
	))
})

test_that("coords with specificity works with as.list and drop=FALSE", {
	obtained <- coords(r.s100b.percent, "best", 
					   ret = c("threshold", "specificity", "accuracy"), 
					   as.list = TRUE, drop = FALSE)
	expect_equal(obtained[[1]], list(
		threshold = 0.205,
		specificity = unname(expected.coords["specificity", 18]) * 100,
		accuracy = unname(expected.coords["accuracy", 18]) * 100
	))
})


test_that("coords with specificity works with as.list and several thresholds", {
	obtained <- coords(r.s100b.percent, c(0.205, 0.51), 
					   ret = c("threshold", "specificity", "accuracy"), 
					   as.list = TRUE, drop = FALSE)
	expect_equal(obtained[[1]], list(
		threshold = 0.205,
		specificity = unname(expected.coords["specificity", 18]) * 100,
		accuracy = unname(expected.coords["accuracy", 18]) * 100
	))
	expect_equal(obtained[[2]], list(
		threshold = 0.51,
		specificity = unname(expected.coords["specificity", 40]) * 100,
		accuracy = unname(expected.coords["accuracy", 40]) * 100
	))
})


test_that("drop works", {
	# First make sure we get matrices with drop = FALSE
	expect_is(coords(r.s100b, 0.51, input = "threshold", ret = c("sensitivity", "specificity"), drop = FALSE, transpose=TRUE), "matrix")
	expect_is(coords(r.s100b, 0.51, input = "threshold", ret = "specificity", drop = FALSE, transpose=TRUE), "matrix")
	expect_is(coords(r.s100b, "local maximas", input = "threshold", ret = "specificity", drop = FALSE, transpose=TRUE), "matrix")	
	expect_is(coords(r.s100b, c(0.51, 0.2), input = "threshold", ret = "specificity", drop = FALSE, transpose=TRUE), "matrix")
	# Look for numeric
	expect_is(coords(r.s100b, 0.51, input = "threshold", ret = c("sensitivity", "specificity"), drop = TRUE, transpose=TRUE), "numeric")
	expect_is(coords(r.s100b, 0.51, input = "threshold", ret = "specificity", drop = TRUE, transpose=TRUE), "numeric")
	expect_is(coords(r.s100b, "local maximas", input = "threshold", ret = "specificity", drop = TRUE, transpose=TRUE), "numeric")	
	expect_is(coords(r.s100b, c(0.51, 0.2), input = "threshold", ret = "specificity", drop = TRUE, transpose=TRUE), "numeric")
})


test_that("as.matrix works", {
	obtained <- coords(r.s100b, c(0.51, 0.205), ret="sensitivity", transpose = FALSE, as.matrix = TRUE)
	expect_equal(obtained, t(expected.coords["sensitivity", c(40, 18), drop = FALSE]))
})


test_that("as.matrix works with drop=TRUE", {
	obtained <- coords(r.s100b, c(0.51, 0.205), ret="sensitivity", transpose = FALSE, as.matrix = TRUE, drop = TRUE)
	expect_equal(obtained, expected.coords["sensitivity", c(40, 18), drop = TRUE])
})


test_that("coords returns the correct basic values ", {
	obtained <- coords(r.s100b, 0.205, 
					   ret = c("t", "tp", "fp", "tn", "fn",
					   		   "sp", "se", "acc",
					   		   "npv", "ppv", "precision", "recall",
					   		   "tpr", "fpr", "tnr", "fnr", "fdr"),
					   transpose=TRUE)
	
	obtained.percent <- coords(r.s100b.percent, 0.205, 
					   ret = c("t", "tp", "fp", "tn", "fn",
					   		"sp", "se", "acc",
					   		"npv", "ppv", "precision", "recall",
					   		"tpr", "fpr", "tnr", "fnr", "fdr"), 
					   transpose=TRUE)
	
	# We assume the following values:
	# tp fp tn fn N
	# 26 14 58 15 113
	
	expected <- c(
		threshold = 0.205,
		tp = 26,
		fp = 14,
		tn = 58,
		fn = 15,
		specificity = 58 / (58 + 14),
		sensitivity = 26 / (26 + 15),
		accuracy = (26 + 58) / 113,
		npv = 58 / (58 + 15),
		ppv = 26 / (26 + 14),
		precision = 26 / (26 + 14),
		recall = 26 / (26 + 15),
		tpr = 26 / (26 + 15),
		fpr = 1 - (58 / (58 + 14)),
		tnr = 58 / (58 + 14),
		fnr = 1 - (26 / (26 + 15)),
		fdr = 14 / (26 + 14)
	)
	
	expect_equal(obtained, expected)
	expect_equal(obtained.percent[1:5], expected[1:5])
	expect_equal(obtained.percent[6:17], expected[6:17]*100)
})


test_that("coords works with smooth.roc and x = 'best'", {
	smooth.s100b <- smooth(r.s100b)
	expect <- structure(c(0.750857175922901, 0.608610567514677, 0.699245574642041, 
						  54.0617166664488, 24.9530332681018, 16.0469667318982, 17.9382833335512, 
						  0.771112992655678, 0.581773544045047, 0.418226455954953, 0.249142824077099, 
						  0.608610567514677, 0.750857175922901, 0.391389432485323, 0.249142824077099, 
						  0.391389432485323, 0.300754425357959, 0.228887007344322, 0.418226455954953, 
						  0.581773544045047, 0.608610567514677, 1.35946774343758, 0.215257834650296
	), .Dim = c(23L, 1L), .Dimnames = list(c("specificity", "sensitivity", 
							"accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "fdr", "fpr", 
							"tpr", "tnr", "fnr", "1-specificity", "1-sensitivity", "1-accuracy", 
							"1-npv", "1-ppv", "precision", "recall", "youden", "closest.topleft"
	), NULL))
	
	
	reduced.cols <- c("specificity", "sensitivity", "youden")
	
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, transpose=TRUE)
	expect_equal(obtained, expect[reduced.cols,])
	
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, drop = FALSE, transpose=TRUE)
	expect_equal(obtained, expect[reduced.cols,, drop=FALSE])
	
	obtained <- coords(smooth.s100b, "best", ret = "all", drop = FALSE, transpose=TRUE)
	expect_equal(obtained, expect)
	
	obtained <- coords(smooth.s100b, "best", ret = "all", transpose=TRUE)
	expect_equal(obtained, expect[, 1])
	
	obtained <- coords(smooth.s100b, "best", ret = "all", drop = FALSE, transpose=TRUE)
	expect_equal(obtained, expect)
	
	obtained <- coords(smooth.s100b, "best", ret = "all", as.list = TRUE)
	expect_equal(obtained, as.list(expect[, 1]))
	expect_equal(names(obtained), rownames(expect))
	
	obtained <- coords(smooth.s100b, "best", ret = "all", as.list = TRUE, drop = FALSE)
	expect_equal(obtained[[1]], as.list(expect[, 1])) # names
	expect_equal(names(obtained[[1]]), rownames(expect))
	
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, as.list = TRUE)
	expect_equal(obtained, as.list(expect[reduced.cols, 1]))
	expect_equal(names(obtained), reduced.cols)
	
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, as.list = TRUE, drop = FALSE)
	expect_equal(obtained[[1]], as.list(expect[reduced.cols, 1])) # names
	expect_equal(names(obtained[[1]]), reduced.cols)
})


test_that("coords works with smooth.roc and transpose = FALSE", {
	smooth.s100b <- smooth(r.s100b)
	expect <- structure(c(0.750857175922901, 0.608610567514677, 0.699245574642041, 
						  54.0617166664488, 24.9530332681018, 16.0469667318982, 17.9382833335512, 
						  0.771112992655678, 0.581773544045047, 0.418226455954953, 0.249142824077099, 
						  0.608610567514677, 0.750857175922901, 0.391389432485323, 0.249142824077099, 
						  0.391389432485323, 0.300754425357959, 0.228887007344322, 0.418226455954953, 
						  0.581773544045047, 0.608610567514677, 1.35946774343758, 0.215257834650296
	), .Dim = c(23L, 1L), .Dimnames = list(c("specificity", "sensitivity", 
											 "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "fdr", "fpr", 
											 "tpr", "tnr", "fnr", "1-specificity", "1-sensitivity", "1-accuracy", 
											 "1-npv", "1-ppv", "precision", "recall", "youden", "closest.topleft"
	), NULL))
	
	
	reduced.cols <- c("specificity", "sensitivity", "youden")
	
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, drop = FALSE, transpose = FALSE)
	expect_equal(obtained, as.data.frame(t(expect[reduced.cols,, drop=FALSE])))
	
	obtained <- coords(smooth.s100b, "best", ret = "all", drop = FALSE, transpose = FALSE)
	expect_equal(obtained, as.data.frame(t(expect)))
	
	# Without drop
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, transpose = FALSE)
	expect_equivalent(obtained, as.data.frame(t(expect[reduced.cols,])))
	
	# drop = TRUE
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, drop = TRUE, transpose = FALSE)
	expect_equal(obtained, as.list(expect[reduced.cols,]))
	
	# With as.matrix
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, transpose = FALSE, as.matrix = TRUE)
	expect_equal(obtained, t(expect[reduced.cols,, drop=FALSE]))
	
	# With matrix and drop = TRUE
	obtained <- coords(smooth.s100b, "best", ret = reduced.cols, transpose = FALSE, as.matrix = TRUE, drop = TRUE)
	expect_equal(obtained, expect[reduced.cols,])
	
	# Default drop with numeric
	obtained <- coords(smooth.s100b, c(0.2, 0.5), input = "specificity", ret="se")
	expect_is(obtained, "data.frame")
	
	# With numeric x
	obtained <- coords(smooth.s100b, c(0.2, 0.5, 0.6), input = "specificity", transpose = FALSE)
	expect_is(obtained, "data.frame")
	expect_equal(dim(obtained), c(3, 2))
	
})


test_that("coords works with smooth.roc and x = numeric", {
	smooth.s100b <- smooth(r.s100b)
	expect <- structure(c(0.5, 0.797749392103789, 0.608032965276596, 36, 32.7077250762554, 
						  8.29227492374464, 36, 0.812782817364406, 0.476041450069183, 0.523958549930817,
						  0.5, 0.797749392103789, 0.5, 0.202250607896211, 0.5, 0.202250607896211, 
						  0.391967034723404, 0.187217182635594, 0.523958549930817, 0.476041450069183, 
						  0.797749392103789, 1.29774939210379, 0.290905308394387, 0.9, 
						  0.412071871553968, 0.722964130386838, 64.8, 16.8949467337127, 
						  24.1050532662873, 7.2, 0.728867456002887, 0.701182157421994, 
						  0.298817842578006, 0.1, 0.412071871553968, 0.9, 0.587928128446032, 
						  0.1, 0.587928128446032, 0.277035869613162, 0.271132543997113, 
						  0.298817842578006, 0.701182157421994, 0.412071871553968, 1.31207187155397, 
						  0.355659484218054), .Dim = c(23L, 2L), .Dimnames = list(c("specificity", 
						  	"sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", 
						  	"fdr", "fpr", "tpr", "tnr", "fnr", "1-specificity", "1-sensitivity", 
						  	"1-accuracy", "1-npv", "1-ppv", "precision", "recall", "youden", 
						  	"closest.topleft"), NULL))
	
	reduced.cols <- c("specificity", "sensitivity", "youden")

	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "sp", ret="all", transpose=TRUE)
	expect_equal(obtained, expect)
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "spe", ret=reduced.cols, transpose=TRUE)
	expect_equal(obtained, expect[reduced.cols,])
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret="all", drop = TRUE, transpose=TRUE)
	expect_equal(obtained, expect[, 2])
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret=reduced.cols, drop = TRUE, transpose=TRUE)
	expect_equal(obtained, expect[reduced.cols, 2])
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret="all", drop = FALSE, transpose=TRUE)
	expect_equal(obtained, expect[, 2, drop=FALSE])
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret=reduced.cols, drop = FALSE, transpose=TRUE)
	expect_equal(obtained, expect[reduced.cols, 2, drop=FALSE])
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "specificity", ret="all", as.list = TRUE, drop = TRUE)
	expect_equal(obtained[[1]], as.list(expect[, 1]))
	expect_equal(obtained[[2]], as.list(expect[, 2]))
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "specificity", ret=reduced.cols, as.list = TRUE, drop = TRUE)
	expect_equal(obtained[[1]], as.list(expect[reduced.cols, 1]))
	expect_equal(obtained[[2]], as.list(expect[reduced.cols, 2]))
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret="all", as.list = TRUE, drop = TRUE)
	expect_equal(obtained, as.list(expect[, 2]))
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret=reduced.cols, as.list = TRUE, drop = TRUE)
	expect_equal(obtained, as.list(expect[reduced.cols, 2]))
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret="all", as.list = TRUE, drop = FALSE)
	expect_equal(obtained[[1]], as.list(expect[, 2]))
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret=reduced.cols, as.list = TRUE, drop = FALSE)
	expect_equal(obtained[[1]], as.list(expect[reduced.cols, 2]))
})


test_that("coords works with smooth.roc and x = numeric and input = 'se'", {
	smooth.s100b <- smooth(r.s100b)
	expect <- structure(c(0.844189345484777, 0.5, 0.719306485618619, 60.781632874904, 
						  20.5, 20.5, 11.218367125096, 0.747790499834687, 0.646313220322748, 
						  0.353686779677252, 0.155810654515223, 0.5, 0.844189345484777, 
						  0.5, 0.155810654515223, 0.5, 0.280693514381381, 0.252209500165313, 
						  0.353686779677252, 0.646313220322748, 0.5, 1.34418934548478, 
						  0.274276960060462, 0.293322024198721, 0.9, 0.513444121613345, 
						  21.1191857423079, 36.9, 4.1, 50.8808142576921, 0.837425361710953, 
						  0.420365205222125, 0.579634794777875, 0.706677975801279, 0.9, 
						  0.293322024198721, 0.1, 0.706677975801279, 0.1, 0.486555878386655, 
						  0.162574638289047, 0.579634794777875, 0.420365205222125, 0.9, 
						  1.19332202419872, 0.509393761482593), .Dim = c(23L, 2L), .Dimnames = list(
						  	c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", 
						  	  "fp", "npv", "ppv", "fdr", "fpr", "tpr", "tnr", "fnr", "1-specificity", 
						  	  "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "precision", 
						  	  "recall", "youden", "closest.topleft"), NULL))
	
	reduced.cols <- c("specificity", "sensitivity", "youden")
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "se", ret="all", transpose=TRUE)
	expect_equal(obtained, expect)
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "se", ret=reduced.cols, transpose=TRUE)
	expect_equal(obtained, expect[reduced.cols,])
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret="all", drop = TRUE, transpose=TRUE)
	expect_equal(obtained, expect[, 2])
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret=reduced.cols, drop = TRUE, transpose=TRUE)
	expect_equal(obtained, expect[reduced.cols, 2])
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret="all", drop = FALSE, transpose=TRUE)
	expect_equal(obtained, expect[, 2, drop=FALSE])
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret=reduced.cols, drop = FALSE, transpose=TRUE)
	expect_equal(obtained, expect[reduced.cols, 2, drop=FALSE])
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "se", ret="all", as.list = TRUE, drop = TRUE)
	expect_equal(obtained[[1]], as.list(expect[, 1]))
	expect_equal(obtained[[2]], as.list(expect[, 2]))
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "se", ret=reduced.cols, as.list = TRUE, drop = TRUE)
	expect_equal(obtained[[1]], as.list(expect[reduced.cols, 1]))
	expect_equal(obtained[[2]], as.list(expect[reduced.cols, 2]))
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret="all", as.list = TRUE, drop = TRUE)
	expect_equal(obtained, as.list(expect[, 2]))
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret=reduced.cols, as.list = TRUE, drop = TRUE)
	expect_equal(obtained, as.list(expect[reduced.cols, 2]))
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret="all", as.list = TRUE, drop = FALSE)
	expect_equal(obtained[[1]], as.list(expect[, 2]))
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret=reduced.cols, as.list = TRUE, drop = FALSE)
	expect_equal(obtained[[1]], as.list(expect[reduced.cols, 2]))
})


test_that("coords with x = 'best' takes partial AUC into account", {
	# with sp
	obtained <- coords(r.s100b.partial1, "b", ret="t", transpose=TRUE)
	expect_equal(unname(obtained), 0.475)
	
	# with se
	obtained <- coords(r.s100b.partial2, "b", ret="t", transpose=TRUE)
	expect_equal(unname(obtained), 0.075)
})

test_that("coords with x = 'best' takes partial AUC into account with smooth.roc", {
	# with sp
	obtained <- coords(smooth(r.s100b.partial1), "b", ret="sp", transpose=TRUE)
	expect_equal(unname(obtained), 0.900608847772859)
	
	obtained <- coords(smooth(r.s100b.partial1), "b", ret=c("se", "se", "youden"), transpose=TRUE)
	expect_equal(as.vector(obtained), c(0.410958904109589, 0.410958904109589, 1.311567751882448))
	
	# with se
	obtained <- coords(smooth(r.s100b.partial2), "b", ret="se", transpose=TRUE)
	expect_equal(unname(obtained), 0.900195694716243)
	
	obtained <- coords(smooth(r.s100b.partial2), "b", ret=c("se", "se", "youden"), transpose=TRUE)
	expect_equal(as.vector(obtained), c(0.900195694716243, 0.900195694716243, 1.193053239288330))
})


test_that("coords with x = 'all' takes partial AUC into account", {
	# with sp
	obtained <- coords(r.s100b.partial1, "all", ret="t", transpose=TRUE)
	expect_equal(length(obtained), 7)
	expect_equal(min(obtained), 0.435)
	
	# with se
	obtained <- coords(r.s100b.partial2, "all", ret="t", transpose=TRUE)
	expect_equal(length(obtained), 5)
	expect_equal(max(obtained), 0.075)
})



test_that("coords with x = 'all' takes partial AUC into account with smooth.roc", {
	# with sp
	obtained <- coords(smooth(r.s100b.partial1), "all", ret="sp", transpose=TRUE)
	expect_equal(length(obtained), 139)
	expect_equal(min(obtained), 0.90060885)
	
	# with se
	obtained <- coords(smooth(r.s100b.partial2), "all", ret="se", transpose=TRUE)
	expect_equal(length(obtained), 46)
	expect_equal(min(obtained), 0.90019569)
})


test_that("coords with x = 'local maximas' takes partial AUC into account", {
	# with sp
	obtained <- coords(r.s100b.partial1, "local maximas", ret="t", transpose=TRUE)
	expect_equal(unname(obtained), c(0.435, 0.475, 0.485))
	
	# with se
	obtained <- coords(r.s100b.partial2, "local maximas", ret="t", transpose=TRUE)
	expect_equal(unname(obtained), c(0.065, 0.075))
})

test_that("invalid best.weights", {
	expect_error(coords(r.s100b, "best", best.weights = 1, transpose=FALSE))
	expect_error(coords(r.s100b, "best", best.weights = 0:1, transpose=FALSE))
	expect_error(coords(r.s100b, "best", best.weights = c(0.1, 0.9), transpose=FALSE), NA)
	expect_error(coords(r.s100b, "best", best.weights = 1:3, transpose=FALSE))
	# with smooth
	expect_error(coords(smooth(r.s100b), "best", best.weights = 1, transpose=FALSE))
	expect_error(coords(smooth(r.s100b), "best", best.weights = 0:1, transpose=FALSE))
	expect_error(coords(smooth(r.s100b), "best", best.weights = c(0.1, 0.9), transpose=FALSE), NA)
	expect_error(coords(smooth(r.s100b), "best", best.weights = 1:3, transpose=FALSE))
})

test_that("invalid best.method", {
	expect_error(coords(r.s100b, "best", best.method = 1, transpose=FALSE))
	expect_error(coords(r.s100b, "best", best.method = "1", transpose=FALSE))
	# with smooth
	expect_error(coords(smooth(r.s100b), "best", best.method = 1, transpose=FALSE))
	expect_error(coords(smooth(r.s100b), "best", best.method = "1", transpose=FALSE))
})

test_that("invalid se/sp", {
	smooth.s100b <- smooth(r.s100b)
	for (inp in c("sens", "spec")) {
		for (r in list(r.s100b, smooth.s100b)) {
			expect_error(coords(r, x=-2, input=inp, transpose=FALSE))
			expect_error(coords(r, x=0, input=inp, transpose=FALSE), NA)
			expect_error(coords(r, x=1, input=inp, transpose=FALSE), NA)
			expect_error(coords(r, x=10, input=inp, transpose=FALSE))
		}
	}
	smooth.s100b.percent <- smooth(r.s100b.percent)
	for (inp in c("sens", "spec")) {
		for (r in list(r.s100b.percent, smooth.s100b.percent)) {
			expect_error(coords(r.s100b.percent, x=-2, input=inp, transpose=FALSE))
			expect_error(coords(r.s100b.percent, x=0, input=inp, transpose=FALSE), NA)
			expect_error(coords(r.s100b.percent, x=10, input=inp, transpose=FALSE), NA)
			expect_error(coords(r.s100b.percent, x=100, input=inp, transpose=FALSE), NA)
			expect_error(coords(r.s100b.percent, x=101, input=inp, transpose=FALSE))
		}
	}
})

test_that("invalid x", {
	expect_error(coords(r.s100b.percent, x=list(1), transpose=FALSE))
	expect_error(coords(r.s100b, x=aSAH, transpose=FALSE))
	expect_error(coords(smooth(r.s100b), x=mean, transpose=FALSE))
	# character but invalid
	expect_error(coords(smooth(r.s100b), x="c", transpose=FALSE))
	expect_error(coords(r.s100b, x="c", transpose=FALSE))
})

test_that("Infinite values work with both directions", {
	# direction = >
	Data <- structure(list(Outcome = c(1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L), Value = c(72L, 65L, 271L, 73L, 87L, 114L, 111L, 47L, 88L, 44L, 121L, 207L, 33L, 138L, 284L, 62L, 120L, 116L, 202L, 172L, 117L, 69L, 102L, 150L, 131L, 77L, 124L, 46L, 579L, 117L, 96L, 83L, 102L)), class = "data.frame", row.names = c(NA, -33L))
	ROC <- roc(Outcome~Value, data=Data, ci=TRUE, direction=">")
	co <- coords(ROC, x=c(-Inf, Inf), transpose = FALSE)
	expect_equivalent(co, data.frame(threshold = c(-Inf, Inf), specificity = c(1, 0), sensitivity = c(0, 1)))
	
	# direction = <
	co <- coords(r.s100b, x=c(-Inf, Inf), transpose = FALSE)
	expect_equivalent(co, data.frame(threshold = c(-Inf, Inf), specificity = c(0, 1), sensitivity = c(1, 0)))
})

test_that("Coords pick the right end of 'flat' bits of the curve, according to direction", {
	# expect_equal(r.s100b$sensitivities[2], 0.975609756097561) # tested elsewhere
	expect_equivalent(
		coords(r.s100b, 0.975609756097561, "se", "sp", transpose = TRUE),
		0.13888888888888889 # and not 0
	)
	expect_equivalent(
		coords(r.s100b, 1, "sp", "se", transpose = TRUE),
		0.2926829268292683 # and not 0
	)
})
