library(pROC)
data(aSAH)

context("coords")

test_that("coords with thresholds works", {
	return.rows <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "lr_pos", "lr_neg", "youden", "closest.topleft")
	obtained <- coords(r.s100b, "all", ret = return.rows)
	expect_equal(obtained, expected.coords[,return.rows])
})

test_that("coords returns all thresholds by default", {
	obtained <- coords(r.s100b)
	expect_equal(obtained, expected.coords[,c("threshold", "specificity", "sensitivity")])
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


test_that("coords with transpose = TRUE works", {
	suppressWarnings({
		return.rows <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "youden", "closest.topleft")
		obtained <- coords(r.s100b, "all", ret = return.rows, transpose = TRUE)
		expect_equal(obtained, t(expected.coords[,return.rows]))
		obtained <- coords(r.s100b, transpose = TRUE)
		expect_equal(obtained, t(expected.coords[,c("threshold", "specificity", "sensitivity")]))
		
		# With drop=TRUE
		obtained <- coords(r.s100b, "all", ret = "se", transpose = TRUE, drop=TRUE)
		expect_is(obtained, "numeric")
		#  Not why drop.data.frame returns a list, skipping
		# obtained <- coords(r.s100b, "best", ret = "all", transpose = FALSE, drop=TRUE)
		
		# With drop=FALSE
		obtained <- coords(r.s100b, "all", ret = "se", transpose = TRUE, drop=FALSE)
		expect_is(obtained, "matrix")
	})
})


test_that("coords with ret='all' works", {
	obtained <- coords(r.s100b, "all", ret = "all")
	expect_equal(dim(obtained), c(51, 26))
	expect_equal(obtained[,colnames(expected.coords)], expected.coords)
})


test_that("coords with ret='all' doesn't accept additional options", {
	expect_error(coords(r.s100b, "all", ret = c("all", "thresholds")))
})


test_that("coords with percent works", {
	return.rows <- "all"
	percent.cols <- c("specificity", "sensitivity", "accuracy", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "youden", "closest.topleft", "fdr", "fpr", "tpr", "tnr", "fnr", "precision", "recall")
	obtained.percent <- coords(r.s100b.percent, "all", ret = return.rows)
	# Adjust for percent
	obtained.percent[,percent.cols] <- obtained.percent[,percent.cols] / 100
	expect_equal(obtained.percent, expected.coords)
})


test_that("coords with local maximas thresholds works", {
	return.rows <- "all"
	obtained <- coords(r.s100b, "local maximas", ret = return.rows)
	#expected.thresholds <- c(-Inf, 0.065, 0.075, 0.085, 0.095, 0.105, 0.115, 0.135, 0.155, 0.205, 0.245, 0.29, 0.325, 0.345, 0.395, 0.435, 0.475, 0.485, 0.51)
	expected.thresholds = c(-Inf, 0x1.0a3d70a3d70a4p-4, 0x1.3333333333334p-4, 
				 0x1.5c28f5c28f5c2p-4, 0x1.851eb851eb852p-4, 0x1.ae147ae147ae2p-4, 
				 0x1.d70a3d70a3d7p-4, 0x1.147ae147ae148p-3, 0x1.3d70a3d70a3d7p-3, 
				 0x1.a3d70a3d70a3ep-3, 0x1.f5c28f5c28f5cp-3, 0x1.28f5c28f5c29p-2, 
				 0x1.4cccccccccccdp-2, 0x1.6147ae147ae14p-2, 0x1.947ae147ae148p-2, 
				 0x1.bd70a3d70a3d7p-2, 0x1.e666666666666p-2, 0x1.f0a3d70a3d70ap-2, 
				 0x1.051eb851eb852p-1)
	expect_equal(as.vector(obtained[,"threshold"]), expected.thresholds)
	expect_equivalent(obtained, expected.coords[expected.coords[,"threshold"] %in% expected.thresholds,])
})


test_that("coords with best threshold works", {
	return.rows <- "all"
	obtained <- coords(r.s100b, "best", ret = return.rows)
	expect_equivalent(obtained, expected.coords[abs(expected.coords[,"threshold"] - 0.205) < 0.001, , drop=FALSE])
})


test_that("coords with arbitrary thresholds works", {
	return.rows <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp",  "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv", "youden", "closest.topleft")
	obtained <- coords(r.s100b, c(0.205, 0.055), input = "threshold", ret = return.rows)
	expect_equivalent(obtained, expected.coords[c(18, 4), return.rows])
})

test_that("coords with arbitrary thresholds at exact data point works", {
	return.rows <- "all"
	expect_equal(sum(aSAH$s100b == 0.05),  3)
	expect_equal(sum(aSAH$s100b == 0.52),  1)
	obtained <- coords(r.s100b, c(0.05, 0.52), input = "threshold", ret = return.rows)
	expect_equivalent(obtained[,-1], expected.coords[c(3, 40), -1])
})

test_that("coords with arbitrary thresholds works with direction=>", {
	obtained <- coords(r.s100b.reversed, c(0.05, 0.055, 0.205, 0.52), input = "threshold", ret = "all")
	expect_equivalent(obtained, expected.coords.reverse)
})


test_that("coords with single arbitrary threshold works", {
	return.rows <- "all"
	obtained <- coords(r.s100b, c(0.205), input = "threshold", ret = return.rows)
	expect_equivalent(obtained, expected.coords[18, , drop=FALSE])
})


test_that("coords with arbitrary thresholds at exact data point works", {
	expect_equal(sum(aSAH$s100b == 0.05),  3)
	expect_equal(sum(aSAH$s100b == 0.52),  1)
	obtained <- coords(r.s100b, c(0.05), input = "threshold", ret = "all")
	expect_equivalent(obtained[,-1], expected.coords[3,-1])
	obtained <- coords(r.s100b, c(0.52), input = "threshold", ret = "all")
	expect_equivalent(obtained[,-1], expected.coords[40, -1])
})


test_that("coords with arbitrary thresholds works with direction=>", {
	obtained <- coords(r.s100b.reversed, c(0.05), input = "threshold", ret = "all")
	expect_equivalent(obtained, expected.coords.reverse[1,])
	obtained <- coords(r.s100b.reversed, c(0.055), input = "threshold", ret = "all")
	expect_equivalent(obtained, expected.coords.reverse[2,])
	obtained <- coords(r.s100b.reversed, c(0.205), input = "threshold", ret = "all")
	expect_equivalent(obtained, expected.coords.reverse[3, ])
	obtained <- coords(r.s100b.reversed, c(0.52), input = "threshold", ret = "all")
	expect_equivalent(obtained, expected.coords.reverse[4,])
})


test_that("coords with sensitivity works", {
	obtained <- coords(r.s100b, seq(0, 1, .1), input = "sensitivity", ret = c("threshold", "specificity", "sensitivity"))
	expect_equal(obtained[,"threshold"], c(Inf, rep(NA, 9), -Inf))
	expect_equal(obtained[,"sensitivity"], seq(0, 1, .1))
	expect_equal(obtained[,"specificity"], c(1, 1, 1, 0.972222222222222, 0.888888888888889, 0.833333333333333, 0.805555555555556, 0.56875, 0.447222222222222, 0.230555555555556, 0))
})


test_that("coords with sensitivity works with percent", {
	obtained <- coords(r.s100b.percent, seq(0, 100, 10), input = "sensitivity", ret = c("threshold", "specificity", "sensitivity"))
	expect_equal(obtained[,"threshold"], c(Inf, rep(NA, 9), -Inf))
	expect_equal(obtained[,"sensitivity"], seq(0, 100, 10))
	expect_equal(obtained[,"specificity"], c(1, 1, 1, 0.972222222222222, 0.888888888888889, 0.833333333333333, 0.805555555555556, 0.56875, 0.447222222222222, 0.230555555555556, 0) * 100)
})


test_that("coords with specificity works", {
	obtained <- coords(r.s100b, seq(0, 1, .1), input = "specificity", ret = c("threshold", "specificity", "sensitivity"))
	expect_equal(obtained[,"threshold"], c(-Inf, rep(NA, 9), 0.51))
	expect_equal(obtained[,"specificity"], seq(0, 1, .1))
	expect_equal(obtained[,"sensitivity"], c(1, 0.975609756097561, 0.921951219512195, 0.879674796747967, 0.823693379790941, 0.774390243902439, 0.675609756097561, 0.655284552845528, 0.634146341463415, 0.390243902439024, 0.292682926829268))
})


test_that("coords with specificity works with percent", {
	obtained <- coords(r.s100b.percent, seq(0, 100, 10), input = "specificity", ret = c("threshold", "specificity", "sensitivity"))
	expect_equal(obtained[,"threshold"], c(-Inf, rep(NA, 9), 0.51))
	expect_equal(obtained[,"specificity"], seq(0, 100, 10))
	expect_equal(obtained[,"sensitivity"], c(1, 0.975609756097561, 0.921951219512195, 0.879674796747967, 0.823693379790941, 0.774390243902439, 0.675609756097561, 0.655284552845528, 0.634146341463415, 0.390243902439024, 0.292682926829268) * 100)
})



test_that("coords with specificity works with as.list", {
	expect_warning(obtained <- coords(r.s100b.percent, "best", ret = c("threshold", "specificity", "accuracy"), as.list = TRUE), "as.list")
	expect_equal(obtained, list(
		threshold = 0.205,
		specificity = unname(expected.coords[18, "specificity"]) * 100,
		accuracy = unname(expected.coords[18, "accuracy"]) * 100
	))
})

test_that("coords with specificity works with as.list and drop=FALSE", {
	expect_warning(obtained <- coords(r.s100b.percent, "best", 
					   ret = c("threshold", "specificity", "accuracy"), 
					   as.list = TRUE, drop = FALSE), "as.list")
	expect_equal(obtained[[1]], list(
		threshold = 0.205,
		specificity = unname(expected.coords[18, "specificity"]) * 100,
		accuracy = unname(expected.coords[18, "accuracy"]) * 100
	))
})


test_that("coords with specificity works with as.list and several thresholds", {
	expect_warning(obtained <- coords(r.s100b.percent, c(0.205, 0.51), 
					   ret = c("threshold", "specificity", "accuracy"), 
					   as.list = TRUE, drop = FALSE), "as.list")
	expect_equal(obtained[[1]], list(
		threshold = 0.205,
		specificity = unname(expected.coords[18, "specificity"]) * 100,
		accuracy = unname(expected.coords[18, "accuracy"]) * 100
	))
	expect_equal(obtained[[2]], list(
		threshold = 0.51,
		specificity = unname(expected.coords[40, "specificity"]) * 100,
		accuracy = unname(expected.coords[40, "accuracy"]) * 100
	))
})


test_that("drop works", {
	suppressWarnings({
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
})


test_that("as.matrix works", {
	suppressWarnings({
		obtained <- coords(r.s100b, c(0.51, 0.205), ret="sensitivity", as.matrix = TRUE)
		expect_equivalent(obtained, as.matrix(expected.coords[c(40, 18), "sensitivity", drop = FALSE]))
	})
})


test_that("as.matrix works with drop=TRUE", {
	suppressWarnings({
		obtained <- coords(r.s100b, c(0.51, 0.205), ret="sensitivity", as.matrix = TRUE, drop = TRUE)
		expect_equal(obtained, expected.coords[c(40, 18), "sensitivity", drop = TRUE])
	})
})


test_that("coords returns the correct basic values ", {
	obtained <- coords(r.s100b, 0.205, 
					   ret = c("t", "tp", "fp", "tn", "fn",
					   		   "sp", "se", "acc",
					   		   "npv", "ppv", "precision", "recall",
					   		   "tpr", "fpr", "tnr", "fnr", "fdr"))
	
	obtained.percent <- coords(r.s100b.percent, 0.205, 
					   ret = c("t", "tp", "fp", "tn", "fn",
					   		"sp", "se", "acc",
					   		"npv", "ppv", "precision", "recall",
					   		"tpr", "fpr", "tnr", "fnr", "fdr"))
	
	# We assume the following values:
	# tp fp tn fn N
	# 26 14 58 15 113
	
	expected <- data.frame(
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
	
	expect_equivalent(obtained, expected)
	expect_equivalent(obtained.percent[,1:5], expected[,1:5])
	expect_equivalent(obtained.percent[,6:17], expected[,6:17]*100)
})


test_that("coords works with smooth.roc and x = 'best' and transpose=TRUE", {
	suppressWarnings({
		smooth.s100b <- smooth(r.s100b)
		expect <- structure(c(0.750857175922901, 0.608610567514677, 0.699245574642041, 
							  54.0617166664488, 24.9530332681018, 16.0469667318982, 17.9382833335512, 
							  0.771112992655678, 0.581773544045047, 0.418226455954953, 0.249142824077099, 
							  0.608610567514677, 0.750857175922901, 0.391389432485323, 0.249142824077099, 
							  0.391389432485323, 0.300754425357959, 0.228887007344322, 0.418226455954953, 
							  0.581773544045047, 0.608610567514677, 2.4428179690470926, 0.52125683157286817,
							  1.35946774343758, 0.215257834650296
		), .Dim = c(25L, 1L), .Dimnames = list(c("specificity", "sensitivity", 
								"accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "fdr", "fpr", 
								"tpr", "tnr", "fnr", "1-specificity", "1-sensitivity", "1-accuracy", 
								"1-npv", "1-ppv", "precision", "recall", "lr_pos", "lr_neg", 
								"youden", "closest.topleft"
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
		
		expect_warning(obtained <- coords(smooth.s100b, "best", ret = "all", as.list = TRUE), "as.list")
		expect_equal(obtained, as.list(expect[, 1]))
		expect_equal(names(obtained), rownames(expect))
		
		expect_warning(obtained <- coords(smooth.s100b, "best", ret = "all", as.list = TRUE, drop = FALSE), "as.list")
		expect_equal(obtained[[1]], as.list(expect[, 1])) # names
		expect_equal(names(obtained[[1]]), rownames(expect))
		
		expect_warning(obtained <- coords(smooth.s100b, "best", ret = reduced.cols, as.list = TRUE), "as.list")
		expect_equal(obtained, as.list(expect[reduced.cols, 1]))
		expect_equal(names(obtained), reduced.cols)
		
		expect_warning(obtained <- coords(smooth.s100b, "best", ret = reduced.cols, as.list = TRUE, drop = FALSE), "as.list")
		expect_equal(obtained[[1]], as.list(expect[reduced.cols, 1])) # names
		expect_equal(names(obtained[[1]]), reduced.cols)
	})
})


test_that("coords works with smooth.roc", {
	suppressWarnings({
		smooth.s100b <- smooth(r.s100b)
		expect <- structure(c(0.750857175922901, 0.608610567514677, 0.699245574642041, 
							  54.0617166664488, 24.9530332681018, 16.0469667318982, 17.9382833335512, 
							  0.771112992655678, 0.581773544045047, 0.418226455954953, 0.249142824077099, 
							  0.608610567514677, 0.750857175922901, 0.391389432485323, 0.249142824077099, 
							  0.391389432485323, 0.300754425357959, 0.228887007344322, 0.418226455954953, 
							  0.581773544045047, 0.608610567514677, 2.4428179690470926, 0.52125683157286817,
							  1.35946774343758, 0.215257834650296
		), .Dim = c(25L, 1L), .Dimnames = list(c("specificity", "sensitivity", 
												 "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "fdr", "fpr", 
												 "tpr", "tnr", "fnr", "1-specificity", "1-sensitivity", "1-accuracy", 
												 "1-npv", "1-ppv", "precision", "recall", "lr_pos", "lr_neg",
												 "youden", "closest.topleft"
		), NULL))
		
		
		reduced.cols <- c("specificity", "sensitivity", "youden")
		
		obtained <- coords(smooth.s100b, "best", ret = reduced.cols)
		expect_equal(obtained, as.data.frame(t(expect[reduced.cols,, drop=FALSE])))
		
		obtained <- coords(smooth.s100b, "best", ret = "all", drop = FALSE)
		expect_equal(obtained, as.data.frame(t(expect)))
		
		# Without drop
		obtained <- coords(smooth.s100b, "best", ret = reduced.cols)
		expect_equivalent(obtained, as.data.frame(t(expect[reduced.cols,])))
		
		# drop = TRUE
		obtained <- coords(smooth.s100b, "best", ret = reduced.cols, drop = TRUE)
		expect_equal(obtained, as.list(expect[reduced.cols,]))
		
		# With as.matrix
		obtained <- coords(smooth.s100b, "best", ret = reduced.cols, as.matrix = TRUE)
		expect_equal(obtained, t(expect[reduced.cols,, drop=FALSE]))
		
		# With matrix and drop = TRUE
		obtained <- coords(smooth.s100b, "best", ret = reduced.cols, as.matrix = TRUE, drop = TRUE)
		expect_equal(obtained, expect[reduced.cols,])
		
		# Default drop with numeric
		obtained <- coords(smooth.s100b, c(0.2, 0.5), input = "specificity", ret="se")
		expect_is(obtained, "data.frame")
		
		# With numeric x
		obtained <- coords(smooth.s100b, c(0.2, 0.5, 0.6), input = "specificity")
		expect_is(obtained, "data.frame")
		expect_equal(dim(obtained), c(3, 2))
	})
	
})


test_that("coords works with smooth.roc and x = numeric", {
	smooth.s100b <- smooth(r.s100b)
	expect <- structure(list(specificity = c(0.5, 0.90000000000000002), sensitivity = c(0.79774939210378937, 
0.41207187155396763), accuracy = c(0.60803296527659623, 0.72296413038683782
), tn = c(36, 64.799999999999997), tp = c(32.707725076255365, 
16.894946733712672), fn = c(8.2922749237446354, 24.105053266287328
), fp = c(36, 7.2000000000000028), npv = c(0.81278281736440605, 
0.72886745600288694), ppv = c(0.47604145006918286, 0.70118215742199363
), fdr = c(0.52395854993081703, 0.29881784257800637), fpr = c(0.5, 
0.099999999999999978), tpr = c(0.79774939210378937, 0.41207187155396763
), tnr = c(0.5, 0.90000000000000002), fnr = c(0.20225060789621063, 
0.58792812844603237), `1-specificity` = c(0.5, 0.099999999999999978
), `1-sensitivity` = c(0.20225060789621063, 0.58792812844603237
), `1-accuracy` = c(0.39196703472340377, 0.27703586961316218), 
    `1-npv` = c(0.18721718263559395, 0.27113254399711306), `1-ppv` = c(0.52395854993081714, 
    0.29881784257800637), precision = c(0.47604145006918286, 
    0.70118215742199363), recall = c(0.79774939210378937, 0.41207187155396763
    ), lr_pos = c(1.5954987842075787, 4.1207187155396774), lr_neg = c(0.40450121579242126, 
    0.65325347605114703), youden = c(1.2977493921037895, 1.3120718715539677
    ), closest.topleft = c(0.29090530839438672, 0.35565948421805432
    )), class = "data.frame", row.names = c(NA, -2L))
	
	reduced.cols <- c("specificity", "sensitivity", "youden")

	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "sp", ret="all")
	expect_equal(obtained, expect)
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "spe", ret=reduced.cols)
	expect_equal(obtained, expect[,reduced.cols])
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret="all")
	expect_equivalent(obtained, expect[2,, drop=FALSE])
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret=reduced.cols)
	expect_equivalent(obtained, expect[2, reduced.cols])
	
	obtained <- coords(smooth.s100b, 0.9, input = "specificity", ret=reduced.cols)
	expect_equivalent(obtained, expect[2, reduced.cols, drop=FALSE])
})


test_that("coords works with smooth.roc and x = numeric and input = 'se'", {
	smooth.s100b <- smooth(r.s100b)
	expect <- structure(list(specificity = c(0.84418934548477731, 0.29332202419872122
), sensitivity = c(0.5, 0.90000000000000002), accuracy = c(0.7193064856186191, 
0.51344412161334452), tn = c(60.781632874903963, 21.119185742307927
), tp = c(20.5, 36.899999999999999), fn = c(20.5, 4.1000000000000014
), fp = c(11.218367125096037, 50.880814257692073), npv = c(0.74779049983468704, 
0.83742536171095305), ppv = c(0.64631322032274796, 0.42036520522212539
), fdr = c(0.35368677967725209, 0.57963479477787461), fpr = c(0.15581065451522269, 
0.70667797580127878), tpr = c(0.5, 0.90000000000000002), tnr = c(0.84418934548477731, 
0.29332202419872122), fnr = c(0.5, 0.10000000000000003), `1-specificity` = c(0.15581065451522269, 
0.70667797580127878), `1-sensitivity` = c(0.5, 0.099999999999999978
), `1-accuracy` = c(0.2806935143813809, 0.48655587838665548), 
    `1-npv` = c(0.25220950016531296, 0.16257463828904695), `1-ppv` = c(0.35368677967725204, 
    0.57963479477787461), precision = c(0.64631322032274796, 
    0.42036520522212539), recall = c(0.5, 0.90000000000000002
    ), lr_pos = c(3.2090231669693039, 1.2735645241802247), lr_neg = c(0.59228418680512263, 
    0.34092223477992739), youden = c(1.3441893454847773, 1.1933220241987212
    ), closest.topleft = c(0.27427696006046209, 0.50939376148259274
    )), class = "data.frame", row.names = c(NA, -2L))
	
	reduced.cols <- c("specificity", "sensitivity", "youden")
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "se", ret="all")
	expect_equal(obtained, expect)
	
	obtained <- coords(smooth.s100b, c(0.5, 0.9), input = "se", ret=reduced.cols)
	expect_equal(obtained, expect[,reduced.cols])
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret="all")
	expect_equivalent(obtained, expect[2,, drop=FALSE])
	
	obtained <- coords(smooth.s100b, 0.9, input = "se", ret=reduced.cols)
	expect_equivalent(obtained, expect[2, reduced.cols, drop=FALSE])
})


test_that("coords with x = 'best' takes partial AUC into account", {
	# with sp
	obtained <- coords(r.s100b.partial1, "b", ret="t")
	expect_equal(obtained$threshold, 0.475)
	
	# with se
	obtained <- coords(r.s100b.partial2, "b", ret="t")
	expect_equal(obtained$threshold, 0.075)
})

test_that("coords with x = 'best' takes partial AUC into account with smooth.roc", {
	# with sp
	obtained <- coords(smooth(r.s100b.partial1), "b", ret="sp")
	expect_equal(obtained$specificity, 0.900608847772859)
	
	obtained <- coords(smooth(r.s100b.partial1), "b", ret=c("se", "se", "youden"))
	expect_equivalent(as.vector(obtained), c(0.410958904109589, 0.410958904109589, 1.311567751882448))
	
	# with se
	obtained <- coords(smooth(r.s100b.partial2), "b", ret="se")
	expect_equal(obtained$sensitivity, 0.900195694716243)
	
	obtained <- coords(smooth(r.s100b.partial2), "b", ret=c("se", "se", "youden"))
	expect_equivalent(as.vector(obtained), c(0.900195694716243, 0.900195694716243, 1.193053239288330))
})


test_that("coords with x = 'all' takes partial AUC into account", {
	# with sp
	obtained <- coords(r.s100b.partial1, "all", ret=c("t", "se", "sp"))
	expect_equal(dim(obtained), c(9, 3))
	expect_equivalent(obtained[1,], c(NA, 0.3902439, 0.9))
	expect_equivalent(obtained[9,], c(NA, 0.292682926, 0.99))
	
	# with se
	obtained <- coords(r.s100b.partial2, "all", ret=c("t", "se", "sp"))
	expect_equal(dim(obtained), c(7, 3))
	expect_equivalent(obtained[1,], c(NA, 0.99, 0.0))
	expect_equivalent(obtained[7,], c(NA, 0.9, 0.230555555555556))
})

test_that("coords with ignore.partial.auc = TRUE ignores partial AUC", {
	# with sp
	obtained <- coords(r.s100b.partial1, "all", ret="t", ignore.partial.auc=TRUE)
	expect_equal(dim(obtained), c(51, 1))

	# with se
	obtained <- coords(r.s100b.partial2, "all", ret="t", ignore.partial.auc=TRUE)
	expect_equal(dim(obtained), c(51, 1))
})


test_that("coords with x = 'all' takes partial AUC into account with smooth.roc", {
	# with sp
	obtained <- coords(smooth(r.s100b.partial1), "all", ret="sp")
	expect_equal(dim(obtained), c(141, 1))
	expect_equal(min(obtained), 0.9)
	expect_equal(max(obtained), 0.99)
	
	# with se
	obtained <- coords(smooth(r.s100b.partial2), "all", ret="se")
	expect_equal(dim(obtained), c(48, 1))
	expect_equal(min(obtained), 0.9)
	expect_equal(max(obtained), 0.99)
})

test_that("coords with ignore.partial.auc = TRUE ignores partial AUC of smooth.roc", {
	# with sp
	obtained <- coords(smooth(r.s100b.partial1), "all", ret="sp", ignore.partial.auc=TRUE)
	expect_equal(dim(obtained), c(514, 1))
	
	# with se
	obtained <- coords(smooth(r.s100b.partial2), "all", ret="se", ignore.partial.auc=TRUE)
	expect_equal(dim(obtained), c(514, 1))
})


test_that("coords with x = 'local maximas' takes partial AUC into account", {
	# with sp
	obtained <- coords(r.s100b.partial1, "local maximas", ret="t")
	expect_equal(obtained$threshold, c(0.435, 0.475, 0.485))
	
	# with se
	obtained <- coords(r.s100b.partial2, "local maximas", ret="t")
	expect_equal(obtained$threshold, c(0.065, 0.075))
})

test_that("invalid best.weights", {
	expect_error(coords(r.s100b, "best", best.weights = 1))
	expect_error(coords(r.s100b, "best", best.weights = 0:1))
	expect_error(coords(r.s100b, "best", best.weights = c(0.1, 0.9)), NA)
	expect_error(coords(r.s100b, "best", best.weights = 1:3))
	# with smooth
	expect_error(coords(smooth(r.s100b), "best", best.weights = 1))
	expect_error(coords(smooth(r.s100b), "best", best.weights = 0:1))
	expect_error(coords(smooth(r.s100b), "best", best.weights = c(0.1, 0.9)), NA)
	expect_error(coords(smooth(r.s100b), "best", best.weights = 1:3))
})

test_that("invalid best.method", {
	expect_error(coords(r.s100b, "best", best.method = 1))
	expect_error(coords(r.s100b, "best", best.method = "1"))
	# with smooth
	expect_error(coords(smooth(r.s100b), "best", best.method = 1))
	expect_error(coords(smooth(r.s100b), "best", best.method = "1"))
})

test_that("invalid se/sp", {
	smooth.s100b <- smooth(r.s100b)
	for (inp in c("sens", "spec")) {
		for (r in list(r.s100b, smooth.s100b)) {
			expect_error(coords(r, x=-2, input=inp))
			expect_error(coords(r, x=0, input=inp), NA)
			expect_error(coords(r, x=1, input=inp), NA)
			expect_error(coords(r, x=10, input=inp))
		}
	}
	smooth.s100b.percent <- smooth(r.s100b.percent)
	for (inp in c("sens", "spec")) {
		for (r in list(r.s100b.percent, smooth.s100b.percent)) {
			expect_error(coords(r.s100b.percent, x=-2, input=inp))
			expect_error(coords(r.s100b.percent, x=0, input=inp), NA)
			expect_error(coords(r.s100b.percent, x=10, input=inp), NA)
			expect_error(coords(r.s100b.percent, x=100, input=inp), NA)
			expect_error(coords(r.s100b.percent, x=101, input=inp))
		}
	}
})

test_that("invalid x", {
	expect_error(coords(r.s100b.percent, x=list(1)))
	expect_error(coords(r.s100b, x=aSAH))
	expect_error(coords(smooth(r.s100b), x=mean))
	# character but invalid
	expect_error(coords(smooth(r.s100b), x="c"))
	expect_error(coords(r.s100b, x="c"))
})

test_that("Infinite values work with both directions", {
	# direction = >
	Data <- structure(list(Outcome = c(1L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L), Value = c(72L, 65L, 271L, 73L, 87L, 114L, 111L, 47L, 88L, 44L, 121L, 207L, 33L, 138L, 284L, 62L, 120L, 116L, 202L, 172L, 117L, 69L, 102L, 150L, 131L, 77L, 124L, 46L, 579L, 117L, 96L, 83L, 102L)), class = "data.frame", row.names = c(NA, -33L))
	ROC <- roc(Outcome~Value, data=Data, ci=TRUE, direction=">")
	co <- coords(ROC, x=c(-Inf, Inf))
	expect_equivalent(co, data.frame(threshold = c(-Inf, Inf), specificity = c(1, 0), sensitivity = c(0, 1)))
	
	# direction = <
	co <- coords(r.s100b, x=c(-Inf, Inf))
	expect_equivalent(co, data.frame(threshold = c(-Inf, Inf), specificity = c(0, 1), sensitivity = c(1, 0)))
})

test_that("Coords pick the right end of 'flat' bits of the curve, according to direction", {
		# expect_equal(r.s100b$sensitivities[2], 0.975609756097561) # tested elsewhere
		expect_equivalent(
			coords(r.s100b, 0.975609756097561, "se", "sp"),
			0.13888888888888889 # and not 0
		)
		expect_equivalent(
			coords(r.s100b, 1, "sp", "se"),
			0.2926829268292683 # and not 0
		)
})

