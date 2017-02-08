library(pROC)
data(aSAH)

# Get static, correct output of the function
source(system.file("extdata", "test-roc-expected.R", package="pROC"))

level.values <- list(
	forward = c("Good", "Poor"),
	reversed = c("Poor", "Good")
)

for (marker in c("ndka", "wfns", "s100b")) {
	for (levels.direction in names(level.values)) {
		for (percent in c(FALSE, TRUE)) {
			for (direction in c("auto", "<", ">")) {
				for (algorithm in 1:4) {
					context(sprintf("'roc' function works with percent = %s, marker = %s, levels.direction = %s, direction = %s and algorithm = %s", percent, marker, levels.direction, direction, algorithm))
					r <- roc(aSAH$outcome, aSAH[[marker]], levels = level.values[[levels.direction]], direction = direction, percent = percent, quiet = TRUE)
					
					test_that("roc.formula produces the same results as roc.default", {
						rf <- roc(as.formula(sprintf("outcome ~ %s", marker)), data = aSAH, levels = level.values[[levels.direction]], direction = direction, percent = percent, quiet = TRUE)
						expect_is(rf, "roc")
						for (item in c("auc", "percent", "sensitivities", "specificities", "thresholds", "direction", "cases", "controls", "fun.sesp")) {
							expect_identical(rf[[item]], r[[item]], label = sprintf("roc(outcome ~ %s, %s, %s, %s, %s)[[\"%s\"]]", marker, levels.direction, percent, direction, algorithm, item))
						}
						for (item in c("original.predictor", "original.response", "predictor", "response", "levels")) {
							expect_identical(unname(rf[[item]]), unname(r[[item]]), label = sprintf("roc(outcome ~ %s, %s, %s, %s, %s)[[\"%s\"]]", marker, levels.direction, percent, direction, algorithm, item))
						}
					})
					
					test_that("roc.default works with control/cases as well", {
						rcs <- roc(controls = r$controls, cases = r$cases, levels = level.values[[levels.direction]], direction = direction, percent = percent, quiet = TRUE)
						expect_is(rcs, "roc")
						for (item in c("auc", "percent", "sensitivities", "specificities", "thresholds", "direction", "cases", "controls", "fun.sesp")) {
							expect_identical(rcs[[item]], r[[item]])
						}
					})
						


					test_that("roc.default produces the expected results", {
						expected.direction <- ifelse(direction == "auto", ifelse(levels.direction == "forward", "<", ">"), direction)
					
						expect_is(r, "roc")
						expect_identical(r$percent, percent)
						expect_identical(r$direction, expected.direction)
						expect_identical(r$levels, level.values[[levels.direction]])
						
						expect_equal(r$thresholds, expected.roc[[marker]][[levels.direction]][[expected.direction]][["thresholds"]])
						expect_equal(r$sensitivities, expected.roc[[marker]][[levels.direction]][[expected.direction]][["sensitivities"]] * ifelse(percent, 100, 1))
						expect_equal(r$specificities, expected.roc[[marker]][[levels.direction]][[expected.direction]][["specificities"]] * ifelse(percent, 100, 1))
					})
				}
			}
		}
	}
}

test_that("roc can't take both response/predictor and case/control", {
	skip("This doesn't throw an error currently.")
	expect_error(roc(aSAH$outcome, aSAH$ndka, controls = aSAH$ndka[aSAH$outcome == "Good"], cases = aSAH$ndka[aSAH$outcome == "Poor"]))
})

# The code below can be used to refresh the "expected.roc" data, just in case...
# expected.roc <- list()
# for (marker in c("ndka", "wfns", "s100b")) {
# 	expected.roc[[marker]] <- list()
# 	for (levels.direction in names(level.values)) {
# 		expected.roc[[marker]][[levels.direction]] <- list()
# 		for (direction in c("<", ">")) {
# 			r <- roc(aSAH$outcome, aSAH[[marker]], levels = level.values[[levels.direction]], direction = direction, percent = FALSE, quiet = TRUE)
# 			if (!isTRUE(percent) && direction != "auto") {
# 				expected.roc[[marker]][[levels.direction]][[direction]] <- r[c("sensitivities", "specificities", "thresholds")]
# 				expected.roc[[marker]][[levels.direction]][[direction]][["auc"]] <- as.numeric(r$auc)
# 			}
# 		}
# 	}
# }
# save("expected.roc", system.file("extdata", "test-roc-expected.R", package="pROC"), file = "dump_roc_expected.R")

