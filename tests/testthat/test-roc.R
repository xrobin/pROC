library(pROC)
data(aSAH)

level.values <- list(
  forward = c("Good", "Poor"),
  reversed = c("Poor", "Good")
)

smooth.methods <- c("binormal", "density", "fitdistr", "logcondens", "logcondens.smooth")

for (marker in c("ndka", "wfns", "s100b")) {
  for (levels.direction in names(level.values)) {
    for (percent in c(FALSE, TRUE)) {
      for (direction in c("auto", "<", ">")) {
        context(sprintf("'roc' function works with percent = %s, marker = %s, levels.direction = %s, direction = %s", percent, marker, levels.direction, direction))
        expected.direction <- ifelse(direction == "auto", ifelse(levels.direction == "forward", "<", ">"), direction)
        r <- roc(aSAH$outcome, aSAH[[marker]], levels = level.values[[levels.direction]], direction = direction, percent = percent, quiet = TRUE)

        test_that("roc.formula produces the same results as roc.default", {
          rf <- roc(as.formula(sprintf("outcome ~ %s", marker)), data = aSAH, levels = level.values[[levels.direction]], direction = direction, percent = percent, quiet = TRUE)
          expect_is(rf, "roc")
          expect_equal(as.numeric(rf$auc), as.numeric(r$auc))
          for (item in c("percent", "sensitivities", "specificities", "thresholds", "direction", "cases", "controls")) {
            expect_identical(rf[[item]], r[[item]], label = sprintf("roc(outcome ~ %s, %s, %s, %s)[[\"%s\"]]", marker, levels.direction, percent, direction, item))
          }
          for (item in c("original.predictor", "original.response", "predictor", "response", "levels")) {
            expect_identical(unname(rf[[item]]), unname(r[[item]]), label = sprintf("roc(outcome ~ %s, %s, %s, %s)[[\"%s\"]]", marker, levels.direction, percent, direction, item))
          }
        })

        test_that("roc.default works with control/cases as well", {
          rcs <- roc(controls = r$controls, cases = r$cases, levels = level.values[[levels.direction]], direction = direction, percent = percent, quiet = TRUE)
          expect_is(rcs, "roc")
          expect_equal(as.numeric(rcs$auc), as.numeric(r$auc))
          for (item in c("percent", "sensitivities", "specificities", "thresholds", "direction", "cases", "controls")) {
            expect_identical(rcs[[item]], r[[item]])
          }
        })

        test_that("roc.default produces the expected results", {
          expect_is(r, "roc")
          expect_identical(r$percent, percent)
          expect_identical(r$direction, expected.direction)
          expect_identical(r$levels, level.values[[levels.direction]])

          expect_equal(r$thresholds, expected.roc[[marker]][[levels.direction]][[expected.direction]][["thresholds"]])
          expect_equal(r$sensitivities, expected.roc[[marker]][[levels.direction]][[expected.direction]][["sensitivities"]] * ifelse(percent, 100, 1))
          expect_equal(r$specificities, expected.roc[[marker]][[levels.direction]][[expected.direction]][["specificities"]] * ifelse(percent, 100, 1))
        })

        if (marker == "wfns") {
          available.smooth.methods <- "binormal"
        } else {
          available.smooth.methods <- smooth.methods
        }
        for (smooth.method in available.smooth.methods) {
          context(sprintf("smooth(roc(...)) works with percent = %s, marker = %s, levels.direction = %s, direction = %s and smooth.method = %s", percent, marker, levels.direction, direction, smooth.method))
          test_that("smoothing a ROC curve produces expected results", {
            skip_slow()
            skip_if(getRversion() < "4.4.0")
            if (smooth.method == "logcondens" || smooth.method == "logcondens.smooth") {
              testthat::skip_if_not_installed("logcondens")
            }
            if (smooth.method == "fitdistr") {
              testthat::skip_if_not_installed("MASS")
            }
            s <- smooth(r, method = smooth.method, 10)
            expect_is(s, "smooth.roc")
            expect_identical(s$percent, percent)
            expect_identical(s$direction, expected.direction)
            expect_equal(s$sensitivities, expected.roc[[marker]][[levels.direction]][[expected.direction]][["smooth"]][[smooth.method]][["sensitivities"]] * ifelse(percent, 100, 1))
            expect_equal(s$specificities, expected.roc[[marker]][[levels.direction]][[expected.direction]][["smooth"]][[smooth.method]][["specificities"]] * ifelse(percent, 100, 1))
          })
          test_that("building curve with smooth=TRUE produces expected results", {
            skip_slow()
            skip_if(getRversion() < "4.4.0")
            context(sprintf("roc(..., smooth=TRUE) works with percent = %s, marker = %s, levels.direction = %s, direction = %s and smooth.method = %s", percent, marker, levels.direction, direction, smooth.method))
            if (smooth.method == "logcondens" || smooth.method == "logcondens.smooth") {
              testthat::skip_if_not_installed("logcondens")
            }
            if (smooth.method == "fitdistr") {
              testthat::skip_if_not_installed("MASS")
            }
            s2 <- roc(aSAH$outcome, aSAH[[marker]],
              levels = level.values[[levels.direction]], direction = direction, percent = percent, quiet = TRUE,
              smooth = TRUE, smooth.n = 10, smooth.method = smooth.method
            )
            expect_is(s2, "smooth.roc")
            expect_identical(s2$percent, percent)
            expect_identical(s2$direction, expected.direction)
            expect_equal(s2$sensitivities, expected.roc[[marker]][[levels.direction]][[expected.direction]][["smooth"]][[smooth.method]][["sensitivities"]] * ifelse(percent, 100, 1))
            expect_equal(s2$specificities, expected.roc[[marker]][[levels.direction]][[expected.direction]][["smooth"]][[smooth.method]][["specificities"]] * ifelse(percent, 100, 1))
          })
        }
      }
    }
  }
}

# dump("expected.roc", file="helper-roc-expected.R")

test_that("roc.default handles NAs", {
  # Generate missing values
  aSAH.missing <- aSAH
  aSAH.missing$ndka[1:20] <- NA
  aSAH.missing$wfns[1:20] <- NA

  # na.rm=FALSE works
  # With NDKA
  expect_true(is.na(roc(aSAH.missing$outcome, aSAH.missing$ndka, na.rm = FALSE)))
  expect_false(is.na(auc(roc(aSAH.missing$outcome, aSAH.missing$ndka, na.rm = TRUE))))
  # With WFNS
  expect_true(is.na(roc(aSAH.missing$outcome, aSAH.missing$wfns, na.rm = FALSE)))
  expect_false(is.na(auc(roc(aSAH.missing$outcome, aSAH.missing$wfns, na.rm = TRUE))))

  # Same as subset
  expect_identical(
    as.numeric(auc(roc(aSAH.missing$outcome, aSAH.missing$ndka, na.rm = TRUE))),
    as.numeric(auc(roc(aSAH[21:113, ]$outcome, aSAH.missing[21:113, ]$ndka)))
  )
  # With ordered
  expect_identical(
    as.numeric(auc(roc(aSAH.missing$outcome, aSAH.missing$wfns, na.rm = TRUE))),
    as.numeric(auc(roc(aSAH[21:113, ]$outcome, aSAH.missing[21:113, ]$wfns)))
  )

  # Also with case/controls
  expect_identical(
    as.numeric(auc(roc(controls = aSAH.missing$ndka[aSAH.missing$outcome == "Good"], cases = aSAH.missing$ndka[aSAH.missing$outcome == "Poor"]))),
    as.numeric(auc(roc(aSAH[21:113, ]$outcome, aSAH.missing[21:113, ]$ndka)))
  )
  # With ordered
  expect_identical(
    as.numeric(auc(roc(controls = aSAH.missing$wfns[aSAH.missing$outcome == "Good"], cases = aSAH.missing$wfns[aSAH.missing$outcome == "Poor"]))),
    as.numeric(auc(roc(aSAH[21:113, ]$outcome, aSAH.missing[21:113, ]$wfns)))
  )
})

test_that("roc.formula behaves", {
  # By this point we've tested the main stuff, so just check a few basic elements
  roc.check.only.items <- c("sensitivities", "specificities", "thresholds", "cases", "controls")

  expect_identical(
    roc(outcome ~ wfns, data = aSAH)[roc.check.only.items],
    roc(aSAH$outcome, aSAH$wfns)[roc.check.only.items]
  )

  # formula without data
  expect_identical(
    roc(aSAH$outcome ~ aSAH$wfns)[roc.check.only.items],
    roc(aSAH$outcome, aSAH$wfns)[roc.check.only.items]
  )

  # formula with data from parent env
  outcome <- aSAH$outcome
  wfns <- aSAH$wfns
  expect_identical(
    roc(outcome ~ wfns)[roc.check.only.items],
    roc(outcome, wfns)[roc.check.only.items]
  )

  expect_identical(
    roc(outcome ~ wfns, data = aSAH, subset = (gender == "Female"))[roc.check.only.items],
    roc(aSAH$outcome[aSAH$gender == "Female"], aSAH$wfns[aSAH$gender == "Female"])[roc.check.only.items]
  )

  # Generate missing values
  aSAH.missing <- aSAH
  aSAH.missing$ndka[1:20] <- NA
  expect_warning(roc1 <- roc(outcome ~ ndka, data = aSAH.missing, na.action = na.omit), "na.omit")
  roc2 <- roc(aSAH[21:113, ]$outcome, aSAH[21:113, ]$ndka)[roc.check.only.items]
  expect_identical(
    roc1[roc.check.only.items],
    roc2[roc.check.only.items]
  )
  # na.fail should fail
  expect_error(roc(outcome ~ ndka, data = aSAH.missing, na.action = na.fail))
  # weights should fail too
  expect_error(roc(outcome ~ ndka, data = aSAH, weights = seq_len(nrow(aSAH)), quiet = TRUE), regexp = "weights are not supported")
  # invalid formula should fail
  expect_error(roc(~ndka, data = aSAH))


  # Both na.action and subset
  expect_warning(roc1 <- roc(outcome ~ ndka, data = aSAH.missing, na.action = na.omit, subset = (gender == "Female")), "na.omit")
  roc2 <- roc(aSAH[21:113, ]$outcome[aSAH[21:113, ]$gender == "Female"], aSAH[21:113, ]$ndka[aSAH[21:113, ]$gender == "Female"])
  expect_identical(
    roc1[roc.check.only.items],
    roc2[roc.check.only.items]
  )
})

test_that("roc can't take both response/predictor and case/control", {
  expect_error(roc(aSAH$outcome, aSAH$ndka, controls = aSAH$ndka[aSAH$outcome == "Good"], cases = aSAH$ndka[aSAH$outcome == "Poor"]))
})

test_that("roc with multiple predictors returns expected ROC curves", {
  roclist <- roc(outcome ~ wfns + ndka + s100b, data = aSAH, quiet = TRUE)
  expect_is(roclist, "list")
  expect_type(roclist, "list")
  expect_length(roclist, 3)
  expect_identical(names(roclist), c("wfns", "ndka", "s100b"))

  expect_equal_roc_formula(roclist$wfns, r.wfns)
  expect_equal_roc_formula(roclist$ndka, r.ndka)
  expect_equal_roc_formula(roclist$s100b, r.s100b)

  attach(aSAH)
  roclist <- roc(outcome ~ wfns + ndka + s100b, quiet = TRUE)
  expect_equal_roc_formula(roclist$wfns, r.wfns)
  expect_equal_roc_formula(roclist$ndka, r.ndka)
  expect_equal_roc_formula(roclist$s100b, r.s100b)
  detach(aSAH)
})

test_that("extra arguments passed to roc with multiple predictors", {
  roclist <- roc(outcome ~ wfns + ndka + s100b,
    data = aSAH, quiet = TRUE,
    percent = TRUE, partial.auc = c(90, 99)
  )

  expect_equal_roc_formula(roclist$wfns, r.wfns.percent.partial1)
  expect_equal_roc_formula(roclist$ndka, r.ndka.percent.partial1)
  expect_equal_roc_formula(roclist$s100b, r.s100b.percent.partial1)
})

test_that("roc works with densitites", {
  range.ndka <- range(aSAH$ndka)
  bw <- bw.nrd0(aSAH$ndka)
  from <- min(aSAH$ndka) - (3 * bw)
  to <- max(aSAH$ndka) + (3 * bw)
  density.controls <- density(aSAH$ndka[aSAH$outcome == "Good"], from = from, to = to, bw = bw)
  density.cases <- density(aSAH$ndka[aSAH$outcome == "Poor"], from = from, to = to, bw = bw)
  density.roc <- roc(density.cases = density.cases$y, density.controls = density.controls$y)
  smoothed.roc <- smooth(r.ndka, method = "density")

  expect_is(density.roc, "smooth.roc")
  expect_equal(density.roc$sensitivities, smoothed.roc$sensitivities)
  expect_equal(density.roc$specificities, smoothed.roc$specificities)
  expect_equal(as.numeric(density.roc$auc), as.numeric(smoothed.roc$auc))
})

test_that("roc.density works with extra arguments", {
  skip_if(getRversion() < "4.4.0")
  range.ndka <- range(aSAH$ndka)
  bw <- bw.nrd0(aSAH$ndka)
  from <- min(aSAH$ndka) - (3 * bw)
  to <- max(aSAH$ndka) + (3 * bw)
  density.controls <- density(aSAH$ndka[aSAH$outcome == "Good"], from = from, to = to, bw = bw)
  density.cases <- density(aSAH$ndka[aSAH$outcome == "Poor"], from = from, to = to, bw = bw)

  density.roc.partial <- roc(
    density.cases = density.cases$y,
    density.controls = density.controls$y,
    partial.auc = c(1, .9), partial.auc.focus = "se",
    partial.auc.correct = TRUE
  )
  expect_equal(as.numeric(density.roc.partial$auc), 0.506196226093)

  density.roc.percent <- roc(
    density.cases = density.cases$y,
    density.controls = density.controls$y,
    percent = TRUE
  )
  expect_equal(as.numeric(density.roc.percent$auc), 60.44728041)
})

test_that("roc doesn't accept density with other arguments", {
  density.controls <- density(aSAH$ndka[aSAH$outcome == "Good"])
  density.cases <- density(aSAH$ndka[aSAH$outcome == "Poor"])
  expect_error(
    roc(aSAH$outcome, aSAH$ndka, density.cases = density.controls, density.controls = density.cases),
    "incompatible"
  )
  expect_error(
    roc(cases = aSAH$ndka, controls = aSAH$ndka, density.cases = density.controls, density.controls = density.cases),
    "incompatible"
  )
})

test_that("roc.data.frame works", {
  r <- roc(aSAH, outcome, s100b, ret = "roc")
  expect_is(r, "roc")
  co <- roc(aSAH, outcome, s100b, ret = "coords")
  expect_equal(dim(co), c(51, 3))
  co <- roc(aSAH, outcome, s100b, ret = "all_coords")
  expect_equal(nrow(co), 51)
  expect_true(nrow(co) >= 22)
})

test_that("roc.data.frame works with quoted names", {
  r <- roc(aSAH, "outcome", "s100b", ret = "roc")
  expect_is(r, "roc")
  co <- roc(aSAH, "outcome", "s100b", ret = "coords")
  expect_equal(dim(co), c(51, 3))
  co <- roc(aSAH, "outcome", "s100b", ret = "all_coords")
  expect_equal(nrow(co), 51)
  expect_true(nrow(co) >= 22)
})

test_that("roc_ works", {
  r <- roc_(aSAH, "outcome", "s100b", ret = "roc")
  expect_is(r, "roc")
  co <- roc_(aSAH, "outcome", "s100b", ret = "coords")
  expect_equal(dim(co), c(51, 3))
  co <- roc_(aSAH, "outcome", "s100b", ret = "all_coords")
  expect_equal(nrow(co), 51)
  expect_true(nrow(co) >= 22)
})


test_that("roc.data.frame reject invalid columns", {
  outcomes <- aSAH$outcome
  suppressWarnings(expect_error(roc(aSAH, outcomes, s100b), "Column")) # Warning about NSE
  suppressWarnings(expect_error(roc(aSAH, "outcomes", "s100b"), "Column")) # Warning about NSE
  expect_error(roc_(aSAH, "outcomes", "s100b"), "Column")
  s100c <- aSAH$s100b
  suppressWarnings(expect_error(roc(aSAH, outcome, s100c), "Column")) # Warning about NSE
  suppressWarnings(expect_error(roc(aSAH, "outcome", "s100c"), "Column")) # Warning about NSE
  expect_error(roc_(aSAH, "outcome", "s100c"), "Column")
})

test_that("roc reject and warns for invalid levels", {
  expect_warning(expect_error(roc(aSAH$gos6, aSAH$s100b), "No case observation"), "levels")
  expect_error(roc(aSAH$gos6, aSAH$s100b, levels = levels(aSAH$gos6)), "levels")
  expect_warning(roc(factor(aSAH$gos6), aSAH$s100b, quiet = TRUE), "levels")

  expect_warning(expect_error(roc(aSAH, gos6, s100b), "No case observation"), "levels")
  expect_error(roc(aSAH, gos6, s100b, levels = levels(aSAH$gos6)), "levels")
  dat <- aSAH
  dat$gos6 <- factor(aSAH$gos6)
  expect_warning(roc(dat, gos6, s100b, quiet = TRUE), "levels")
})

test_that("roc reject and warns for invalid predictors", {
  expect_error(roc(aSAH$outcome, as.character(aSAH$wfns)), "Predictor")
  expect_warning(roc(aSAH$outcome, as.matrix(aSAH$ndka)), "Deprecated")
  expect_warning(roc(as.matrix(aSAH$outcome), aSAH$ndka), "Deprecated")

  expect_error(roc(aSAH$outcome[1:100], aSAH$ndka), "length")
  expect_error(roc(aSAH$outcome[1:100], aSAH$ndka[1:50]), "length")
})

test_that("roc reject requires cases & controls", {
  expect_error(roc(aSAH$outcome[aSAH$outcome == "Good"], aSAH$ndka[aSAH$outcome == "Good"]), "case")
  expect_error(roc(aSAH$outcome[aSAH$outcome == "Poor"], aSAH$ndka[aSAH$outcome == "Poor"]), "control")

  expect_error(roc(aSAH[aSAH$outcome == "Good", ], outcome, ndka), "case")
  expect_error(roc(aSAH[aSAH$outcome == "Poor", ], outcome, ndka), "control")
})


test_that("roc works with ordered predictor", {
  wfns2 <- ordered(as.numeric(aSAH$wfns) + 2)
  r <- roc(aSAH$outcome, wfns2)
  expect_equal(r$thresholds, c(-Inf, 3.5, 4.5, 5.5, 6.5, Inf))

  levels(wfns2) <- letters[1:5]
  # TODO: this behavior should be fixed, see issue #63.
  # For now ensure basic behavior with warning is at least consistent.
  expect_warning(r <- roc(aSAH$outcome, wfns2))
  expect_equal(r$thresholds, c(-Inf, 1.5, 2.5, 3.5, 4.5, Inf))
  # In reality we want to say something like  c(-Inf, "a", "b", "c", "d", "e", Inf)
})


test_that("roc works with spaces in formula", {
  test_dat <- aSAH[c("outcome", "wfns", "age", "ndka")]
  colnames(test_dat)[3] <- "A ge"

  # All data
  r <- roc(outcome ~ ., data = test_dat, quiet = TRUE)
  expect_equal(length(r), 3)
  expect_equal(names(r), c("wfns", "A ge", "ndka"))
  expect_equal(r$`A ge`$predictor, aSAH$age)

  # One column
  r <- roc(outcome ~ `A ge`, data = test_dat, quiet = TRUE)
  expect_equal(r$predictor, aSAH$age)

  # Two columns
  r <- roc(outcome ~ `A ge` + ndka, data = test_dat, quiet = TRUE)
  expect_equal(names(r), c("A ge", "ndka"))
  expect_equal(r$`A ge`$predictor, aSAH$age)

  # Two columns, different order
  r <- roc(outcome ~ ndka + `A ge`, data = test_dat, quiet = TRUE)
  expect_equal(names(r), c("ndka", "A ge"))
  expect_equal(r$`A ge`$predictor, aSAH$age)
})


test_that("roc works with `with` and formula", {
  # This is the correct way to do things... with on the formula
  formW <- with(aSAH, as.formula("outcome ~ wfns"))
  roc(formW)

  # This must fail gracefully
  form <- as.formula("outcome ~ wfns")
  expect_error(roc(form), "outcome")

  # Now this can work with #111
  with(aSAH, roc(form))

  # Wrapping functions can mess things up...
  wrap_roc <- function(formula) {
    roc(formula)
  }
  with(aSAH, wrap_roc(form))

  wrap_roc2 <- function(formula) {
    with(aSAH, roc(formula))
  }
  wrap_roc2(form)
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
# 			for (smooth.method in available.smooth.methods) {
# 				s <- smooth(r, method=smooth.method, 10)
# 				expected.roc[[marker]][[levels.direction]][[expected.direction]][["smooth"]][[smooth.method]][["sensitivities"]] <- s$sensitivities
# 				expected.roc[[marker]][[levels.direction]][[expected.direction]][["smooth"]][[smooth.method]][["specificities"]] <- s$specificities
# 			}
# 		}
# 	}
# }
# save("expected.roc", system.file("extdata", "test-roc-expected.R", package="pROC"), file = "dump_roc_expected.R")
