library(pROC)
data(aSAH)

context("ordered predictors")

text.levels <- aSAH$wfns
levels(text.levels) <- c("very low", "low", "medium", "high", "very high")

test_that("issue 63: text ordered levels are kept as thresholds", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  expect_true(is.ordered(r$predictor))
  expect_true(is.ordered(r$cases))
  expect_true(is.ordered(r$controls))
  expect_true(is.ordered(r$thresholds))
  expect_equal(r$direction, "<")
  expect_equal(as.character(r$thresholds), c("very low", "low", "medium", "high", "very high", "Inf"))
  expect_equal(levels(r$predictor), c("very low", "low", "medium", "high", "very high"))
  expect_false("Inf" %in% levels(r$predictor))
})

test_that("auto direction on ordered uses level ranks, not a converted predictor", {
  # even-length groups: median() does not work on ordered values
  pred <- ordered(c("a", "a", "b", "b", "c", "c"), levels = c("a", "b", "c"))
  resp <- factor(c("ctl", "ctl", "ctl", "case", "case", "case"), levels = c("ctl", "case"))
  expect_error(median(pred[resp == "ctl"]), "numeric")
  r <- roc(resp, pred, quiet = TRUE)
  expect_equal(r$direction, "<")
  expect_true(is.ordered(r$cases))
  expect_true(is.ordered(r$controls))
})

test_that("an observation at the threshold is classified as positive", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE, direction = "<")
  co <- coords(r, "medium", input = "threshold")
  expect_equal(
    co$sensitivity,
    mean(r$cases >= ordered("medium", levels = levels(text.levels)))
  )
  expect_equal(
    co$specificity,
    mean(r$controls < ordered("medium", levels = levels(text.levels)))
  )
  # Inclusive: medium itself is positive
  expect_true(mean(r$cases >= ordered("medium", levels = levels(text.levels))) >=
    mean(r$cases > ordered("medium", levels = levels(text.levels))))
})

test_that("direction > uses <= at the threshold", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE, direction = ">")
  co <- coords(r, "medium", input = "threshold")
  expect_equal(
    co$sensitivity,
    mean(r$cases <= ordered("medium", levels = levels(text.levels)))
  )
  expect_equal(
    co$specificity,
    mean(r$controls > ordered("medium", levels = levels(text.levels)))
  )
})

test_that("unused defined levels are kept as thresholds", {
  keep <- text.levels != "medium" | is.na(text.levels)
  pred <- text.levels[keep]
  resp <- aSAH$outcome[keep]
  expect_true("medium" %in% levels(pred))
  expect_false("medium" %in% as.character(pred))
  r <- roc(resp, pred, quiet = TRUE, direction = "<")
  expect_true("medium" %in% as.character(r$thresholds))
  co_med <- coords(r, "medium", input = "threshold")
  co_high <- coords(r, "high", input = "threshold")
  expect_equal(co_med$sensitivity, co_high$sensitivity)
  expect_equal(co_med$specificity, co_high$specificity)
})

test_that("coords returns an ordered threshold column", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  co <- coords(r, "all")
  expect_true(is.data.frame(co))
  expect_true(is.ordered(co$threshold))
  expect_equal(levels(co$threshold), levels(r$thresholds))
  co_one <- coords(r, "low", input = "threshold")
  expect_true(is.ordered(co_one$threshold))
  expect_equal(levels(co_one$threshold), levels(r$thresholds))
  expect_equal(as.character(co_one$threshold), "low")
  co_ord <- coords(r, ordered("low", levels = levels(text.levels)), input = "threshold")
  expect_true(is.ordered(co_ord$threshold))
  expect_equal(levels(co_ord$threshold), levels(r$thresholds))
  expect_equal(as.character(co_ord$threshold), "low")
})

test_that("numeric x is rejected for ordered threshold lookup", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  expect_error(coords(r, 2, input = "threshold"), "Numeric")
  expect_error(ci.thresholds(r, thresholds = 2, boot.n = 3), "Numeric")
})

test_that("cases and controls require identical levels", {
  controls <- aSAH$wfns[aSAH$outcome == "Good"]
  cases <- aSAH$wfns[aSAH$outcome == "Poor"]
  r <- roc(controls = controls, cases = cases, quiet = TRUE)
  expect_equal(as.character(r$predictor), c(as.character(controls), as.character(cases)))
  expect_error(
    roc(controls = controls, cases = as.numeric(cases), quiet = TRUE),
    "both be numeric or both be ordered"
  )

  cases_reordered <- ordered(as.character(cases), levels = rev(levels(cases)))
  expect_error(
    roc(controls = controls, cases = cases_reordered, quiet = TRUE),
    "Levels of cases and controls differ"
  )

  cases_extra <- ordered(as.character(cases), levels = c(levels(cases), "extra"))
  expect_error(
    roc(controls = controls, cases = cases_extra, quiet = TRUE),
    "Levels of cases and controls differ"
  )

  ctrl_dropped <- droplevels(controls[as.integer(controls) <= 2])
  expect_true(length(base::levels(ctrl_dropped)) < length(base::levels(cases)))
  expect_error(
    roc(controls = ctrl_dropped, cases = cases, quiet = TRUE),
    "Levels of cases and controls differ"
  )
})

test_that("roc_utils_c_ordered does not drop to integer codes", {
  a <- aSAH$wfns[1:5]
  b <- aSAH$wfns[6:10]
  combined <- pROC:::roc_utils_c_ordered(a, b)
  expect_true(is.ordered(combined))
  expect_equal(levels(combined), levels(aSAH$wfns))
  expect_equal(as.character(combined), as.character(c(as.character(a), as.character(b))))
})

test_that("ci.thresholds accepts a level label", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  ci <- ci.thresholds(r, thresholds = "medium", boot.n = 5)
  expect_s3_class(ci, "ci.thresholds")
  expect_equal(rownames(ci$sensitivity), "medium")
})

test_that("plot print.thres does not error with ordered thresholds", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plot(r, print.thres = "best"), NA)
  expect_error(plot(r, print.thres = "medium"), NA)
})

test_that("coords accepts Inf sentinel labels", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE, direction = "<")
  co <- coords(r, "Inf", input = "threshold")
  expect_equal(as.character(co$threshold), "Inf")
  expect_equal(co$specificity, 1)
  expect_equal(co$sensitivity, 0)

  rgt <- roc(aSAH$outcome, text.levels, quiet = TRUE, direction = ">")
  co_ninf <- coords(rgt, "-Inf", input = "threshold")
  expect_equal(as.character(co_ninf$threshold), "-Inf")
  expect_equal(co_ninf$specificity, 1)
  expect_equal(co_ninf$sensitivity, 0)
})

test_that("ci.coords returns NA for ordered threshold CIs", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  expect_warning(
    ci <- ci.coords(r, x = "medium", input = "threshold",
                    ret = c("threshold", "specificity", "sensitivity"), boot.n = 5),
    "not available for ordered"
  )
  expect_true(all(is.na(ci$threshold)))
  expect_true(all(is.finite(ci$specificity)))
  expect_true(all(is.finite(ci$sensitivity)))
})

test_that("only binormal smoothing is allowed for ordered predictors", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  expect_error(smooth(r, method = "density"), "binormal")
  expect_error(smooth(r, method = "fitdistr"), "binormal")
  expect_s3_class(smooth(r, method = "binormal"), "smooth.roc")
})

test_that("wfns SE/SP match as.numeric(wfns) after the label shift", {
  r_ord <- r.wfns
  r_num <- r.wfns.num
  expect_equal(as.numeric(r_ord$auc), as.numeric(r_num$auc))
  expect_equal(r_ord$sensitivities, r_num$sensitivities)
  expect_equal(r_ord$specificities, r_num$specificities)
  expect_equal(as.character(r_ord$thresholds), c("1", "2", "3", "4", "5", "Inf"))
  expect_equal(r_num$thresholds, c(-Inf, 1.5, 2.5, 3.5, 4.5, Inf))
})

test_that("a level named Inf is ranked by position, not by label", {
  resp <- factor(c("ctl", "ctl", "ctl", "case", "case", "case"), levels = c("ctl", "case"))
  vals <- c("high", "high", "mid", "low", "low", "mid")
  ok <- ordered(vals, levels = c("low", "mid", "high"))
  inf_mid <- ordered(sub("^mid$", "Inf", vals), levels = c("low", "Inf", "high"))
  r_ok <- roc(resp, ok, quiet = TRUE, direction = ">")
  r_inf <- roc(resp, inf_mid, quiet = TRUE, direction = ">")
  expect_equal(as.numeric(r_inf$auc), as.numeric(r_ok$auc))
  expect_equal(r_inf$sensitivities, r_ok$sensitivities)
  expect_equal(r_inf$specificities, r_ok$specificities)
  # Direction "<" needs the Inf sentinel, so the same level name is rejected.
  expect_error(roc(resp, inf_mid, quiet = TRUE, direction = "<"), "collides")
})

test_that("a colliding Inf sentinel is rejected for the matching direction", {
  resp <- factor(c("ctl", "ctl", "ctl", "case", "case", "case"), levels = c("ctl", "case"))
  expect_error(
    roc(resp, ordered(c("a", "a", "a", "Inf", "Inf", "Inf"), levels = c("a", "Inf")),
        direction = "<", quiet = TRUE),
    "collides"
  )
  expect_error(
    roc(resp, ordered(c("-Inf", "-Inf", "-Inf", "z", "z", "z"), levels = c("-Inf", "z")),
        direction = ">", quiet = TRUE),
    "collides"
  )
})

test_that("keyword prefixes do not hijack level labels", {
  lv <- c("a", "b", "c")
  resp <- factor(rep(c("ctl", "case"), each = 5), levels = c("ctl", "case"))
  pred <- ordered(c("a", "a", "b", "b", "c", "b", "c", "c", "c", "b"), levels = lv)
  r <- roc(resp, pred, quiet = TRUE)
  co_a <- coords(r, "a", input = "threshold")
  expect_equal(as.character(co_a$threshold), "a")
  expect_equal(nrow(co_a), 1L)
  co_b <- coords(r, "b", input = "threshold")
  expect_equal(as.character(co_b$threshold), "b")
  expect_equal(nrow(co_b), 1L)
  co_all <- coords(r, "all")
  expect_equal(nrow(co_all), length(r$thresholds))
  co_bes <- coords(r, "bes")
  expect_equal(nrow(co_bes), nrow(coords(r, "best")))
})

test_that("an exact keyword that is also a level warns and keeps the keyword", {
  lv <- c("low", "best", "high")
  resp <- factor(rep(c("ctl", "case"), each = 4), levels = c("ctl", "case"))
  pred <- ordered(c("low", "low", "best", "high", "best", "high", "high", "best"), levels = lv)
  r <- roc(resp, pred, quiet = TRUE)
  expect_warning(co <- coords(r, "best"), "both a keyword and a threshold level")
  expect_equal(nrow(co), nrow(suppressWarnings(coords(r, "bes"))))
  co_level <- coords(r, ordered("best", levels = lv), input = "threshold")
  expect_equal(as.character(co_level$threshold), "best")
  expect_equal(nrow(co_level), 1L)
})

test_that("plot print.thres accepts a vector of labels and an ordered value", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plot(r, print.thres = c("low", "high")), NA)
  expect_error(plot(r, print.thres = ordered("medium", levels = levels(text.levels))), NA)
})

test_that("roc(smooth=TRUE) with a non-binormal method is rejected for ordered predictors", {
  expect_error(
    roc(aSAH$outcome, text.levels, quiet = TRUE, smooth = TRUE, smooth.method = "density"),
    "binormal"
  )
  expect_error(
    roc(aSAH$outcome, text.levels, quiet = TRUE, smooth = TRUE, smooth.method = "logcondens"),
    "binormal"
  )
  expect_s3_class(
    roc(aSAH$outcome, text.levels, quiet = TRUE, smooth = TRUE, smooth.method = "binormal"),
    "smooth.roc"
  )
})

test_that("transpose with numeric-only ret stays numeric on an ordered ROC", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  suppressWarnings({
    obtained <- coords(r, "all", ret = c("specificity", "sensitivity"), transpose = TRUE)
    obtained_m <- coords(r, "all", ret = c("specificity", "sensitivity"), as.matrix = TRUE)
  })
  expect_type(obtained, "double")
  expect_type(obtained_m, "double")
})

test_that("ci.coords default ret omits threshold on an ordered ROC", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE)
  ci <- ci.coords(r, "best", boot.n = 5, best.policy = "random")
  expect_false("threshold" %in% names(ci))
  expect_true("sensitivity" %in% names(ci))
  expect_true("specificity" %in% names(ci))
})

test_that("unknown and missing threshold labels have clear messages", {
  r <- roc(aSAH$outcome, text.levels, quiet = TRUE, direction = "<")
  expect_error(coords(r, "-Inf", input = "threshold"), "Inf\" sentinel")
  expect_error(coords(r, c("low", NA), input = "threshold"), "Missing values")
})

test_that("coords on a numeric ROC restores the match.arg error", {
  expect_error(coords(r.ndka, "foo"), "arg")
})
