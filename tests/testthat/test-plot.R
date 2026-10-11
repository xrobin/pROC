context("plot")

# Tests powered by vdiffr.
# To update the reference with vdiffr:
# > library(vdiffr)
# > source("tests/testthat.R")
# > manage_cases()

test_that("plot draws correctly", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")
  test_basic_plot <- function() plot(r)
  # S100b
  r <- r.s100b
  expect_doppelganger("basic-s100b", test_basic_plot)

  r <- r.ndka
  expect_doppelganger("basic-ndka", test_basic_plot)

  r <- r.wfns
  expect_doppelganger("basic-wfns", test_basic_plot)
})

test_that("legacy.axis works correctly", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")
  r <- r.s100b
  test_legacy.axis_plot <- function() plot(r, legacy.axes = TRUE)
  expect_doppelganger("legacy.axes", test_legacy.axis_plot)
})

test_that("Advanced screenshot 1 works correctly", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")
  test_advanced_screenshot_1 <- function() {
    # Two partial-AUC regions (focus sp and focus se) on the same curve, built
    # from separate auc() objects and drawn via panel.first instead of the old
    # two-call plot.roc(..., add=TRUE, type="n") idiom. This changes the
    # z-order versus the old idiom: the old two calls drew the identity line
    # above the SP polygon but below the SE polygon (inconsistently, since the
    # second call's output always lands on top of the first call's entire
    # output, curve included); panel.first draws all polygons once, under a
    # single curve draw, so the identity line now sits consistently above
    # both. Reviewed side-by-side via vdiffr::snapshot_review() and accepted
    # as an improvement.
    auc_sp <- auc(r.s100b.percent, partial.auc = c(100, 90), partial.auc.correct = TRUE)
    auc_se <- auc(r.s100b.percent,
      partial.auc = c(100, 90), partial.auc.correct = TRUE,
      partial.auc.focus = "se"
    )
    plot(r.s100b.percent,
      main = "Partial AUC (pAUC)",
      panel.first = {
        polygon_max_auc(auc_se, col = "#00860022")
        polygon_max_auc(auc_sp, col = "#1c61b622")
        polygon_auc(auc_se, col = "#008600")
        polygon_auc(auc_sp, col = "#1c61b6")
      }
    )
    text(auc_sp,
      xy = c(50, 50), col = "#1c61b6",
      pattern = "Corrected pAUC (100-90%% SP):\n%.1f%%"
    )
    text(auc_se,
      xy = c(50, 40), col = "#008600",
      pattern = "Corrected pAUC (100-90%% SE):\n%.1f%%"
    )
  }
  expect_doppelganger("advanced.screenshot.1", test_advanced_screenshot_1)
})


test_that("Advanced screenshot 2 works correctly", {
  skip_if_not_installed("vdiffr")
  skip_slow()
  skip_if(getRversion() < "4.1")
  test_advanced_screenshot_2 <- function() {
    RNGkind(sample.kind = "Rejection")
    set.seed(42) # For reproducible CI
    suppressMessages(rocobj <- plot.roc(aSAH$outcome, aSAH$s100b,
      main = "Confidence intervals", percent = TRUE,
      ci = TRUE # compute AUC (of AUC by default)
    ))
    text(rocobj) # print the AUC (will contain the CI)

    ciobj <- ci.se(rocobj, # CI of sensitivity
      specificities = seq(0, 100, 5)
    ) # over a select set of specificities

    plot(ciobj, type = "shape", col = "#1c61b6AA") # plot as a blue shape
    plot(ci(rocobj, of = "thresholds", thresholds = "best")) # add one threshold
  }
  expect_doppelganger("advanced.screenshot.2", test_advanced_screenshot_2)
})


test_that("Advanced screenshot 3 works correctly", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.4.0")
  test_advanced_screenshot_3 <- function() {
    plot(r.s100b.percent, main = "Smoothing")

    lines(smooth(r.s100b.percent), # smoothing (default: binormal)
      col = "#1c61b6"
    )

    lines(smooth(r.s100b.percent, method = "density"), # density smoothing
      col = "#008600"
    )

    lines(
      smooth(r.s100b.percent,
        method = "fitdistr", # fit a distribution
        density = "lognormal"
      ), # let the distribution be log-normal
      col = "#840000"
    )

    legend("bottomright", legend = c("Empirical", "Binormal", "Density", "Fitdistr\n(Log-normal)"), col = c("black", "#1c61b6", "#008600", "#840000"), lwd = 2)
  }
  expect_doppelganger("advanced.screenshot.3", test_advanced_screenshot_3)
})



test_that("Advanced screenshot 4 works correctly", {
  skip_if_not_installed("vdiffr")
  skip_slow()
  skip_if(getRversion() < "4.1")
  test_advanced_screenshot_4 <- function() {
    RNGkind(sample.kind = "Rejection")
    set.seed(42) # For reproducible CI
    suppressMessages(rocobj <- plot.roc(aSAH$outcome, aSAH$s100b,
      main = "Confidence intervals of specificity/sensitivity", percent = TRUE,
      ci = TRUE, of = "se", # ci of sensitivity
      specificities = seq(0, 100, 5), # on a select set of specificities
      ci.type = "shape", ci.col = "#1c61b6AA" # plot the CI as a blue shape
    ))

    plot(ci.sp(rocobj, sensitivities = seq(0, 100, 5)), # ci of specificity
      type = "bars"
    ) # print this one as bars
  }
  expect_doppelganger("advanced.screenshot.4", test_advanced_screenshot_4)
})



test_that("Advanced screenshot 5 works correctly", {
  skip_if_not_installed("vdiffr")
  skip_slow()
  skip_if(getRversion() < "4.1")
  test_advanced_screenshot_5 <- function() {
    RNGkind(sample.kind = "Rejection")
    set.seed(42) # For reproducible CI
    suppressMessages(plot.roc(aSAH$outcome, aSAH$s100b,
      main = "Confidence interval of a threshold", percent = TRUE,
      ci = TRUE, of = "thresholds", # compute AUC (of threshold)
      thresholds = "best", # select the (best) threshold
      print.thres = "best" # also highlight this threshold on the plot
    ))
  }
  expect_doppelganger("advanced.screenshot.5", test_advanced_screenshot_5)
})


test_that("Advanced screenshot 6 works correctly", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")
  test_advanced_screenshot_6 <- function() {
    plot(r.s100b.percent, main = "Statistical comparison", col = "#1c61b6")
    lines(r.ndka.percent, col = "#008600")
    testobj <- roc.test(r.s100b.percent, r.ndka.percent)
    text(50, 50, labels = paste("p-value =", format.pval(testobj$p.value)), adj = c(0, .5))
    legend("bottomright", legend = c("S100B", "NDKA"), col = c("#1c61b6", "#008600"), lwd = 2)
  }
  expect_doppelganger("advanced.screenshot.6", test_advanced_screenshot_6)
})


test_that("plot and lines work with formula and subset", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")
  test_plot_formula <- function() {
    suppressMessages(plot.roc(outcome ~ ndka, data = aSAH, subset = gender == "Female", col = "red"))
    suppressMessages(lines.roc(outcome ~ ndka, data = aSAH))
    suppressMessages(lines.roc(outcome ~ ndka, data = aSAH, subset = gender == "Male", col = "blue"))
  }
  expect_doppelganger("plot_formula", test_plot_formula)
})


test_that("PR curve with CI works", {
  skip_if_not_installed("vdiffr")
  skip_slow()
  skip_if(getRversion() < "4.1")
  test_plot_pr <- function() {
    RNGkind(sample.kind = "Rejection")
    set.seed(42) # For reproducible CI

    co <- coords(r.s100b, x = "all", input = "recall", ret = c("recall", "precision"))
    ci <- ci.coords(r.s100b, x = seq(0, 1, .1), input = "recall", ret = "precision")
    plot(co, type = "l", ylim = c(0, 1))
    testthat::expect_warning(plot(ci, type = "shape"), "Low definition shape")
    plot(ci, type = "bars")
    lines(co)
  }
  expect_doppelganger("plot_pr", test_plot_pr)
})

test_that("plot.roc and lines.roc formula methods keep the variable names", {
  pdf(NULL)
  on.exit(dev.off())
  dataline <- function(x) grep("^Data:", capture.output(print(x)), value = TRUE)
  expected <- dataline(roc(outcome ~ s100b, data = aSAH, quiet = TRUE))
  expect_identical(expected, "Data: s100b in 72 controls (outcome Good) < 41 cases (outcome Poor).")
  expect_identical(dataline(plot.roc(outcome ~ s100b, data = aSAH, quiet = TRUE)), expected)
  expect_identical(dataline(lines.roc(outcome ~ s100b, data = aSAH, quiet = TRUE)), expected)
})

test_that("plot.ci.coords draws visible ticks for a single x", {
  skip_if(getRversion() < "4.1")
  set.seed(42)
  ci <- ci.coords(r.s100b, x = 0.5, input = "specificity", ret = "sensitivity", boot.n = 10, progress = "none")
  expect_doppelganger("plot-ci-coords-single-x", function() {
    plot(r.s100b)
    plot(ci)
  })
})

test_that("plot.ci.coords gives a clear error for non-numeric x", {
  set.seed(42)
  pdf(NULL)
  on.exit(dev.off())
  plot(r.s100b)
  ci <- ci.coords(r.s100b, x = "best", ret = "sensitivity", best.policy = "random", boot.n = 10, progress = "none")
  expect_error(plot(ci), "numeric 'x'")
})

test_that("plot with ci = TRUE on a curve without CI gives a clear error", {
  set.seed(42)
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plot(r.s100b, ci = TRUE), "has no confidence interval")
  expect_error(plot(smooth(r.s100b), ci = TRUE), "has no confidence interval")
})

test_that("lines.roc default/formula return the NA/NaN of roc() with a warning, like plot.roc", {
  pdf(NULL)
  on.exit(dev.off())
  plot(r.s100b)
  d <- data.frame(
    y = c(0, 0, 0, 1, 1, 1),
    xna = c(1, 2, NA, 4, 5, 6),
    xinf = c(1, 2, Inf, 4, 5, 6)
  )
  expect_warning(r <- withVisible(lines.roc(d$y, d$xna, na.rm = FALSE, quiet = TRUE)), "Nothing plotted")
  expect_identical(r, list(value = NA, visible = FALSE))
  expect_warning(r <- withVisible(lines.roc(y ~ xna, d, na.rm = FALSE, quiet = TRUE)), "Nothing plotted")
  expect_identical(r, list(value = NA, visible = FALSE))
  expect_warning(r <- withVisible(lines.roc(d$y, d$xinf, quiet = TRUE)), "Nothing plotted")
  expect_identical(r, list(value = NaN, visible = FALSE))
  expect_warning(r <- withVisible(lines.roc(y ~ xinf, d, quiet = TRUE)), "Nothing plotted")
  expect_identical(r, list(value = NaN, visible = FALSE))
  expect_warning(r <- withVisible(plot.roc(d$y, d$xinf, quiet = TRUE)), "Nothing plotted")
  expect_identical(r, list(value = NaN, visible = FALSE))
})

test_that("percent print.thres labels keep 3 decimals for the threshold", {
  skip_if(getRversion() < "4.1")
  r.percent <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE, quiet = TRUE)
  expect_doppelganger("plot-percent-print-thres-best", function() plot(r.percent, print.thres = "best"))
  expect_doppelganger("plot-percent-print-thres-local-maximas", function() plot(r.percent, print.thres = "local maximas"))
})

test_that("CI shapes are drawn in increasing order whatever order the CI was computed in", {
  skip_if(getRversion() < "4.1")
  set.seed(42)
  ci_sp <- ci.sp(r.s100b, sensitivities = c(0.75, 0.25, 0.5), boot.n = 10, progress = "none")
  ci_se <- ci.se(r.s100b, specificities = c(0.75, 0.25, 0.5), boot.n = 10, progress = "none")
  ci_co <- ci.coords(r.s100b, x = c(0.75, 0.25, 0.5), input = "specificity", ret = "sensitivity", boot.n = 10, progress = "none")
  expect_doppelganger("plot-ci-sp-shape-unordered-sensitivities", function() {
    plot(r.s100b)
    expect_warning(plot(ci_sp, type = "shape", no.roc = TRUE), "Low definition shape")
  })
  expect_doppelganger("plot-ci-se-shape-unordered-specificities", function() {
    plot(r.s100b)
    expect_warning(plot(ci_se, type = "shape", no.roc = TRUE), "Low definition shape")
  })
  expect_doppelganger("plot-ci-coords-shape-unordered-x", function() {
    plot(r.s100b)
    expect_warning(plot(ci_co, type = "shape"), "Low definition shape")
  })
})

test_that("print.auc does not print the CI of another AUC specification", {
  skip_if(getRversion() < "4.1")
  r.ci <- roc(aSAH$outcome, aSAH$s100b, ci = TRUE, quiet = TRUE)
  expect_doppelganger("plot-print-auc-partial-auc-without-stored-ci", function() {
    plot(r.ci, print.auc = TRUE, reuse.auc = FALSE, partial.auc = c(1, 0.8))
  })
  # Same specification: the stored CI is still printed
  expect_doppelganger("plot-print-auc-reuse-auc-false-with-stored-ci", function() {
    plot(r.ci, print.auc = TRUE, reuse.auc = FALSE)
  })
  expect_doppelganger("plot-print-auc-with-stored-ci", function() plot(r.ci, print.auc = TRUE))
})

test_that("legacy.axes labels are not padded", {
  skip_if(getRversion() < "4.1")
  r.percent <- roc(c(0, 0, 1, 1), c(1, 2, 3, 4), percent = TRUE, quiet = TRUE)
  expect_doppelganger("plot-legacy-axes-percent", function() plot(r.percent, legacy.axes = TRUE))
  # Same number of decimals on all labels: see the "legacy.axes" snapshot
})

test_that("plot.roc draws the sensitivity of a ci.coords with plot.ci.coords", {
  skip_if(getRversion() < "4.1")
  set.seed(42)
  r <- r.s100b
  r$ci <- ci.coords(r.s100b, x = c(0.5, 0.9), input = "specificity", ret = "sensitivity", boot.n = 10, progress = "none")
  expect_doppelganger("plot-roc-ci-coords-sensitivity", function() plot(r))
  # Other coordinates are skipped, even before the sensitivity
  r$ci <- ci.coords(r.s100b, x = c(0.5, 0.9), input = "specificity", ret = c("precision", "sensitivity"), boot.n = 10, progress = "none")
  expect_warning(
    expect_doppelganger("plot-roc-ci-coords-sensitivity-not-precision", function() plot(r)),
    "not 'precision' with input 'specificity'"
  )
})

test_that("plot.roc skips with a warning a ci.coords that does not map onto the ROC axes", {
  set.seed(42)
  pdf(NULL)
  on.exit(dev.off())
  r <- r.s100b
  # CI of the specificity at given sensitivities
  r$ci <- ci.coords(r.s100b, x = c(0.5, 0.9), input = "sensitivity", ret = "specificity", boot.n = 10, progress = "none")
  expect_warning(expect_identical(plot(r), r), "not 'specificity' with input 'sensitivity'")
  expect_warning(plot(r, ci.type = "shape"), "not 'specificity' with input 'sensitivity'")
  # CI at a threshold
  r$ci <- ci.coords(r.s100b, x = 0.2, boot.n = 10, progress = "none")
  expect_warning(plot(r), "not 'threshold', 'specificity', 'sensitivity' with input 'threshold'")
  # CI at the best threshold (perfect separation: a single best threshold in every resample)
  r <- roc(rep(0:1, each = 10), c(1:10, 11:20), quiet = TRUE)
  r$ci <- expect_warning(ci.coords(r, x = "best", boot.n = 10, progress = "none"), "AUC == 1")
  expect_warning(expect_identical(plot(r), r), "with input 'threshold'")
  # Nothing is drawn
  skip_if(getRversion() < "4.1")
  r <- r.s100b
  r$ci <- ci.coords(r.s100b, x = 0.2, boot.n = 10, progress = "none")
  expect_warning(
    expect_doppelganger("plot-roc-ci-coords-threshold-skipped", function() plot(r)),
    "with input 'threshold'"
  )
})

test_that("panel.first and panel.last are drawn on smoothed curves from roc() and plot.roc()", {
  pdf(NULL)
  on.exit(dev.off())
  n <- 0
  s <- roc(aSAH$outcome, aSAH$s100b,
    smooth = TRUE, plot = TRUE, quiet = TRUE,
    panel.first = {
      n <- n + 1
      abline(h = 0.5)
    },
    panel.last = text(0.5, 0.5, "last")
  )
  expect_equal(n, 1)
  expect_s3_class(s, "smooth.roc")
  expect_null(s$smoothing.args$panel.first)
  expect_null(s$smoothing.args$panel.last)
  skip_if(getRversion() < "4.1")
  expect_doppelganger("roc-smooth-plot-panel-first-panel-last", function() {
    roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, plot = TRUE, quiet = TRUE, panel.first = abline(h = 0.5), panel.last = text(0.5, 0.5, "last"))
  })
  expect_doppelganger("plot-roc-default-smooth-panel-last", function() {
    plot.roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, quiet = TRUE, panel.last = text(0.5, 0.5, "last"))
  })
})

test_that("plot.roc and lines.roc default methods keep the variable names of smoothed curves", {
  pdf(NULL)
  on.exit(dev.off())
  dataline <- function(x) grep("^Data:", capture.output(print(x)), value = TRUE)
  expected <- "Data: aSAH$s100b in 72 controls (aSAH$outcome Good) < 41 cases (aSAH$outcome Poor)."
  expect_identical(dataline(roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, quiet = TRUE)), expected)
  expect_identical(dataline(plot.roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, quiet = TRUE)), expected)
  expect_identical(dataline(lines.roc(aSAH$outcome, aSAH$s100b, smooth = TRUE, quiet = TRUE)), expected)
})
