context("plot.auc")

# polygon area via the shoelace formula, used to cross-check that the
# partial-AUC polygon vertices (the riskiest part of polygon_auc) enclose
# exactly the AUC's own area. Only valid for uncorrected partial AUCs:
# partial.auc.correct rescales the AUC value (McClish), so it is no longer
# the geometric area of the drawn shape.
polygon_area <- function(x, y) {
  abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
}

# trace()'s tracer always runs in the global environment, so a plain local
# `captured <<- ...` would silently create/update a *different* global
# variable instead of the caller's local one. Go through an explicit,
# dedicated global binding instead.
capture_polygon_call <- function(expr) {
  assign(".test_captured_polygon", NULL, envir = globalenv())
  trace(graphics::polygon,
    tracer = quote(assign(".test_captured_polygon", list(x = x, y = y), envir = globalenv())),
    print = FALSE
  )
  on.exit({
    suppressMessages(try(untrace(graphics::polygon), silent = TRUE))
    rm(".test_captured_polygon", envir = globalenv())
  })
  force(expr)
  get(".test_captured_polygon", envir = globalenv())
}

test_that("text(auc), text(ci.auc) and text(roc) match ggroc_auc_label()", {
  pdf(NULL)
  on.exit(dev.off())

  plot(r.s100b)
  expect_identical(text(r.s100b$auc), ggroc_auc_label(r.s100b$auc))
  expect_identical(text(r.s100b), ggroc_auc_label(r.s100b$auc, r.s100b$ci))

  plot(r.s100b.percent)
  expect_identical(text(r.s100b.percent$auc), ggroc_auc_label(r.s100b.percent$auc))

  plot(r.s100b.partial1)
  expect_identical(text(r.s100b.partial1$auc), ggroc_auc_label(r.s100b.partial1$auc))

  roc_ci <- roc(aSAH$outcome, aSAH$s100b, ci = TRUE, boot.n = 100, quiet = TRUE)
  plot(roc_ci)
  expect_identical(text(roc_ci), ggroc_auc_label(roc_ci$auc, roc_ci$ci))
  expect_identical(text(roc_ci$ci), ggroc_auc_label(roc_ci$auc, roc_ci$ci))
})

test_that("text() ignores a non-ci.auc confidence interval", {
  pdf(NULL)
  on.exit(dev.off())

  roc_se_ci <- roc(aSAH$outcome, aSAH$s100b, quiet = TRUE)
  roc_se_ci$ci <- ci.se(roc_se_ci, boot.n = 100)
  plot(roc_se_ci)
  expect_identical(text(roc_se_ci), ggroc_auc_label(roc_se_ci$auc))
})

test_that("polygon_auc partial polygon area matches the (uncorrected) partial AUC: focus sp", {
  pdf(NULL)
  on.exit(dev.off())

  plot(r.s100b.partial1)
  captured <- capture_polygon_call(polygon_auc(r.s100b.partial1))
  expect_equal(polygon_area(captured$x, captured$y), as.numeric(r.s100b.partial1$auc))
})

test_that("polygon_auc partial polygon area matches the (uncorrected) partial AUC: focus se, percent", {
  pdf(NULL)
  on.exit(dev.off())

  plot(r.s100b.percent.partial2)
  captured <- capture_polygon_call(polygon_auc(r.s100b.percent.partial2))
  # percent-scale coordinates: area is 100x the (already /100) auc value
  expect_equal(polygon_area(captured$x, captured$y), as.numeric(r.s100b.percent.partial2$auc) * 100)
})

test_that("polygon_auc/polygon_max_auc/text work on 'auc', 'roc' and 'smooth.roc'", {
  pdf(NULL)
  on.exit(dev.off())

  aucobj <- r.s100b$auc
  plot(r.s100b)
  expect_identical(polygon_auc(aucobj), aucobj)
  expect_identical(polygon_max_auc(aucobj), aucobj)
  expect_identical(polygon_auc(r.s100b), r.s100b)
  expect_identical(polygon_max_auc(r.s100b), r.s100b)

  sm <- smooth(r.s100b)
  plot(sm)
  expect_error(polygon_auc(sm), NA)
  expect_error(polygon_max_auc(sm), NA)
  expect_error(text(sm), NA)
})

test_that("panel.first is evaluated once, under the identity line and the curve", {
  pdf(NULL)
  on.exit(dev.off())

  count <- 0
  plot.roc(aSAH$outcome, aSAH$s100b, quiet = TRUE, panel.first = {
    count <- count + 1
  })
  expect_equal(count, 1)
})

test_that("panel.first is not forwarded through '...': no call when combined with add=TRUE", {
  pdf(NULL)
  on.exit(dev.off())

  plot.roc(aSAH$outcome, aSAH$s100b, quiet = TRUE)
  count <- 0
  plot.roc(aSAH$outcome, aSAH$ndka,
    quiet = TRUE, add = TRUE,
    panel.first = {
      count <- count + 1
    }
  )
  expect_equal(count, 0)
})

test_that("panel.first works through plot.roc(response, predictor, ...)", {
  pdf(NULL)
  on.exit(dev.off())

  count <- 0
  plot.roc(aSAH$outcome, aSAH$s100b, quiet = TRUE, panel.first = {
    count <- count + 1
  })
  expect_equal(count, 1)
})

test_that("panel.last is evaluated once, on top of everything plot.roc draws", {
  pdf(NULL)
  on.exit(dev.off())

  count <- 0
  plot.roc(aSAH$outcome, aSAH$s100b,
    quiet = TRUE, print.auc = TRUE, print.thres = TRUE,
    panel.last = {
      count <- count + 1
    }
  )
  expect_equal(count, 1)
})

test_that("panel.last is not forwarded through '...'", {
  pdf(NULL)
  on.exit(dev.off())

  plot.roc(aSAH$outcome, aSAH$s100b, quiet = TRUE)
  count <- 0
  plot.roc(aSAH$outcome, aSAH$ndka,
    quiet = TRUE, add = TRUE,
    panel.last = {
      count <- count + 1
    }
  )
  expect_equal(count, 1) # unlike panel.first, panel.last fires on every call, add=TRUE included
})

test_that("panel.last works through plot.roc(response, predictor, ...)", {
  pdf(NULL)
  on.exit(dev.off())

  count <- 0
  plot.roc(aSAH$outcome, aSAH$s100b, quiet = TRUE, panel.last = {
    count <- count + 1
  })
  expect_equal(count, 1)
})

test_that("polygon_auc/polygon_max_auc/text error on a 'multiclass.auc'", {
  pdf(NULL)
  on.exit(dev.off())

  m.auc <- suppressWarnings(multiclass.roc(aSAH$gos6, aSAH$s100b, quiet = TRUE))$auc
  plot(r.s100b)
  expect_error(polygon_auc(m.auc), "multiclass.auc")
  expect_error(polygon_max_auc(m.auc), "multiclass.auc")
  expect_error(text(m.auc), "multiclass.auc")
})

test_that("polygon_auc/polygon_max_auc/text error when no device is open", {
  skip_if(dev.cur() != 1, "a device is already open outside this test")
  expect_error(polygon_auc(r.s100b$auc), "plot")
  expect_error(polygon_max_auc(r.s100b$auc), "plot")
  expect_error(text(r.s100b$auc), "plot")
})

test_that("plot.roc(auc.polygon=TRUE, max.auc.polygon=TRUE, print.auc=TRUE) matches the add-on functions via panel.first", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")

  all_in_one <- function() {
    plot(r.s100b.percent.partial1,
      reuse.auc = TRUE, auc.polygon = TRUE, max.auc.polygon = TRUE, print.auc = TRUE
    )
  }
  composed <- function() {
    plot(r.s100b.percent.partial1,
      panel.first = {
        polygon_max_auc(r.s100b.percent.partial1)
        polygon_auc(r.s100b.percent.partial1)
      }
    )
    text(r.s100b.percent.partial1)
  }

  f1 <- tempfile(fileext = ".svg")
  f2 <- tempfile(fileext = ".svg")
  vdiffr:::write_svg(all_in_one, f1, title = "all_in_one")
  vdiffr:::write_svg(composed, f2, title = "composed")
  expect_identical(readLines(f1), readLines(f2))
})

test_that("the target documentation example renders without error", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")
  skip_slow()

  test_target_example <- function() {
    roc4 <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE, quiet = TRUE)
    auc4 <- auc(roc4, partial.auc = c(100, 90), partial.auc.correct = TRUE, partial.auc.focus = "sens")
    set.seed(42)
    ci4 <- ci(auc4, boot.n = 100, conf.level = 0.9, boot.stratified = FALSE)
    plot(roc4,
      grid = TRUE, print.thres = TRUE,
      panel.first = {
        polygon_max_auc(auc4)
        polygon_auc(auc4)
      }
    )
    text(ci4)
  }
  expect_doppelganger("plot.auc.target.example", test_target_example)
})

test_that("polygon_ci() never redraws the curve, unlike plot.ci(type='shape', no.roc=FALSE)", {
  pdf(NULL)
  on.exit(dev.off())

  ci_se <- ci.se(r.s100b, specificities = seq(0, 1, .05), boot.n = 20)

  assign(".test_call_count", 0L, envir = globalenv())
  trace(pROC:::plot.roc.roc,
    tracer = quote(assign(".test_call_count", get(".test_call_count", envir = globalenv()) + 1L, envir = globalenv())),
    print = FALSE
  )
  on.exit({
    suppressMessages(try(untrace(pROC:::plot.roc.roc), silent = TRUE))
    rm(".test_call_count", envir = globalenv())
  })

  plot(r.s100b) # counts too, since the trace is already active; reset after
  assign(".test_call_count", 0L, envir = globalenv())
  plot(ci_se, type = "shape")
  expect_equal(get(".test_call_count", envir = globalenv()), 1) # default (no.roc=FALSE) redraws the curve once

  plot(r.s100b)
  assign(".test_call_count", 0L, envir = globalenv())
  polygon_ci(ci_se)
  expect_equal(get(".test_call_count", envir = globalenv()), 0) # polygon_ci() never redraws
})

test_that("polygon_ci() preserves custom curve styling when used via panel.first", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")

  # Tracing graphics::lines.default does not actually intercept lines()
  # calls (confirmed even for a plain, non-pROC lines() call) -- some base
  # generics dispatch in a way R-level trace() can't see. Inspect the
  # rendered SVG instead: a single <polyline> in the curve's own color
  # means it was drawn exactly once, with no default-styled redraw.
  ci_se <- ci.se(r.s100b, specificities = seq(0, 1, .05), boot.n = 20)
  f <- tempfile(fileext = ".svg")
  vdiffr:::write_svg(function() {
    plot(r.s100b, col = "red", lwd = 4, panel.first = polygon_ci(ci_se))
  }, f, title = "t")
  polylines <- grep("<polyline", readLines(f), value = TRUE, fixed = TRUE)
  expect_length(polylines, 1)
  expect_match(polylines, "#FF0000", fixed = TRUE) # red, preserved
})

test_that("polygon_ci works on ci.se and ci.sp", {
  pdf(NULL)
  on.exit(dev.off())

  ci_se <- ci.se(r.s100b, specificities = seq(0, 1, .05), boot.n = 20)
  ci_sp <- ci.sp(r.s100b, sensitivities = seq(0, 1, .05), boot.n = 20)
  plot(r.s100b)
  expect_identical(polygon_ci(ci_se), ci_se)
  expect_identical(polygon_ci(ci_sp), ci_sp)
})

test_that("polygon_ci errors when no device is open", {
  skip_if(dev.cur() != 1, "a device is already open outside this test")
  ci_se <- ci.se(r.s100b, specificities = seq(0, 1, .05), boot.n = 20)
  expect_error(polygon_ci(ci_se), "plot")
})

test_that("plot.ci.se/plot.ci.sp/plot.ci.thresholds/plot.ci.coords error clearly when no device is open", {
  skip_if(dev.cur() != 1, "a device is already open outside this test")
  ci_se <- ci.se(r.s100b, specificities = seq(0, 1, .05), boot.n = 20)
  ci_sp <- ci.sp(r.s100b, sensitivities = seq(0, 1, .05), boot.n = 20)
  ci_th <- ci.thresholds(r.s100b, thresholds = "best", boot.n = 20)
  ci_co <- ci.coords(r.s100b, x = seq(0, 1, .2), input = "recall", ret = "precision", boot.n = 20)
  expect_error(plot(ci_se, type = "shape"), "plot")
  expect_error(plot(ci_sp, type = "bars"), "plot")
  expect_error(plot(ci_th), "plot")
  expect_error(plot(ci_co, type = "bars"), "plot")
})

test_that("polygon_ci via panel.first draws the curve once, unlike plot.ci(type='shape') drawn after", {
  skip_if_not_installed("vdiffr")
  skip_if(getRversion() < "4.1")

  ci_se <- ci.se(r.s100b, specificities = seq(0, 1, .1), boot.n = 20)

  after <- function() {
    plot(r.s100b)
    suppressWarnings(plot(ci_se, type = "shape", col = "gainsboro"))
  }
  panel_first <- function() {
    plot(r.s100b, panel.first = suppressWarnings(polygon_ci(ci_se, col = "gainsboro")))
  }

  f1 <- tempfile(fileext = ".svg")
  f2 <- tempfile(fileext = ".svg")
  vdiffr:::write_svg(after, f1, title = "after")
  vdiffr:::write_svg(panel_first, f2, title = "panel_first")
  count_polylines <- function(f) sum(grepl("<polyline", readLines(f), fixed = TRUE))
  # "after" redraws the curve on top of the shape (one extra <polyline>);
  # panel_first draws it exactly once. Everything else is identical.
  expect_equal(count_polylines(f1), count_polylines(f2) + 1)
})
