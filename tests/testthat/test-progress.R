library(pROC)
data(aSAH)

# Progress bars were lost in 1.19 when plyr was dropped: '.progress' was an
# argument to plyr's laply(), so it went with it. Now that every bootstrap
# iterates through one dispatcher, the bar is back without the dependency.

B <- 20

# Capture what a call prints, so the bar can be tested without a terminal.
progress_output <- function(expr) {
  # invisible() so that the returned ci object is not printed into the capture:
  # only what the bar itself writes should end up here.
  out <- utils::capture.output(invisible(force(expr)), type = "output")
  paste(out, collapse = "")
}

test_that("progress = TRUE draws a bar", {
  drawn <- progress_output(ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = TRUE))
  expect_match(drawn, "=")
  expect_match(drawn, "100%")
})


test_that("progress = FALSE draws nothing", {
  expect_identical(
    progress_output(ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = FALSE)),
    ""
  )
})


test_that("the bar does not change the result", {
  set.seed(42)
  with.bar <- suppressWarnings(utils::capture.output(
    quiet <- ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = TRUE)
  ))
  set.seed(42)
  without.bar <- ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = FALSE)
  expect_equal(as.numeric(quiet), as.numeric(without.bar))
})


test_that("every bootstrap entry point accepts progress", {
  # The point of the dispatcher: one place to add this, not two dozen.
  expect_match(progress_output(ci.se(r.wfns, boot.n = B, progress = TRUE)), "100%")
  expect_match(progress_output(ci.sp(r.wfns, boot.n = B, progress = TRUE)), "100%")
  expect_match(progress_output(ci.thresholds(r.wfns, boot.n = B, progress = TRUE)), "100%")
  expect_match(progress_output(ci.coords(r.wfns, x = 0.5, input = "specificity",
                                         boot.n = B, progress = TRUE)), "100%")
  expect_match(progress_output(var(r.wfns, method = "bootstrap", boot.n = B,
                                   progress = TRUE)), "100%")
  expect_match(progress_output(cov(r.wfns, r.ndka, method = "bootstrap", boot.n = B,
                                   progress = TRUE)), "100%")
  expect_match(progress_output(suppressWarnings(
    roc.test(r.wfns, r.ndka, method = "bootstrap", boot.n = B, progress = TRUE))), "100%")
})


test_that("smoothed bootstraps accept progress too", {
  s <- smooth(r.ndka)
  expect_match(progress_output(ci.auc(s, boot.n = B, progress = TRUE)), "100%")
  expect_match(progress_output(ci.se(s, boot.n = B, progress = TRUE)), "100%")
})


test_that("'progress' is no longer deprecated", {
  expect_no_warning(ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = FALSE))
  expect_no_warning(utils::capture.output(
    ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = TRUE)
  ))
})


test_that("the pre-1.19 plyr spellings are still understood", {
  # Old scripts and .Rprofile settings should keep working rather than error.
  expect_false(roc_utils_normalise_progress("none"))
  expect_true(roc_utils_normalise_progress("text"))
  expect_true(roc_utils_normalise_progress("win"))
  expect_true(roc_utils_normalise_progress("tk"))
  # The option used to hold a list, not a name.
  expect_true(roc_utils_normalise_progress(list(name = "text", width = NA)))
  expect_false(roc_utils_normalise_progress(list(name = "none")))
  # And the supported forms.
  expect_true(roc_utils_normalise_progress(TRUE))
  expect_false(roc_utils_normalise_progress(FALSE))
  expect_false(roc_utils_normalise_progress(NULL))
})


test_that("a character progress value still draws a bar", {
  expect_match(
    progress_output(ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = "text")),
    "100%"
  )
  expect_identical(
    progress_output(ci.auc(r.wfns, method = "bootstrap", boot.n = B, progress = "none")),
    ""
  )
})


test_that("the pROCProgress option sets the default", {
  withr_option <- getOption("pROCProgress")
  on.exit(options(pROCProgress = withr_option))

  options(pROCProgress = TRUE)
  expect_match(progress_output(ci.auc(r.wfns, method = "bootstrap", boot.n = B)), "100%")

  options(pROCProgress = FALSE)
  expect_identical(progress_output(ci.auc(r.wfns, method = "bootstrap", boot.n = B)), "")

  # Unset, the default follows interactive(), which is FALSE under testthat.
  options(pROCProgress = NULL)
  expect_identical(progress_output(ci.auc(r.wfns, method = "bootstrap", boot.n = B)), "")
})


test_that("attaching pROC no longer removes the pROCProgress option", {
  # 1.19's .onAttach deleted it and warned; it is a supported option again.
  old <- getOption("pROCProgress")
  on.exit(options(pROCProgress = old))
  options(pROCProgress = TRUE)
  suppressMessages(pROC:::.onAttach(NULL, "pROC"))
  expect_true(getOption("pROCProgress"))
})
