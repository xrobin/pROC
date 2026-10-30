library(pROC)
data(aSAH)

context("roc.utils")

test_that("roc_utils_thr_idx finds correc thresholds with direction=<", {
  obtained <- pROC:::roc_utils_thr_idx(r.s100b, c(-Inf, 0.205, 0.055, Inf))
  expect_equal(obtained, c(1, 18, 4, 51))
})

test_that("roc_utils_thr_idx finds correc thresholds with direction=>", {
  obtained <- pROC:::roc_utils_thr_idx(r.s100b, c(Inf, -Inf, 0.05, 0.055, 0.52, 0.205))
  expect_equal(obtained, c(51, 1, 3, 4, 40, 18))
})

test_that("roc_utils_calc_coords works", {
  obtained <- pROC:::roc_utils_calc_coords(r.s100b, -1:-4, c(1, .5, .1, 0), c(0, .5, .9, 1), c(12, .9))
  expect_equal(as.matrix(obtained), expected_roc_utils_calc_coords)
})

test_that("roc_utils_calc_coords works with percent", {
  obtained <- pROC:::roc_utils_calc_coords(r.s100b.percent, -1:-4, c(100, 50, 10, 0), c(0, 50, 90, 100), c(12, .9))
  expect_equal(as.matrix(obtained), expected_roc_utils_calc_coords.percent)
})

test_that("roc_utils_match_coords_input_args works", {
  expect_equal(pROC:::roc_utils_match_coords_input_args("t"), "threshold")
  expect_equal(pROC:::roc_utils_match_coords_input_args("threshold"), "threshold")
  expect_equal(pROC:::roc_utils_match_coords_input_args("fp"), "fp")
  expect_equal(pROC:::roc_utils_match_coords_input_args("1-se"), "1-sensitivity")
  for (coord in names(which(pROC:::coord.is.monotone))) {
    expect_equal(pROC:::roc_utils_match_coords_input_args(coord), coord)
  }

  # Errors
  # t with threshold=False
  expect_error(pROC:::roc_utils_match_coords_input_args("t", threshold = FALSE))
  # all only for ret
  expect_error(pROC:::roc_utils_match_coords_input_args("all"))
  # Only one allowed
  expect_error(pROC:::roc_utils_match_coords_input_args(c("specificity", "sensitivity")))
  # Invalid arg
  expect_error(pROC:::roc_utils_match_coords_input_args("blah"))
  # Not monotone
  expect_error(pROC:::roc_utils_match_coords_input_args("npe"))
  expect_error(pROC:::roc_utils_match_coords_input_args("accuracy"))
})


test_that("roc_utils_match_coords_ret_args works", {
  expect_equal(pROC:::roc_utils_match_coords_ret_args("t"), "threshold")
  expect_equal(pROC:::roc_utils_match_coords_ret_args("threshold"), "threshold")
  expect_equal(pROC:::roc_utils_match_coords_ret_args("fp"), "fp")
  expect_equal(pROC:::roc_utils_match_coords_ret_args("1-se"), "1-sensitivity")
  expect_equal(pROC:::roc_utils_match_coords_ret_args("npe"), "1-npv")
  for (coord in pROC:::roc.utils.valid.coords) {
    expect_equal(pROC:::roc_utils_match_coords_ret_args(coord), coord)
  }
  expect_equal(pROC:::roc_utils_match_coords_ret_args(pROC:::roc.utils.valid.coords), pROC:::roc.utils.valid.coords)

  # Errors
  # t with threshold=False
  expect_error(pROC:::roc_utils_match_coords_ret_args("t", threshold = FALSE))
  # Invalid arg
  expect_error(pROC:::roc_utils_match_coords_ret_args("blah"))
  # The following should be invalid but somehow it seems valid to say:
  # match.arg(c("sensitivity", "blah"), "sensitivity", TRUE)
  # and the extra 'blah' arg is ignored by match.arg.
  # Ignoring for now
  # expect_error(pROC:::roc_utils_match_coords_ret_args(c("sensitivity", "blah")))
})

test_that("roc_utils_bug saves the data in the temporary directory and stops", {
  # The non-interactive path; an interactive session asks before saving
  skip_if(interactive())
  pattern <- "^pROC_bug.*\\.RData$"
  before <- list.files(tempdir(), pattern = pattern, full.names = TRUE)
  wd.before <- list.files(getwd(), pattern = pattern)
  on.exit(unlink(setdiff(list.files(tempdir(), pattern = pattern, full.names = TRUE), before)))
  expect_error(
    pROC:::roc_utils_bug("internal check failed", x = 1:3, direction = "<"),
    "pROC: internal check failed. Data to reproduce it saved in .*pROC_bug.*\\.RData. Please report this bug to <https://github.com/xrobin/pROC/issues>"
  )
  new <- setdiff(list.files(tempdir(), pattern = pattern, full.names = TRUE), before)
  expect_length(new, 1)
  saved <- new.env()
  load(new, envir = saved)
  expect_setequal(ls(saved), c("x", "direction", "sessionInfo"))
  expect_identical(saved$x, 1:3)
  expect_identical(saved$direction, "<")
  expect_s3_class(saved$sessionInfo, "sessionInfo")
  expect_identical(list.files(getwd(), pattern = pattern), wd.before)
})
