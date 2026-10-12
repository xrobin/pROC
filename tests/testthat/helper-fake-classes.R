# Emulate a vctrs-classed numeric response such as haven::labelled (as read
# from SPSS/Stata/SAS files) without depending on haven: subsetting keeps
# the class, and == against a character vector errors like vctrs does.
# Returns the constructor; the methods are removed at the end of 'frame'.
local_fake_labelled <- function(frame = parent.frame()) {
  skip_if_no_mocked_s3_method()
  local_mocked_s3_method("Ops", "pROC_fake_labelled", function(e1, e2) {
    if (is.character(e1) || (!missing(e2) && is.character(e2))) {
      stop("Can't combine `x` <character> and `y` <double>.")
    }
    if (missing(e2)) {
      return(get(.Generic)(unclass(e1)))
    }
    get(.Generic)(unclass(e1), unclass(e2))
  }, frame = frame)
  local_mocked_s3_method("[", "pROC_fake_labelled", function(x, ...) {
    structure(unclass(x)[...], class = class(x))
  }, frame = frame)
  function(x) structure(x, class = "pROC_fake_labelled")
}

# local_mocked_s3_method() is not in every testthat 3.2.x (missing in 3.2.1):
# test for the function itself rather than for a version
skip_if_no_mocked_s3_method <- function() {
  if (!exists("local_mocked_s3_method", envir = asNamespace("testthat"), inherits = FALSE)) {
    skip("testthat::local_mocked_s3_method() is not available")
  }
}

# local_mocked_bindings() needs testthat >= 3.1.7
skip_if_no_mocked_bindings <- function() {
  if (!exists("local_mocked_bindings", envir = asNamespace("testthat"), inherits = FALSE)) {
    skip("testthat::local_mocked_bindings() is not available")
  }
}
