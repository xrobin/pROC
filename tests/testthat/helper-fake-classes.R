# Emulate a vctrs-classed numeric response such as haven::labelled (as read
# from SPSS/Stata/SAS files) without depending on haven: subsetting keeps
# the class, and == against a character vector errors like vctrs does.
# Returns the constructor; the methods are removed at the end of 'frame'.
local_fake_labelled <- function(frame = parent.frame()) {
  skip_if_not_installed("testthat", "3.2.0") # local_mocked_s3_method
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
