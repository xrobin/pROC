# Wrapper for test_that that disables/hides the progress bar.
# This makes the output of devtools::test() more readable
# Cleans up: restores the previous progress at the end.
test_that_no_progress <- function(test, code) {
	old_progress <- getOption("pROCProgress")
	options(pROCProgress = list(name = "none"))
	tryCatch({
			test_that(test, code)
		},
		finally = function() {
			options(pROCProgress = old_progress)
		}
	)
}
