# Skip slow tests
skip_slow <- function(message = "Slow test skipped") {
	if (exists("run_slow_tests", envir=.GlobalEnv)) {
		if (! get("run_slow_tests", envir=.GlobalEnv)) {
			skip(message)
		}
	}
	else if (! identical(Sys.getenv("RUN_SLOW_TESTS"), "true")) {
		skip(message)
	}
}
