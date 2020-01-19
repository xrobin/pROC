# Skip slow tests
skip_slow <- function() {
	skip_if_not(exists("run_slow_tests") && run_slow_tests, message = "Slow test skipped")
}
