remove.calls.recursive <- function(x) {
	if (is.null(x)) return(NULL)
	attr(x, "roc") <- remove.calls.recursive(attr(x, "roc"))
	attr(x, "auc") <- remove.calls.recursive(attr(x, "auc"))
	attr(x, "ci") <- remove.calls.recursive(attr(x, "ci"))
	if (!is.list(x)) return(x)
	x$roc <- remove.calls.recursive(x$roc)
	x$auc <- remove.calls.recursive(x$auc)
	x$ci <- remove.calls.recursive(x$ci)
	x$call <- NULL
	return(x)
}

expect_equal_ignore_call <- function(x, y, ...) {
	x <- remove.calls.recursive(x)
	y <- remove.calls.recursive(y)
	expect_equal(x, y, ...)
}
