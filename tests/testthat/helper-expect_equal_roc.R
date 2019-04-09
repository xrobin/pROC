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

remove.response.names.recursive <- function(x) {
	if (is.null(x)) return(NULL)
	attr(x, "roc") <- remove.response.names.recursive(attr(x, "roc"))
	attr(x, "auc") <- remove.response.names.recursive(attr(x, "auc"))
	attr(x, "ci") <- remove.response.names.recursive(attr(x, "ci"))
	if (!is.list(x)) return(x)
	x$roc <- remove.response.names.recursive(x$roc)
	x$auc <- remove.response.names.recursive(x$auc)
	x$ci <- remove.response.names.recursive(x$ci)
	names(x$response) <- NULL
	names(x$original.response) <- NULL
	return(x)
}


expect_equal_ignore_call <- function(x, y, ...) {
	x <- remove.calls.recursive(x)
	y <- remove.calls.recursive(y)
	expect_equal(x, y, ...)
}

expect_equal_roc_formula <- function(x, y, ...) {
	# roc.formula adds names to response and original.response
	# this expectation ignores them, as well as the call
	x <- remove.calls.recursive(x)
	x <- remove.response.names.recursive(x)
	y <- remove.calls.recursive(y)
	y <- remove.response.names.recursive(y)
	expect_equal(x, y, ...)
}
