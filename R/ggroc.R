
#' Returns the coords as a data.frame in the right ordering for ggplot2 
get.coords.for.ggplot <- function(roc) {
	df <- coords(roc, "all", transpose = FALSE)
	return(df[rev(seq(nrow(df))),])
}

get.aes.for.ggplot <- function(roc, legacy.axes, extra_aes = c()) {
	# Prepare the aesthetics
	if(roc$percent) {
		if (legacy.axes) {
			aes_list <- list(x = "1-specificity", y = "sensitivity")
			xlims <- ggplot2::scale_x_continuous(lim=c(0, 100))		
		}
		else {
			aes_list <- list(x = "specificity", y = "sensitivity")
			xlims <- ggplot2::scale_x_reverse(lim=c(100, 0))
		}
	}
	else {
		if (legacy.axes) {
			aes_list <- list(x = "1-specificity", y = "sensitivity")
			xlims <- ggplot2::scale_x_continuous(lim=c(0, 1))
		}
		else {
			aes_list <- list(x = "specificity", y = "sensitivity")
			xlims <- ggplot2::scale_x_reverse(lim=c(1, 0))
		}
	}
	# Add extra aes
	for (ae in extra_aes) {
		aes_list[[ae]] <- "name"
	}
	aes <- do.call(ggplot2::aes_string, aes_list)
	
	return(list(aes=aes, xlims=xlims))
}

load.ggplot2 <- function() {
	if (! isNamespaceLoaded("ggplot2")) {
		message('You may need to call library(ggplot2) if you want to add layers, etc.')
	}
	load.suggested.package("ggplot2")
}

ggroc <- function(data, ...) {
	UseMethod("ggroc")
}

ggroc.roc <- function(data, legacy.axes = FALSE, ...) {
	load.ggplot2()
	# Get the roc data with coords
	df <- get.coords.for.ggplot(data)

	# Prepare the aesthetics
	aes <- get.aes.for.ggplot(data, legacy.axes)

	# Do the plotting
	ggplot2::ggplot(df) + ggplot2::geom_line(aes$aes, ...) + aes$xlims
}

ggroc.smooth.roc <- ggroc.roc

ggroc.list <- function(data, aes = c("colour", "alpha", "linetype", "size", "group"), legacy.axes = FALSE, ...) {
	load.ggplot2()
	if (missing(aes)) {
		aes <- "colour"
	}
	aes <- sub("color", "colour", aes)
	aes <- match.arg(aes, several.ok = TRUE)

	# Make sure data is a list and every element is a roc object
	if (! all(sapply(data, methods::is, "roc") | sapply(data, methods::is, "smooth.roc"))) {
		stop("All elements in 'data' must be 'roc' objects.")
	}
	
	# Make sure percent is consistent
	percents <- sapply(data, `[[`, "percent")
	if (!(all(percents) || all(!percents))) {
		stop("ROC curves use percent inconsistently and cannot be plotted together")
	}
	
	# Make sure the data is a named list
	if (is.null(names(data))) {
		names(data) <- seq(data)
	}
	# Make sure names are unique:
	if (any(duplicated(names(data)))) {
		stop("Names of 'data' must be unique")
	}
	
	# Get the coords
	coord.dfs <- sapply(data, get.coords.for.ggplot, simplify = FALSE)

	# Add a "name" colummn
	for (i in seq_along(coord.dfs)) {
		coord.dfs[[i]]$name <- names(coord.dfs)[i]
	}
	
	# Make a big data.frame
	coord.dfs <- do.call(rbind, coord.dfs)
    coord.dfs$name <- factor(coord.dfs$name, as.vector(names(data)))
	
	# Prepare the aesthetics
	aes.ggplot <- get.aes.for.ggplot(data[[1]], legacy.axes, aes)

	# Do the plotting
	ggplot2::ggplot(coord.dfs, aes.ggplot$aes) + ggplot2::geom_line(...) + aes.ggplot$xlims
}
