ggplot.roc <- function(data, ...) {
	# Get the roc data with coords
	df <- as.data.frame(t(coords(data, "all")), row.names = NA)

	# Prepare the aesthetics
	if(data$percent) {
		aes <- ggplot2::aes_string(x = "100 - specificity", y = "sensitivity")
	}
	else {
		aes <- ggplot2::aes_string(x = "1 - specificity", y = "sensitivity")
	}

	# Do the plotting
	ggplot(df[rev(seq(nrow(df))),]) + ggplot2::geom_line(aes, ...)
	# Or with ggvis:
	# ggvis(df[rev(seq(nrow(df))),], ~1-specificity, ~sensitivity) %>% layer_lines()
}