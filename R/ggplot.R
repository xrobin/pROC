ggplot.roc <- function(data, ...) {
	# Get the roc data with coords
	df <- as.data.frame(t(coords(data, "all")), row.names = NA)

	# Prepare the aesthetics
	if(data$percent) {
		aes <- ggplot2::aes_string(x = "specificity", y = "sensitivity")
		xlims <- scale_x_reverse(lim=c(100, 0))
	}
	else {
		aes <- ggplot2::aes_string(x = "specificity", y = "sensitivity")
		xlims <- scale_x_reverse(lim=c(1, 0))
	}

	# Do the plotting
	ggplot(df[rev(seq(nrow(df))),]) + ggplot2::geom_line(aes, ...) + xlims
		
	# Or with ggvis:
	# ggvis(df[rev(seq(nrow(df))),], ~1-specificity, ~sensitivity) %>% layer_lines()
}