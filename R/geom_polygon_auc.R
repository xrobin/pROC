geom_polygon_auc <- function(data, ...) {
	UseMethod("geom_polygon_auc")
}

geom_polygon_auc.auc <- function(data, legacy.axes = FALSE, ...) {
	# Get the roc data with coords
	roc <- attr(data, "roc")
	roc$auc <- data
	df <- get.coords.for.ggplot(roc, ignore.partial.auc = FALSE)
	
	# Add bottom-right point
	partial.auc <- attr(data, "partial.auc")
	one.or.hundred <- ifelse(attr(data, "percent"), 100, 1)
	if (legacy.axes) {
		if (identical(partial.auc, FALSE)) {
			df[nrow(df) + 1, ] <- c(NA, one.or.hundred, 0, one.or.hundred)
		}
		else if (attr(data, "partial.auc.focus") == "sensitivity") {
			df[nrow(df) + c(1, 2), ] <- c(NA, NA, one.or.hundred, one.or.hundred, partial.auc, one.or.hundred, one.or.hundred)
		}
		else { # partial.auc.focus == "specificity"
			df[nrow(df) + c(1, 2), ] <- c(NA, NA, rev(partial.auc), 0, 0, one.or.hundred - rev(partial.auc))
		}
	}
	else {
		if (identical(partial.auc, FALSE)) {
			df[nrow(df) + 1, ] <- c(NA, 0, 0, 0)
		}
		else if (attr(data, "partial.auc.focus") == "sensitivity") {
			df[nrow(df) + c(1, 2), ] <- c(NA, NA, 0, 0, partial.auc, 0, 0)
		}
		else { # partial.auc.focus == "specificity"
			df[nrow(df) + c(1, 2), ] <- c(NA, NA, rev(partial.auc), 0, 0, one.or.hundred - rev(partial.auc))
		}
	}
	
	# Prepare the aesthetics
	aes <- get.aes.for.ggplot(attr(data, "roc"), legacy.axes)
	
	# Do the plotting
	ggplot2::geom_polygon(aes$aes, data=df, ...)
}

geom_polygon_auc.roc <- function(data, ...) {
	geom_polygon_auc(data$auc, ...)
}

geom_polygon_auc.smooth.roc <- geom_polygon_auc.roc
