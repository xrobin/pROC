geom_roc_threshold <- function(data, ...) {
  UseMethod("geom_roc_threshold")
}

geom_roc_threshold.roc <- function(data,
                                   thresholds = "best",
                                   pattern = NULL,
                                   hjust = -0.05,
                                   vjust = 1.25,
                                   size = 2,
                                   ...) {
  load.ggplot2()
  # A smoothed curve has no thresholds: only "best" can be marked (as with plot.roc's print.thres)
  if (methods::is(data, "smooth.roc") && !identical(coords_special_x(thresholds, roc = data), "best")) {
    stop("Only 'best' 'thresholds' can be marked on a smoothed ROC curve: it has no thresholds.")
  }
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  percent <- data$percent
  co <- coords(data, thresholds)
  # Only mark real thresholds: with a partial AUC, coords(data, "all") also
  # returns the interpolated limits of the window, with an NA threshold
  if (!methods::is(data, "smooth.roc") && !is.null(co)) {
    co <- co[!is.na(co$threshold), , drop = FALSE]
  }
  co <- ggroc_add_one_or_hundred(co, percent)
  if (is.null(co$threshold)) {
    co$threshold <- NA_real_
  }
  co$label <- roc_utils_threshold_labels(pattern, co$threshold, co$specificity, co$sensitivity, percent)
  new_ggroc_layer(
    behind = FALSE,
    make_layers = function(plot) {
      xcol <- ggroc_x_col(plot)
      .data <- rlang::.data
      point_aes <- ggplot2::aes(
        x = .data[[xcol]],
        y = .data[["sensitivity"]]
      )
      text_aes <- ggplot2::aes(
        x = .data[[xcol]],
        y = .data[["sensitivity"]],
        label = .data[["label"]]
      )
      list(
        do.call(
          ggplot2::geom_point,
          c(
            list(data = co, mapping = point_aes, size = size, inherit.aes = FALSE),
            extras
          )
        ),
        do.call(
          ggplot2::geom_text,
          c(
            list(
              data = co, mapping = text_aes,
              hjust = hjust, vjust = vjust, inherit.aes = FALSE
            ),
            extras
          )
        )
      )
    }
  )
}

geom_roc_threshold.smooth.roc <- geom_roc_threshold.roc
