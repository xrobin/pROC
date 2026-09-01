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
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  percent <- data$percent
  if (is.null(pattern)) {
    pattern <- if (percent) {
      "%.1f (%.1f%%, %.1f%%)"
    } else {
      "%.3f (%.3f, %.3f)"
    }
  }
  co <- coords(data, thresholds)
  co <- ggroc_add_one_or_hundred(co, percent)
  if (is.null(co$threshold)) {
    co$threshold <- NA_real_
  }
  co$label <- sprintf(pattern, co$threshold, co$specificity, co$sensitivity)
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
