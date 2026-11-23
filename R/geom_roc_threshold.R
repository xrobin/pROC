geom_roc_threshold <- function(data, ...) {
  UseMethod("geom_roc_threshold")
}

geom_roc_threshold.roc <- function(data,
                                   thresholds = "best",
                                   pattern = NULL,
                                   hjust = -0.05,
                                   vjust = 1.25,
                                   size = 2,
                                   best.method = c("youden", "closest.topleft"),
                                   best.weights = c(1, 0.5),
                                   mapping = NULL,
                                   ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  percent <- data$percent
  co <- roc_utils_marked_thresholds(data, thresholds, best.method, best.weights)
  co <- ggroc_add_one_or_hundred(co, percent)
  co$label <- roc_utils_threshold_labels(pattern, co$threshold, co$specificity, co$sensitivity, percent)
  new_ggroc_layer(
    behind = FALSE,
    make_layers = function(plot, ...) {
      xcol <- ggroc_x_col(plot)
      .data <- rlang::.data
      point_aes <- ggroc_add_mapping(ggplot2::aes(
        x = .data[[xcol]],
        y = .data[["sensitivity"]]
      ), mapping)
      text_aes <- ggroc_add_mapping(ggplot2::aes(
        x = .data[[xcol]],
        y = .data[["sensitivity"]],
        label = .data[["label"]]
      ), mapping)
      list(
        do.call(
          ggplot2::geom_point,
          c(
            list(data = co, mapping = point_aes, inherit.aes = FALSE),
            if (is.null(mapping$size)) list(size = size),
            extras
          )
        ),
        do.call(
          ggplot2::geom_text,
          c(
            list(data = co, mapping = text_aes, inherit.aes = FALSE),
            if (is.null(mapping$hjust)) list(hjust = hjust),
            if (is.null(mapping$vjust)) list(vjust = vjust),
            extras
          )
        )
      )
    }
  )
}

geom_roc_threshold.smooth.roc <- geom_roc_threshold.roc

geom_roc_threshold.auc <- function(data, ...) {
  # The curve, with the (partial) AUC specification of the auc object
  roc <- attr(data, "roc")
  roc$auc <- data
  geom_roc_threshold(roc, ...)
}

geom_roc_threshold.list <- function(data, ...) {
  ggroc_layer_list(data, geom_roc_threshold, ...)
}
