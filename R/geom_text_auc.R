geom_text_auc <- function(data, ...) {
  UseMethod("geom_text_auc")
}

geom_text_auc.auc <- function(data, ci = NULL, x = NULL, y = NULL, pattern = NULL, hjust = 0, vjust = 1, ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  percent <- attr(data, "percent")
  one <- if (percent) 100 else 1
  if (is.null(x)) {
    x <- one / 2
  }
  if (is.null(y)) {
    y <- one / 2
  }
  label <- ggroc_auc_label(data, ci, pattern)
  new_ggroc_layer(
    behind = FALSE,
    make_layers = function(plot) {
      list(do.call(
        ggplot2::annotate,
        c(
          list("text", x = x, y = y, label = label, hjust = hjust, vjust = vjust),
          extras
        )
      ))
    }
  )
}

geom_text_auc.roc <- function(data, ...) {
  geom_text_auc(data$auc, ci = data$ci, ...)
}

geom_text_auc.smooth.roc <- geom_text_auc.roc
