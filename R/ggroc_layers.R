# Shared ggplot2 helpers for ggroc layers.
# Layers emit columns: specificity, sensitivity, 1-specificity.
# Background layers are inserted immediately before the ROC GeomLine.

ggroc_layer_geom_class <- function(layer) {
  class(layer$geom)[[1]]
}

ggroc_roc_line_index <- function(layers) {
  classes <- vapply(layers, ggroc_layer_geom_class, character(1))
  i <- match("GeomLine", classes)
  if (is.na(i)) {
    length(layers) + 1L
  } else {
    i
  }
}

# Move the last layer to just before the first GeomLine in the remaining layers.
ggroc_insert_last_before_roc_line <- function(layers) {
  n <- length(layers)
  if (n <= 1L) {
    return(layers)
  }
  new_layer <- layers[[n]]
  rest <- layers[-n]
  roc_i <- ggroc_roc_line_index(rest)
  append(rest, list(new_layer), after = roc_i - 1L)
}

ggroc_plot_mapping <- function(plot) {
  mapping <- plot@mapping
  if (is.null(mapping$x)) {
    for (layer in plot@layers) {
      if (!is.null(layer$mapping$x)) {
        return(layer$mapping)
      }
    }
  }
  mapping
}

ggroc_legacy_axes_from_plot <- function(plot) {
  mapping <- ggroc_plot_mapping(plot)
  if (is.null(mapping$x)) {
    return(FALSE)
  }
  grepl("1-specificity", rlang::as_label(mapping$x), fixed = TRUE)
}

ggroc_x_col <- function(plot) {
  if (ggroc_legacy_axes_from_plot(plot)) {
    "1-specificity"
  } else {
    "specificity"
  }
}

ggroc_xy_aes <- function(plot, extra = list()) {
  xcol <- ggroc_x_col(plot)
  .data <- rlang::.data
  aes_list <- c(
    list(
      x = ggplot2::expr(.data[[xcol]]),
      y = ggplot2::expr(.data[["sensitivity"]])
    ),
    extra
  )
  do.call(ggplot2::aes, aes_list)
}

ggroc_add_one_or_hundred <- function(df, percent) {
  one <- if (percent) 100 else 1
  df[["1-specificity"]] <- one - df[["specificity"]]
  df
}

ggroc_errorbar_width <- function(percent) {
  0.02 * if (percent) 100 else 1
}

new_ggroc_layer <- function(make_layers, behind = FALSE) {
  structure(
    list(make_layers = make_layers, behind = behind),
    class = "ggroc_layer"
  )
}

ggplot_add.ggroc_layer <- function(object, plot, ...) {
  layers <- object$make_layers(plot)
  if (inherits(layers, "Layer") || inherits(layers, "LayerInstance")) {
    layers <- list(layers)
  }
  behind <- isTRUE(object$behind)
  for (layer in layers) {
    plot <- ggplot2::ggplot_add(layer, plot, ...)
    if (behind) {
      plot@layers <- ggroc_insert_last_before_roc_line(plot@layers)
    }
  }
  plot
}

ggroc_axis_labs <- function(roc, legacy.axes) {
  percent <- roc$percent
  xlab <- if (percent) {
    if (legacy.axes) "100 - Specificity (%)" else "Specificity (%)"
  } else {
    if (legacy.axes) "1 - Specificity" else "Specificity"
  }
  ylab <- if (percent) "Sensitivity (%)" else "Sensitivity"
  ggplot2::labs(x = xlab, y = ylab)
}

ggroc_auc_label <- function(auc, ci = NULL, pattern = NULL) {
  percent <- attr(auc, "percent")
  partial.auc <- attr(auc, "partial.auc")
  has.ci <- !is.null(ci) && methods::is(ci, "ci.auc")
  if (is.null(pattern)) {
    pattern <- ifelse(identical(partial.auc, FALSE), "AUC: ", "Partial AUC: ")
    pattern <- paste0(pattern, ifelse(percent, "%.1f%%", "%.3f"))
    if (has.ci) {
      pattern <- paste0(
        pattern, " (",
        ifelse(percent, "%.1f%%", "%.3f"), "\u2013",
        ifelse(percent, "%.1f%%", "%.3f"), ")"
      )
    }
  }
  if (has.ci) {
    sprintf(pattern, auc, ci[1], ci[3])
  } else {
    sprintf(pattern, auc)
  }
}
