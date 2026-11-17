geom_polygon_max_auc <- function(data, ...) {
  UseMethod("geom_polygon_max_auc")
}

geom_polygon_max_auc.auc <- function(data, fill = NULL, colour = NA, mapping = NULL, ...) {
  load.ggplot2()
  extras <- size.to.linewidth(list(...))
  names(extras) <- sub("color", "colour", names(extras))
  if (is.null(extras$fill)) {
    extras$fill <- fill
  }
  if (is.null(extras$colour) && is.null(mapping$colour)) {
    extras$colour <- colour
  }
  percent <- attr(data, "percent")
  partial.auc <- attr(data, "partial.auc")
  partial.auc.focus <- attr(data, "partial.auc.focus")
  new_ggroc_layer(
    behind = TRUE,
    make_layers = function(plot, ...) {
      one <- if (percent) 100 else 1
      if (identical(partial.auc, FALSE)) {
        df <- data.frame(
          specificity = c(one, one, 0, 0),
          sensitivity = c(0, one, one, 0)
        )
      } else if (identical(partial.auc.focus, "sensitivity")) {
        df <- data.frame(
          specificity = c(0, one, one, 0),
          sensitivity = c(partial.auc[2], partial.auc[2], partial.auc[1], partial.auc[1])
        )
      } else {
        df <- data.frame(
          specificity = c(partial.auc[2], partial.auc[2], partial.auc[1], partial.auc[1]),
          sensitivity = c(0, one, one, 0)
        )
      }
      df <- ggroc_add_one_or_hundred(df, percent)
      aes <- ggroc_add_mapping(ggroc_xy_aes(plot), mapping)
      layer <- do.call(
        ggplot2::geom_polygon,
        c(
          list(mapping = aes, data = df),
          extras,
          list(inherit.aes = FALSE)
        )
      )
      list(ggroc_theme_defaults(layer, if (is.null(extras$fill) && is.null(mapping$fill)) list(fill = ggroc_theme_mix(0.8))))
    }
  )
}

geom_polygon_max_auc.roc <- function(data, ...) {
  roc_utils_stop_if_no_auc(data, "data")
  geom_polygon_max_auc(data$auc, ...)
}

geom_polygon_max_auc.smooth.roc <- geom_polygon_max_auc.roc

geom_polygon_max_auc.list <- function(data, ...) {
  ggroc_layer_list(data, geom_polygon_max_auc, ...)
}
