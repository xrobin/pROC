geom_polygon_max_auc <- function(data, ...) {
  UseMethod("geom_polygon_max_auc")
}

geom_polygon_max_auc.auc <- function(data, fill = "#EEEEEE", colour = NA, ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  if (is.null(extras$fill)) {
    extras$fill <- fill
  }
  if (is.null(extras$colour)) {
    extras$colour <- colour
  }
  percent <- attr(data, "percent")
  partial.auc <- attr(data, "partial.auc")
  partial.auc.focus <- attr(data, "partial.auc.focus")
  new_ggroc_layer(
    behind = TRUE,
    make_layers = function(plot) {
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
      aes <- ggroc_xy_aes(plot)
      list(do.call(
        ggplot2::geom_polygon,
        c(
          list(mapping = aes, data = df),
          extras,
          list(inherit.aes = FALSE)
        )
      ))
    }
  )
}

geom_polygon_max_auc.roc <- function(data, ...) {
  geom_polygon_max_auc(data$auc, ...)
}

geom_polygon_max_auc.smooth.roc <- geom_polygon_max_auc.roc
