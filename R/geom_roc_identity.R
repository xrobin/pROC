geom_roc_identity <- function(data, ...) {
  UseMethod("geom_roc_identity")
}

geom_roc_identity.roc <- function(data, colour = "darkgrey", ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  if (is.null(extras$colour)) {
    extras$colour <- colour
  }
  percent <- data$percent
  new_ggroc_layer(
    behind = TRUE,
    make_layers = function(plot) {
      one <- if (percent) 100 else 1
      if (ggroc_legacy_axes_from_plot(plot)) {
        df <- data.frame(x = 0, xend = one, y = 0, yend = one)
      } else {
        df <- data.frame(x = one, xend = 0, y = 0, yend = one)
      }
      list(do.call(
        ggplot2::geom_segment,
        c(
          list(
            data = df,
            mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
            inherit.aes = FALSE
          ),
          extras
        )
      ))
    }
  )
}

geom_roc_identity.smooth.roc <- geom_roc_identity.roc
