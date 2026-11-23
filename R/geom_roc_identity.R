geom_roc_identity <- function(data, ...) {
  UseMethod("geom_roc_identity")
}

geom_roc_identity.roc <- function(data, colour = NULL, mapping = NULL, ...) {
  load.ggplot2()
  extras <- size.to.linewidth(list(...))
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
      .data <- rlang::.data
      layer <- do.call(
        ggplot2::geom_segment,
        c(
          list(
            data = df,
            mapping = ggroc_add_mapping(ggplot2::aes(
              x = .data$x,
              y = .data$y,
              xend = .data$xend,
              yend = .data$yend
            ), mapping),
            inherit.aes = FALSE
          ),
          extras
        )
      )
      list(ggroc_theme_defaults(layer, if (is.null(extras$colour) && is.null(mapping$colour)) list(colour = ggroc_theme_mix(0.66))))
    }
  )
}

geom_roc_identity.smooth.roc <- geom_roc_identity.roc
