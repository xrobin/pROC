geom_polygon_ci <- function(data, ...) {
  UseMethod("geom_polygon_ci")
}

geom_polygon_ci.ci.se <- function(data, fill = "gainsboro", colour = NA, ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  if (is.null(extras$fill)) {
    extras$fill <- fill
  }
  if (is.null(extras$colour)) {
    extras$colour <- colour
  }
  roc <- attr(data, "roc")
  percent <- roc$percent
  one <- if (percent) 100 else 1
  sp <- as.numeric(attr(data, "specificities"))
  if (length(data[, 1]) < 15) {
    warning("Low definition shape.")
  }
  df <- data.frame(
    specificity = c(0, sp, one, rev(sp)),
    sensitivity = c(one, data[, 1], 0, rev(data[, 3]))
  )
  df <- ggroc_add_one_or_hundred(df, percent)
  new_ggroc_layer(
    behind = TRUE,
    make_layers = function(plot) {
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

geom_polygon_ci.ci.sp <- function(data, fill = "gainsboro", colour = NA, ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  if (is.null(extras$fill)) {
    extras$fill <- fill
  }
  if (is.null(extras$colour)) {
    extras$colour <- colour
  }
  roc <- attr(data, "roc")
  percent <- roc$percent
  one <- if (percent) 100 else 1
  se <- as.numeric(attr(data, "sensitivities"))
  if (length(data[, 1]) < 15) {
    warning("Low definition shape.")
  }
  df <- data.frame(
    specificity = c(one, data[, 1], 0, rev(data[, 3])),
    sensitivity = c(0, se, one, rev(se))
  )
  df <- ggroc_add_one_or_hundred(df, percent)
  new_ggroc_layer(
    behind = TRUE,
    make_layers = function(plot) {
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
