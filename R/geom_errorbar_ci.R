geom_errorbar_ci <- function(data, ...) {
  UseMethod("geom_errorbar_ci")
}

geom_errorbar_ci.ci.se <- function(data, ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  roc <- attr(data, "roc")
  percent <- roc$percent
  if (is.null(extras$width)) {
    extras$width <- ggroc_errorbar_width(percent)
  }
  sp <- as.numeric(attr(data, "specificities"))
  df <- data.frame(
    specificity = sp,
    ymin = data[, 1],
    ymax = data[, 3]
  )
  df <- ggroc_add_one_or_hundred(df, percent)
  new_ggroc_layer(
    behind = FALSE,
    make_layers = function(plot) {
      xcol <- ggroc_x_col(plot)
      .data <- rlang::.data
      mapping <- ggplot2::aes(
        x = .data[[xcol]],
        ymin = .data$ymin,
        ymax = .data$ymax
      )
      list(do.call(
        ggplot2::geom_errorbar,
        c(
          list(mapping = mapping, data = df),
          extras,
          list(inherit.aes = FALSE)
        )
      ))
    }
  )
}

geom_errorbar_ci.ci.sp <- function(data, ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  roc <- attr(data, "roc")
  percent <- roc$percent
  one <- if (percent) 100 else 1
  if (is.null(extras$width)) {
    extras$width <- ggroc_errorbar_width(percent)
  }
  se <- as.numeric(attr(data, "sensitivities"))
  df <- data.frame(
    sensitivity = se,
    xmin = data[, 1],
    xmax = data[, 3]
  )
  new_ggroc_layer(
    behind = FALSE,
    make_layers = function(plot) {
      bars <- df
      if (ggroc_legacy_axes_from_plot(plot)) {
        tmp <- one - bars$xmax
        bars$xmax <- one - bars$xmin
        bars$xmin <- tmp
      }
      .data <- rlang::.data
      mapping <- ggplot2::aes(
        y = .data$sensitivity,
        xmin = .data$xmin,
        xmax = .data$xmax
      )
      list(do.call(
        ggplot2::geom_errorbar,
        c(
          list(mapping = mapping, data = bars),
          extras,
          list(inherit.aes = FALSE, orientation = "y")
        )
      ))
    }
  )
}

geom_errorbar_ci.ci.thresholds <- function(data, ...) {
  load.ggplot2()
  extras <- list(...)
  names(extras) <- sub("color", "colour", names(extras))
  roc <- attr(data, "roc")
  percent <- roc$percent
  one <- if (percent) 100 else 1
  if (is.null(extras$width)) {
    extras$width <- ggroc_errorbar_width(percent)
  }
  df <- data.frame(
    specificity = data$specificity[, 2],
    sensitivity = data$sensitivity[, 2],
    sp_lower = data$specificity[, 1],
    sp_upper = data$specificity[, 3],
    se_lower = data$sensitivity[, 1],
    se_upper = data$sensitivity[, 3]
  )
  df <- ggroc_add_one_or_hundred(df, percent)
  new_ggroc_layer(
    behind = FALSE,
    make_layers = function(plot) {
      xcol <- ggroc_x_col(plot)
      .data <- rlang::.data
      if (ggroc_legacy_axes_from_plot(plot)) {
        tmp <- one - df$sp_upper
        df$h_xmax <- one - df$sp_lower
        df$h_xmin <- tmp
      } else {
        df$h_xmin <- df$sp_lower
        df$h_xmax <- df$sp_upper
      }
      v_aes <- ggplot2::aes(
        x = .data[[xcol]],
        ymin = .data$se_lower,
        ymax = .data$se_upper
      )
      h_aes <- ggplot2::aes(
        y = .data$sensitivity,
        xmin = .data$h_xmin,
        xmax = .data$h_xmax
      )
      list(
        do.call(
          ggplot2::geom_errorbar,
          c(
            list(mapping = v_aes, data = df),
            extras,
            list(inherit.aes = FALSE)
          )
        ),
        do.call(
          ggplot2::geom_errorbar,
          c(
            list(mapping = h_aes, data = df),
            extras,
            list(inherit.aes = FALSE, orientation = "y")
          )
        )
      )
    }
  )
}
