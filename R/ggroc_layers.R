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

# Theme-based default colours, as ggplot2's own geoms have: 'defaults' is a
# named list of from_theme() expressions that replace those of the layer's
# geom. In the geom rather than the mapping, so that they stay out of the
# scales and legend titles.
ggroc_theme_defaults <- function(layer, defaults) {
  if (length(defaults) > 0) {
    # ggproto() finds its parent lazily: take it out of the layer first
    geom <- layer$geom
    default_aes <- geom$default_aes
    default_aes[names(defaults)] <- defaults
    layer$geom <- ggplot2::ggproto(NULL, geom, default_aes = default_aes)
  }
  layer
}

# 'amount' of the way from the theme's ink to its paper
ggroc_theme_mix <- function(amount) {
  bquote(ggplot2::from_theme(scales::col_mix(ink, paper, .(amount))))
}

# The layer's own mapping plus the user's 'mapping', which wins
ggroc_add_mapping <- function(mapping, user) {
  for (aesthetic in names(user)) {
    mapping[[aesthetic]] <- user[[aesthetic]]
  }
  mapping
}

# The default colour of the polygons and error bars of a curve: the colour
# that ggroc maps to the curves of a list ('named' layers only), otherwise
# the fixed colour of the ROC line, otherwise the theme's ink. A colour
# given by the user, set in 'extras' or mapped in 'mapping', wins. Fills
# ('aesthetic = "fill"') take the colour with 'alpha'.
# Returns the mapping to add, the updated extras and the theme defaults
# for ggroc_theme_defaults().
ggroc_curve_colour <- function(plot, named, extras, mapping, aesthetic = c("colour", "fill"), alpha = 0.25) {
  aesthetic <- match.arg(aesthetic)
  out <- list(mapping = list(), extras = extras, defaults = list())
  if (!is.null(extras[[aesthetic]]) || !is.null(mapping[[aesthetic]])) {
    return(out)
  }
  plot_colour <- ggroc_plot_mapping(plot)$colour
  if (named && !is.null(plot_colour) && rlang::as_label(plot_colour) == "name") {
    .data <- rlang::.data
    if (aesthetic == "colour") {
      out$mapping$colour <- ggplot2::expr(.data$name)
    } else {
      # The fill is computed from the mapped colour, which then becomes the
      # border colour of 'extras'. A colour set in 'extras' would drop the
      # mapping.
      border <- if (is.null(extras$colour)) NA else extras$colour
      out$mapping$colour <- bquote(ggplot2::stage(.data$name, after_scale = .(border)))
      out$mapping$fill <- bquote(ggplot2::after_scale(ggplot2::alpha(colour, .(alpha))))
      out$extras$colour <- NULL
    }
    if (is.null(extras$show.legend)) {
      out$extras$show.legend <- FALSE
    }
  } else {
    roc_i <- ggroc_roc_line_index(plot@layers)
    line_colour <- if (roc_i <= length(plot@layers)) plot@layers[[roc_i]]$aes_params$colour
    if (aesthetic == "colour") {
      out$extras$colour <- line_colour
    } else if (is.null(line_colour)) {
      out$defaults$fill <- bquote(ggplot2::from_theme(ggplot2::alpha(ink, .(alpha))))
    } else {
      out$extras$fill <- ggplot2::alpha(line_colour, alpha)
    }
  }
  out
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

# The layers use ggplot2 >= 4.0.0 (S7 objects, `plot@layers`). Suggests
# versions are not enforced by requireNamespace(), so check here.
ggroc_ggplot2_version <- function() {
  getNamespaceVersion("ggplot2")
}

ggplot_add.ggroc_layer <- function(object, plot, ...) {
  version <- ggroc_ggplot2_version()
  if (package_version(version) < "4.0.0") {
    stop(sprintf("pROC's ggplot2 layers (geom_*) require ggplot2 >= 4.0.0, but ggplot2 %s is loaded.", version), call. = FALSE)
  }
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

# List method shared by the layers: one layer per curve, named as in
# ggroc.list, with a 'name' column so that facet_*(~ name) draws each layer
# in its own panel only. make_layers(plot, named = TRUE) tells the layer
# that it will have the 'name' column.
ggroc_layer_list <- function(data, layer_fun, ...) {
  data <- ggroc_name_list(data)
  layers <- lapply(data, layer_fun, ...)
  new_ggroc_layer(
    behind = any(vapply(layers, function(l) isTRUE(l$behind), logical(1))),
    make_layers = function(plot) {
      unlist(lapply(names(layers), function(name) {
        lapply(layers[[name]]$make_layers(plot, named = TRUE), function(layer) {
          layer$data$name <- factor(name, levels = names(layers))
          layer
        })
      }), recursive = FALSE)
    }
  )
}
