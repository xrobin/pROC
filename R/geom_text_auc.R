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
  ggroc_text_layer(data.frame(x = x, y = y, label = label, vjust = vjust), hjust, extras)
}

geom_text_auc.roc <- function(data, ci = NULL, ...) {
  roc_utils_stop_if_no_auc(data, "data")
  geom_text_auc(data$auc, ci = ci, ...)
}

geom_text_auc.smooth.roc <- geom_text_auc.roc

geom_text_auc.ci.auc <- function(data, ...) {
  geom_text_auc(attr(data, "auc"), ci = data, ...)
}

# One text layer for all the curves, so that the labels sharing a panel
# can be stacked
geom_text_auc.list <- function(data, ci = NULL, pattern = NULL, name = FALSE, ...) {
  if (!is.null(ci)) {
    stop("'ci' cannot be shared by the curves of a list: pass a list of 'ci.auc' objects as 'data' instead.")
  }
  if (length(data) == 0) {
    return(new_ggroc_layer(make_layers = function(plot) list()))
  }
  data <- ggroc_name_list(data)
  texts <- lapply(data, function(d) geom_text_auc(d, pattern = pattern, ...)$text)
  df <- do.call(rbind, lapply(texts, `[[`, "df"))
  df$name <- factor(names(data), levels = names(data))
  if (name) {
    df$label <- paste(names(data), df$label)
  }
  ggroc_text_layer(df, texts[[1]]$hjust, texts[[1]]$extras)
}

# The text layer of geom_text_auc. The labels that share a panel are stacked
# 1.5 lines apart below the first one, so that facets show one label per
# panel at (x, y) and a plot without facets shows them all.
# Unless 'colour' is given, the labels take the colour of their curve, as
# plot.roc(print.auc = TRUE) does: the ggroc colour mapping of the curves
# for a list, or the fixed colour of the ROC line.
ggroc_text_layer <- function(df, hjust, extras) {
  layer <- new_ggroc_layer(
    behind = FALSE,
    make_layers = function(plot) {
      stat <- ggplot2::ggproto("StatRocStackText", ggplot2::Stat,
        compute_panel = function(data, scales) {
          if (nrow(data) > 1) {
            if (!is.numeric(data$vjust)) {
              stop("Labels sharing a panel are stacked from a numeric 'vjust'.")
            }
            data$vjust <- data$vjust + 1.5 * (seq_len(nrow(data)) - 1)
          }
          data
        }
      )
      .data <- rlang::.data
      mapping <- ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, vjust = .data$vjust)
      if (is.null(extras$colour)) {
        plot_colour <- ggroc_plot_mapping(plot)$colour
        if (!is.null(df$name) && !is.null(plot_colour) && rlang::as_label(plot_colour) == "name") {
          mapping$colour <- ggplot2::expr(.data$name)
          if (is.null(extras$show.legend)) {
            extras$show.legend <- FALSE
          }
        } else {
          roc_i <- ggroc_roc_line_index(plot@layers)
          if (roc_i <= length(plot@layers)) {
            extras$colour <- plot@layers[[roc_i]]$aes_params$colour
          }
        }
      }
      list(do.call(
        ggplot2::geom_text,
        c(
          list(mapping = mapping, data = df, stat = stat, hjust = hjust, inherit.aes = FALSE),
          extras
        )
      ))
    }
  )
  layer$text <- list(df = df, hjust = hjust, extras = extras)
  layer
}
