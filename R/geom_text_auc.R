geom_text_auc <- function(data, ...) {
  UseMethod("geom_text_auc")
}

geom_text_auc.auc <- function(data, ci = NULL, x = NULL, y = NULL, pattern = NULL, hjust = 0, vjust = 1, mapping = NULL, ...) {
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
  label <- roc_utils_auc_label(data, ci, pattern)
  ggroc_text_layer(data.frame(x = x, y = y, label = label, vjust = vjust), hjust, extras, mapping)
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
  ggroc_text_layer(df, texts[[1]]$hjust, texts[[1]]$extras, texts[[1]]$mapping)
}

# The text layer of geom_text_auc. The labels that share a panel are stacked
# 1.5 lines apart below the first one, so that facets show one label per
# panel at (x, y) and a plot without facets shows them all.
ggroc_text_layer <- function(df, hjust, extras, mapping) {
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
      aes <- ggroc_add_mapping(
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, vjust = .data$vjust),
        mapping
      )
      list(do.call(
        ggplot2::geom_text,
        c(
          list(mapping = aes, data = df, stat = stat, inherit.aes = FALSE),
          if (is.null(mapping$hjust)) list(hjust = hjust),
          extras
        )
      ))
    }
  )
  layer$text <- list(df = df, hjust = hjust, extras = extras, mapping = mapping)
  layer
}
