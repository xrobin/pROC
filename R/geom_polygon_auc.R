geom_polygon_auc <- function(data, ...) {
  UseMethod("geom_polygon_auc")
}

# Extra vertices that close the AUC polygon. Set columns by name so empirical
# coords (with threshold) and smooth.roc coords (without) both close to (0, 0).
ggroc_auc_polygon_close <- function(df, specificity, sensitivity, one) {
  extra <- df[rep(NA_integer_, length(specificity)), , drop = FALSE]
  extra[["specificity"]] <- specificity
  extra[["sensitivity"]] <- sensitivity
  extra[["1-specificity"]] <- one - specificity
  rbind(df, extra)
}

geom_polygon_auc.auc <- function(data, ...) {
  load.ggplot2()
  extras <- list(...)
  roc <- attr(data, "roc")
  roc$auc <- data
  new_ggroc_layer(
    behind = TRUE,
    make_layers = function(plot) {
      la <- ggroc_legacy_axes_from_plot(plot)
      df <- get.coords.for.ggplot(roc, ignore.partial.auc = FALSE)

      partial.auc <- attr(data, "partial.auc")
      one.or.hundred <- ifelse(attr(data, "percent"), 100, 1)
      if (identical(partial.auc, FALSE)) {
        df <- ggroc_auc_polygon_close(df, 0, 0, one.or.hundred)
      } else if (attr(data, "partial.auc.focus") == "sensitivity") {
        df <- ggroc_auc_polygon_close(df, c(0, 0), partial.auc, one.or.hundred)
      } else {
        df <- ggroc_auc_polygon_close(df, rev(partial.auc), c(0, 0), one.or.hundred)
      }

      aes <- get.aes.for.ggplot(roc, la)
      list(do.call(
        ggplot2::geom_polygon,
        c(list(mapping = aes$aes, data = df), extras, list(inherit.aes = FALSE))
      ))
    }
  )
}

geom_polygon_auc.roc <- function(data, ...) {
  geom_polygon_auc(data$auc, ...)
}

geom_polygon_auc.smooth.roc <- geom_polygon_auc.roc
