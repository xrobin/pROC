layer_geom_class <- function(layer) {
  class(layer$geom)[[1]]
}

layer_geom_classes <- function(plot) {
  unname(vapply(plot@layers, layer_geom_class, character(1)))
}
