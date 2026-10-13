# Capture the graphics primitives drawn by an expression, from the display list
# of a fresh null pdf device. Returns a list of list(name = "text", args = list(...)).
capture_display_list <- function(expr) {
  pdf(NULL)
  on.exit(dev.off())
  dev.control("enable")
  force(expr)
  lapply(recordPlot()[[1]], function(e) {
    list(name = sub("^C_", "", e[[2]][[1]]$name), args = as.list(e[[2]])[-1])
  })
}

display_list_calls <- function(dl, name) {
  Filter(function(e) e$name == name, dl)
}

# Labels drawn with text()
display_list_texts <- function(dl) {
  unlist(lapply(display_list_calls(dl, "text"), function(e) e$args[[2]]))
}

# Vertices of the polygons drawn with polygon()
display_list_polygons <- function(dl) {
  lapply(display_list_calls(dl, "polygon"), function(e) list(x = e$args[[1]], y = e$args[[2]]))
}

# Labels passed to axis(), one element per axis() call
display_list_axis_labels <- function(dl) {
  lapply(display_list_calls(dl, "axis"), function(e) e$args[[3]])
}
