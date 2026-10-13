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
