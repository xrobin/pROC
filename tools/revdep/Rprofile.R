# Headless node: Xvfb has no XLFD Helvetica, so base graphics fail under
# bitmapType="Xlib". Use cairo/fontconfig instead.
options(bitmapType = "cairo")
if (capabilities("cairo")) try(grDevices::X11.options(type = "cairo"), silent = TRUE)
