#!/usr/bin/env Rscript
## Install the LaTeX packages a run turned out to need.
##
## TinyTeX ships a minimal TeX Live, so vignettes fail one missing .sty at a
## time. This scans the check logs of a run, maps the missing files to TeX Live
## packages and installs them. Needs no root: TinyTeX lives in $REVDEP_WORK.
##   usage: install_missing_tex.R [run_dir]
## Re-run step 3 afterwards; repeat until it reports nothing missing.

WORK <- Sys.getenv("REVDEP_WORK", unset = file.path("/scratch", Sys.getenv("USER"), "pROC_revdeps"))
run_dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else
  readLines(file.path(WORK, "last_run_dir.txt"))[1]
.libPaths(c(file.path(WORK, "Library"), .Library))

logs <- Sys.glob(file.path(run_dir, c("old", "new"), "*.Rcheck", "00check.log"))
if (!length(logs)) stop("no check logs under ", run_dir)
txt <- unlist(lapply(logs, readLines, warn = FALSE))
hits <- regmatches(txt, regexpr("File `[^']+' not found", txt))
files <- unique(sub("^File `", "", sub("' not found$", "", hits)))

## Missing font metrics are reported differently and would otherwise be missed:
##   ! Font U/bbm/m/n/10.95=bbm10 at 10.95pt not loadable: Metric (TFM) file not found.
## The .sty and the fonts are often separate TeX Live packages (bbm-macros vs bbm).
fonts <- regmatches(txt, regexpr("Font [A-Z]+/[^/]+/[^ ]*=([A-Za-z0-9]+)", txt))
fonts <- unique(sub(".*=", "", fonts))
if (length(fonts)) {
  message(sprintf("missing font metrics: %s", paste(fonts, collapse = ", ")))
  files <- c(files, paste0(fonts, ".tfm"))
}

if (!length(files)) {
  message("no missing LaTeX files in ", run_dir)
  quit(save = "no")
}
message(sprintf("missing: %s", paste(files, collapse = ", ")))
pkgs <- tinytex::parse_packages(files = files, quiet = c(TRUE, TRUE, TRUE))
if (!length(pkgs)) {
  message("could not map those files to TeX Live packages; install by hand")
  quit(save = "no", status = 1)
}
message(sprintf("installing: %s", paste(pkgs, collapse = ", ")))
tinytex::tlmgr_install(pkgs)

message("done. Re-run step 3, then this again, until it reports nothing missing.")
