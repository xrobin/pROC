#!/usr/bin/env Rscript
## Step 3 — check every reverse dependency twice: against the CRAN release of
## pROC (old) and against this working tree (new).
##
## Both sides get byte-identical configuration. That is what makes the
## comparison mean anything: absolute check results are noisy, the delta is not.
##
## Two traps this guards against, both of which invalidated an earlier run:
##
##  1. reverse = <character vector>. If you pass a list(which=...) instead,
##     check_packages_in_dir recomputes the reverse dependencies of *every*
##     tarball in the directory. On a restart the revdeps are already sitting
##     there, so 206 packages silently becomes 2138.
##  2. pROC must not be visible anywhere except the side under test. The
##     EasyBuild R-bundle-CRAN ships one, and the shared library holds one
##     after step 2 — either will be picked up in preference to the tarball.

options(warn = 1)
WORK <- Sys.getenv("REVDEP_WORK", unset = file.path("/scratch", Sys.getenv("USER"), "pROC_revdeps"))
LIB  <- Sys.getenv("REVDEP_LIB",  unset = file.path(WORK, "Library"))
NCPUS <- as.integer(Sys.getenv("NCPUS", "8"))
SIDES <- strsplit(Sys.getenv("REVDEP_SIDES", "old,new"), ",")[[1]]
req  <- readRDS(file.path(WORK, "requirements.rds"))
## REVDEP_ONLY restricts the run to named packages. Use it to smoke-test the
## pipeline before committing to the full set, or to re-check one regression.
only <- Sys.getenv("REVDEP_ONLY")
if (nzchar(only)) {
  sel <- trimws(strsplit(only, ",")[[1]])
  unknown <- setdiff(sel, req$tocheck)
  if (length(unknown)) stop("not reverse dependencies: ", paste(unknown, collapse = ", "))
  req$tocheck <- sel
  message(sprintf("REVDEP_ONLY: restricted to %d package(s): %s",
                  length(sel), paste(sel, collapse = ", ")))
}
## Locate the package source tree. Searching upwards rather than assuming a
## fixed depth means these scripts keep working wherever they are placed, and
## the presence of man/ distinguishes a source checkout from an installed copy
## (which has help/ instead). REVDEP_PKG_ROOT overrides it.
find_pkg_root <- function(start) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:8) {
    if (file.exists(file.path(d, "DESCRIPTION")) && dir.exists(file.path(d, "man"))) return(d)
    p <- dirname(d)
    if (identical(p, d)) break
    d <- p
  }
  stop("cannot find the package source root above ", start,
       " -- set REVDEP_PKG_ROOT to the pROC checkout")
}
PKG_ROOT <- Sys.getenv("REVDEP_PKG_ROOT")
if (!nzchar(PKG_ROOT)) {
  PKG_ROOT <- find_pkg_root(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
}
options(repos = c(req$repos, req$extra_repos), timeout = 900)

## --- pROC must not be on the shared library path -----------------------
shared_proc <- file.path(LIB, req$pkg)
if (file.exists(shared_proc)) {
  keep <- file.path(WORK, sprintf("removed_%s_from_shared_Library", req$pkg))
  unlink(keep, recursive = TRUE)
  ok <- file.rename(shared_proc, keep)
  message(sprintf("moved %s out of the shared library (%s)", req$pkg,
                  if (ok) "ok" else "FAILED"))
}
stopifnot(!file.exists(shared_proc))

## --- the two tarballs --------------------------------------------------
stamp   <- format(Sys.time(), "%Y%m%d-%H%M%S")
run_dir <- file.path(WORK, "runs", paste0("run-", stamp))
dirs <- setNames(file.path(run_dir, SIDES), SIDES)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
message("run directory: ", run_dir)

if ("old" %in% SIDES) {
  message("downloading CRAN ", req$pkg, " ", req$cran_version, " ...")
  dl <- download.packages(req$pkg, destdir = dirs[["old"]], type = "source",
                          repos = req$repos[["CRAN"]])
  stopifnot(nrow(dl) == 1L)
}
if ("new" %in% SIDES) {
  message("building ", req$pkg, " ", req$local_version, " from ", PKG_ROOT, " ...")
  xvfb <- Sys.getenv("XVFB_SERVERARGS")
  build <- if (nzchar(Sys.which("xvfb-run")))
    sprintf("xvfb-run -a %s R CMD build %s",
            if (nzchar(xvfb)) sprintf("-s %s", shQuote(xvfb)) else "", shQuote(PKG_ROOT))
  else sprintf("R CMD build %s", shQuote(PKG_ROOT))
  owd <- setwd(dirs[["new"]]); on.exit(setwd(owd), add = TRUE)
  writeLines(system(paste(build, "2>&1"), intern = TRUE))
  setwd(owd)
  tb <- Sys.glob(file.path(dirs[["new"]], sprintf("%s_*.tar.gz", req$pkg)))
  stopifnot(length(tb) == 1L)
  message("built: ", basename(tb))
}

## --- identical check configuration on both sides -----------------------
## FORCE_SUGGESTS false: a Suggests that is not on CRAN would otherwise abort
##   the check at "checking package dependencies" with an ERROR, so the package
##   would tell us nothing. CRAN's own machines set this for the same reason.
## DONTTEST true: \donttest examples are frequently where a package actually
##   calls into pROC. They are the most informative part of the check, so they
##   are bounded with timeouts rather than switched off.
check_env <- c(
  "_R_CHECK_FORCE_SUGGESTS_=false",
  "_R_CHECK_DONTTEST_EXAMPLES_=true",
  sprintf("_R_CHECK_ELAPSED_TIMEOUT_=%s",          Sys.getenv("REVDEP_TIMEOUT",       "45m")),
  sprintf("_R_CHECK_INSTALL_ELAPSED_TIMEOUT_=%s",  Sys.getenv("REVDEP_TIMEOUT_INST",  "20m")),
  sprintf("_R_CHECK_EXAMPLES_ELAPSED_TIMEOUT_=%s", Sys.getenv("REVDEP_TIMEOUT_EX",    "20m")),
  sprintf("_R_CHECK_TESTS_ELAPSED_TIMEOUT_=%s",    Sys.getenv("REVDEP_TIMEOUT_TESTS", "20m")),
  sprintf("_R_CHECK_ONE_TEST_ELAPSED_TIMEOUT_=%s", Sys.getenv("REVDEP_TIMEOUT_ONE",   "10m")),
  sprintf("_R_CHECK_ONE_VIGNETTE_ELAPSED_TIMEOUT_=%s", Sys.getenv("REVDEP_TIMEOUT_VIG", "10m")),
  "_R_CHECK_BUILD_VIGNETTES_ELAPSED_TIMEOUT_=20m")

## check_args/check_env as a 2-element list: [[1]] for the package itself,
## [[2]] for its reverse dependencies. --as-cran on the revdeps would add
## policy NOTEs unrelated to pROC; they would cancel in the diff but cost time.
##
## --no-manual on the reverse dependencies: "checking PDF version of manual"
## tests the formatting of each package's OWN Rd files, which nothing in pROC
## can affect, so it is pure cost. Vignettes ARE built: that is where reverse
## dependencies actually exercise pROC. LaTeX for them comes from TinyTeX (see
## 00_env.sh); without it, 11 of 213 packages failed on vignettes alone.
check_args <- list("--as-cran", "--no-manual")

## ---- X display -------------------------------------------------------
## This R has no cairo bitmap devices (grDevices C_cairoProps(2) is FALSE),
## so png() can only go through X11, and capabilities("X11") is therefore a
## RUNTIME question: it is TRUE only when DISPLAY points at a live server.
## Any check whose examples or vignettes draw a PNG -- which includes every
## knitr .Rmd vignette with a plot -- fails with "unable to start device PNG"
## if the display is missing.
##
## check_packages_in_dir(xvfb=) starts a server itself, but it proved
## unreliable here: it leaves the server running when the R process is killed
## (a stale one survived 19 days), and a run produced X11 failures in 4
## packages on one side and 9 on the other. That asymmetry looks exactly like
## a regression and is not one. So start one server, PROVE a PNG can be drawn
## through it, and hand the same DISPLAY to both sides explicitly.
xvfb_args <- Sys.getenv("XVFB_SERVERARGS", "-screen 0 1280x1024x24")

start_xvfb <- function() {
  for (attempt in 1:40) {
    n <- 90L + attempt
    if (file.exists(sprintf("/tmp/.X%d-lock", n))) next
    ## Capture the pid so it can be killed precisely at the end. Killing by
    ## pattern is both unreliable and dangerous, and leaving it running leaks
    ## a server: two were found still alive days after their run was killed.
    pid <- suppressWarnings(system(
      sprintf("Xvfb :%d %s >/dev/null 2>&1 & echo $!", n, xvfb_args), intern = TRUE))
    Sys.sleep(3)
    ok <- system2("Rscript", c("--vanilla", "-e", shQuote(
        'q(status = if (capabilities("X11") && !inherits(try({png(tempfile(fileext=".png")); plot(1); dev.off()}, silent=TRUE), "try-error")) 0L else 1L)')),
      env = sprintf("DISPLAY=:%d", n), stdout = NULL, stderr = NULL)
    if (identical(ok, 0L)) {
      message(sprintf("Xvfb on :%d (pid %s), png() verified", n, pid))
      return(list(display = n, pid = pid))
    }
    if (length(pid)) try(system2("kill", pid, stdout = NULL, stderr = NULL), silent = TRUE)
    message(sprintf("  display :%d did not work, trying another", n))
  }
  stop("could not start an Xvfb on which png() works")
}

xserver <- start_xvfb()
check_env <- c(check_env, sprintf("DISPLAY=:%d", xserver$display))
## on.exit() is useless here: at the top level of a script there is no function
## to exit from, so it never runs and the server is leaked. Shut it down
## explicitly, and from an on.exit INSIDE a function for the error path.
stop_xserver <- function() {
  try(system2("kill", xserver$pid, stdout = NULL, stderr = NULL), silent = TRUE)
  unlink(sprintf("/tmp/.X%d-lock", xserver$display))
}

run_sides <- function() {
on.exit(stop_xserver(), add = TRUE)
for (side in SIDES) {
  message(sprintf("\n=== checking %s side: %d reverse dependencies ===", side, length(req$tocheck)))
  message(format(Sys.time()))
  res <- tools::check_packages_in_dir(
    dir        = dirs[[side]],
    reverse    = req$tocheck,          # explicit vector; never recomputed
    check_args = check_args,
    check_env  = list(check_env, check_env),
    Ncpus      = NCPUS,
    clean      = FALSE,
    xvfb       = FALSE)   # we start and verify our own, above
  saveRDS(res, file.path(run_dir, sprintf("result_%s.rds", side)))
  message(sprintf("=== %s side finished %s ===", side, format(Sys.time())))
}
}
run_sides()
stop_xserver()
writeLines(run_dir, file.path(WORK, "last_run_dir.txt"))
message("\nrun directory: ", run_dir)
