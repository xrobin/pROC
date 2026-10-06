#!/usr/bin/env Rscript
## Step 2 — make the shared library hold every package the check needs,
## at a current version, and actually loadable.
##
## "Installed" and "usable" are different claims. A concurrent or interrupted
## install leaves packages that are present in the library yet fail to load,
## and R CMD check then reports that as the dependant package's problem. So
## every package is load-tested in a fresh process, and anything that fails is
## reinstalled. The cycle repeats until it stops improving.
##
## pROC itself is deliberately NOT pinned here: it is present during installs
## (revdeps that Depend on it need it), and step 3 removes it from the shared
## library so each side can supply its own.

options(warn = 1)
WORK  <- Sys.getenv("REVDEP_WORK", unset = file.path("/scratch", Sys.getenv("USER"), "pROC_revdeps"))
LIB   <- Sys.getenv("REVDEP_LIB",  unset = file.path(WORK, "Library"))
NCPUS <- as.integer(Sys.getenv("NCPUS", "8"))
MAXIT <- as.integer(Sys.getenv("REVDEP_MAX_PASSES", "3"))
req   <- readRDS(file.path(WORK, "requirements.rds"))
dir.create(LIB, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(LIB, .Library))
options(repos = c(req$repos, req$extra_repos), Ncpus = NCPUS, timeout = 900)

## Packages that are symlinks into a shared EasyBuild bundle rather than real
## installations. These must be rebuilt. A bundle's compiled packages are built
## against the bundle's own versions of their dependencies; as soon as this
## library updates one of those dependencies, the bundle's binaries can break
## with an undefined symbol and there is no way to predict which. (Observed:
## updating RcppParallel 5.1.11 -> 6.2.1 removed tbb::task from the bundled
## TBB, which broke the bundle's Rfast.so and every install that loaded it.)
symlinked <- function() {
  d <- list.files(LIB, full.names = TRUE)
  basename(d[!is.na(Sys.readlink(d)) & nzchar(Sys.readlink(d))])
}

audit <- function() {
  ip <- installed.packages(lib.loc = c(LIB, .Library), noCache = TRUE)
  ip <- ip[!duplicated(rownames(ip)), , drop = FALSE]
  present <- intersect(req$needed, rownames(ip))
  missing <- setdiff(req$needed, rownames(ip))
  ## Look the wanted version up positionally. ifelse() would evaluate the
  ## whole ap_all[present, ] subscript regardless of the condition, which
  ## errors as soon as a package came from the Archive or another repository
  ## and so has no row there.
  have <- ip[present, "Version"]
  want <- rep(NA_character_, length(present)); names(want) <- present
  inap <- intersect(present, rownames(req$ap_all))
  want[inap] <- req$ap_all[inap, "Version"]
  cmp <- !is.na(want)
  outdated <- present[cmp][package_version(have[cmp]) < package_version(want[cmp])]
  linked <- intersect(present, symlinked())
  list(present = present, missing = missing, outdated = outdated, linked = linked)
}

load_test <- function(pkgs) {
  if (!length(pkgs)) return(character())
  message(sprintf("  load-testing %d packages ...", length(pkgs)))
  one <- function(p) {
    r <- system2("Rscript", c("--vanilla", "-e", shQuote(sprintf(
      '.libPaths(c("%s", .Library)); suppressMessages(loadNamespace("%s"))', LIB, p))),
      stdout = TRUE, stderr = TRUE)
    st <- attr(r, "status")
    if (is.null(st) || st == 0L) NA_character_ else p
  }
  res <- unlist(parallel::mclapply(pkgs, one, mc.cores = min(NCPUS + 4L, 16L)))
  res[!is.na(res)]
}

archive_tarball <- function(p) {
  idx <- tryCatch(readLines(sprintf("https://cran.r-project.org/src/contrib/Archive/%s/", p),
                            warn = FALSE), error = function(e) character())
  tars <- sort(unlist(unique(regmatches(idx, gregexpr(sprintf("%s_[0-9.-]+\\.tar\\.gz", p), idx)))))
  if (length(tars)) sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s", p, tail(tars, 1)) else NA_character_
}

install_archived <- function(pkgs, depth = 2L) {
  ## Archived packages are reachable only by explicit tarball URL, and their
  ## own dependencies are often archived too (hdi needs scalreg, LPKsample
  ## needs LPGraph). Resolve one or two levels of that, then give up: these
  ## are Suggests of a single reverse dependency, and a failure costs a NOTE.
  if (depth <= 0L || !length(pkgs)) return(invisible(NULL))
  for (p in pkgs) {
    u <- archive_tarball(p)
    if (is.na(u)) { message(sprintf("  %s: no archive tarball found", p)); next }
    tf <- file.path(tempdir(), basename(u))
    if (!file.exists(tf)) try(download.file(u, tf, quiet = TRUE), silent = TRUE)
    if (!file.exists(tf)) next
    ## pre-install any dependency that is itself only in the Archive
    deps <- tryCatch({
      td <- tempfile(); dir.create(td)
      untar(tf, files = file.path(p, "DESCRIPTION"), exdir = td)
      d <- read.dcf(file.path(td, p, "DESCRIPTION"))[1, ]
      f <- intersect(c("Depends", "Imports", "LinkingTo"), names(d))
      x <- unlist(lapply(f, function(k) trimws(sub("\\(.*$", "", strsplit(d[[k]], ",")[[1]]))))
      x[nzchar(x) & x != "R"]
    }, error = function(e) character())
    ## Recurse into archived dependencies even when they are already present:
    ## an earlier pass may have installed one without ITS dependencies, and
    ## that only shows up as a load failure several levels up
    ## (LPKsample -> LPGraph -> PMA).
    need <- setdiff(deps, rownames(req$ap_all))
    if (length(need)) {
      message(sprintf("  %s: archived dependencies %s", p, paste(need, collapse = ", ")))
      install_archived(need, depth - 1L)
    }
    ## Install the repository dependencies by name rather than skipping the
    ## ones already present, so install.packages resolves THEIR dependencies.
    from_repo <- intersect(deps, rownames(req$ap_all))
    if (length(from_repo)) try(install.packages(from_repo, lib = LIB, dependencies = NA), silent = TRUE)
    message(sprintf("  %s: installing %s", p, basename(u)))
    try(install.packages(tf, lib = LIB, repos = NULL, type = "source"), silent = TRUE)
  }
}

for (pass in seq_len(MAXIT)) {
  a <- audit()
  ## Drop symlinks first, so install.packages sees them as missing and builds
  ## a real copy instead of deciding the version on the far end is good enough.
  ## Remove every symlink, not only the ones in the required set. The rest are
  ## equally capable of breaking when a shared dependency is updated, they are
  ## visible on the library path during checks, and they would not exist on
  ## another machine -- so leaving them makes the run unreproducible.
  all_linked <- symlinked()
  if (length(all_linked)) {
    message(sprintf("  removing %d symlinks into shared bundles (%d of them required)",
                    length(all_linked), length(intersect(all_linked, req$needed))))
    unlink(file.path(LIB, all_linked))
    a <- audit()
  }
  ## Only load-test once nothing is missing. Testing earlier flags every
  ## package whose dependency is merely not installed yet, which is true but
  ## useless: it would reinstall most of the library to learn nothing.
  ## Compare against what is still obtainable: the permanently unavailable
  ## packages would otherwise defer the load test on every pass forever, and
  ## anything it found would never be reinstalled.
  pending <- setdiff(a$missing, req$unavailable)
  broken <- if (length(pending)) {
    message(sprintf("  (load test deferred: %d still to install)", length(pending)))
    character()
  } else load_test(a$present)
  targets <- setdiff(unique(c(a$missing, a$outdated, broken)), req$unavailable)
  message(sprintf("\n== pass %d: missing %d, outdated %d, unloadable %d -> %d to install",
                  pass, length(a$missing), length(a$outdated), length(broken), length(targets)))
  if (!length(targets)) { message("library is complete, current and loadable"); break }

  from_repos <- setdiff(intersect(targets, c(rownames(req$ap_all), req$elsewhere)), req$archived)
  from_arch  <- intersect(targets, req$archived)
  if (length(from_repos)) install.packages(from_repos, lib = LIB, dependencies = NA)
  if (length(from_arch))  install_archived(from_arch)
}

a <- audit()
broken <- load_test(a$present)
report <- list(stamp = Sys.time(), needed = length(req$needed),
               present = length(a$present), missing = a$missing,
               outdated = a$outdated, unloadable = broken,
               symlinked = a$linked, unavailable = req$unavailable)
saveRDS(report, file.path(WORK, "library_report.rds"))
writeLines(a$missing,  file.path(WORK, "library_missing.txt"))
writeLines(broken,     file.path(WORK, "library_unloadable.txt"))

message(sprintf("\n=== library status ===\nrequired   : %d\npresent    : %d\nmissing    : %d\noutdated   : %d\nunloadable : %d\nunobtainable (expected): %d",
        length(req$needed), length(a$present), length(a$missing),
        length(a$outdated), length(broken), length(req$unavailable)))
if (length(a$missing)) message("MISSING: ",  paste(a$missing, collapse = ", "))
if (length(broken))    message("UNLOADABLE: ", paste(broken, collapse = ", "))
