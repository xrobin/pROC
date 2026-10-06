#!/usr/bin/env Rscript
## Step 1 — what has to be checked, and what has to be installed to check it.
##
## Recomputed from live repository metadata every run: the reverse-dependency
## set drifts as packages are added to and archived from CRAN, so it must never
## be hard-coded or reused from a previous run.
##
## Writes into $REVDEP_WORK: tocheck.txt, needed.txt, needed_archived.txt,
## needed_elsewhere.txt, needed_unavailable.txt, requirements.rds

options(timeout = 900, warn = 1)
WORK <- Sys.getenv("REVDEP_WORK", unset = file.path("/scratch", Sys.getenv("USER"), "pROC_revdeps"))
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
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)

CRAN <- "https://cloud.r-project.org"
BIOC <- Sys.getenv("REVDEP_BIOC", unset = "3.23")
repos <- c(CRAN = CRAN,
           BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", BIOC),
           BioCann  = sprintf("https://bioconductor.org/packages/%s/data/annotation", BIOC),
           BioCexp  = sprintf("https://bioconductor.org/packages/%s/data/experiment", BIOC))

pkg  <- as.character(read.dcf(file.path(PKG_ROOT, "DESCRIPTION"), fields = "Package"))
lver <- as.character(read.dcf(file.path(PKG_ROOT, "DESCRIPTION"), fields = "Version"))

message("fetching repository metadata ...")
ap_cran <- available.packages(repos = CRAN,  type = "source")
ap_all  <- available.packages(repos = repos, type = "source")
stopifnot(pkg %in% rownames(ap_cran))
cran_ver <- ap_cran[pkg, "Version"]
message(sprintf("%s: CRAN %s vs local %s", pkg, cran_ver, lver))

## --- what to check -----------------------------------------------------
## CRAN's own practice: reverse Depends/Imports/LinkingTo/Suggests, not
## recursive, CRAN only (we can only check what we can fetch as source).
TOCHECK <- sort(unique(unlist(tools::package_dependencies(
  pkg, db = ap_cran, reverse = TRUE, which = "most"))))
message(sprintf("packages to check: %d", length(TOCHECK)))

## --- what to install ---------------------------------------------------
## R CMD check needs each package's Suggests as well as its strong deps,
## then the recursive strong closure of all of that.
lvl1 <- sort(unique(c(TOCHECK, unlist(tools::package_dependencies(
  TOCHECK, db = ap_all, which = c("Depends","Imports","LinkingTo","Suggests"),
  recursive = FALSE)))))
closure <- sort(unique(c(lvl1, unlist(tools::package_dependencies(
  lvl1, db = ap_all, which = "strong", recursive = TRUE)))))

desc <- read.dcf(file.path(PKG_ROOT, "DESCRIPTION"))[1, ]
own <- unlist(lapply(intersect(c("Depends","Imports","LinkingTo","Suggests"), names(desc)),
  function(k) { x <- trimws(sub("\\(.*$", "", strsplit(desc[[k]], ",")[[1]])); x[nzchar(x) & x != "R"] }))
own <- unique(c(own, unlist(tools::package_dependencies(own, db = ap_all,
                                                        which = "strong", recursive = TRUE))))

base_pkgs <- rownames(installed.packages(lib.loc = .Library, priority = "base"))
NEEDED <- sort(setdiff(unique(c(closure, own)), c(base_pkgs, pkg, "R")))
message(sprintf("packages to install: %d", length(NEEDED)))

## --- where each one can be obtained ------------------------------------
in_repos <- intersect(NEEDED, rownames(ap_all))
absent   <- setdiff(NEEDED, in_repos)

## Some absent ones were on CRAN and are in the Archive; some are declared
## by their dependant via Additional_repositories (that field is NOT in the
## PACKAGES index, so it must be read from the package database).
message("reading CRAN package database for Additional_repositories ...")
db <- tools::CRAN_package_db(); db <- db[!duplicated(db$Package), ]; rownames(db) <- db$Package
extra_repos <- unique(unlist(strsplit(
  db[intersect(TOCHECK, rownames(db)), "Additional_repositories"], "[, ]+")))
extra_repos <- extra_repos[!is.na(extra_repos) & nzchar(extra_repos)]
## Drop Bioconductor URLs declared here: the release is pinned above, and some
## packages still name an older one (3.22), which would otherwise be merged
## into options(repos) and could supply mismatched builds.
extra_repos <- grep("bioconductor\\.org", extra_repos, value = TRUE, invert = TRUE)

## Probe each declared repository separately, so one dead URL neither aborts
## the step nor buries the output in warnings.
elsewhere <- character(); live_repos <- character()
for (r in extra_repos) {
  ap_extra <- tryCatch(suppressWarnings(available.packages(repos = r, type = "source")),
                       error = function(e) NULL)
  if (is.null(ap_extra) || !nrow(ap_extra)) {
    message(sprintf("  additional repository unreachable, ignored: %s", r)); next
  }
  live_repos <- c(live_repos, r)
  elsewhere <- union(elsewhere, intersect(absent, rownames(ap_extra)))
}
extra_repos <- live_repos
if (length(extra_repos))
  message(sprintf("additional repositories in use: %s", paste(extra_repos, collapse = ", ")))

archived <- Filter(function(p) {
  u <- sprintf("https://cran.r-project.org/src/contrib/Archive/%s/", p)
  !inherits(try(suppressWarnings(readLines(u, n = 1, warn = FALSE)), silent = TRUE), "try-error")
}, setdiff(absent, elsewhere))

unavailable <- setdiff(absent, c(elsewhere, archived))

message(sprintf("  in CRAN/Bioc  : %d", length(in_repos)))
message(sprintf("  CRAN Archive  : %d  (%s)", length(archived), paste(archived, collapse = ", ")))
message(sprintf("  other repos   : %d  (%s)", length(elsewhere), paste(elsewhere, collapse = ", ")))
message(sprintf("  unobtainable  : %d  (%s)", length(unavailable), paste(unavailable, collapse = ", ")))

writeLines(TOCHECK,     file.path(WORK, "tocheck.txt"))
writeLines(NEEDED,      file.path(WORK, "needed.txt"))
writeLines(archived,    file.path(WORK, "needed_archived.txt"))
writeLines(elsewhere,   file.path(WORK, "needed_elsewhere.txt"))
writeLines(unavailable, file.path(WORK, "needed_unavailable.txt"))
saveRDS(list(pkg = pkg, local_version = lver, cran_version = cran_ver,
             repos = repos, extra_repos = extra_repos, bioc = BIOC,
             ap_cran = ap_cran, ap_all = ap_all,
             tocheck = TOCHECK, needed = NEEDED, archived = archived,
             elsewhere = elsewhere, unavailable = unavailable,
             stamp = Sys.time()),
        file.path(WORK, "requirements.rds"))
message(sprintf("\nwritten to %s", WORK))
