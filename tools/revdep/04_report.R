#!/usr/bin/env Rscript
## Step 4 — compare the two sides and report every package that got worse.
##
## The comparison is PER CHECK ITEM, not on the overall status. Comparing only
## the overall status hides regressions: a package that already reports ERROR
## on both sides for an environmental reason (a vignette needing LaTeX, say)
## stays "ERROR" when the new side additionally breaks in its examples, and the
## real regression caused by pROC would never be reported. Item-level diffing
## catches that, and also makes the report say exactly what degraded.
##
## Severity ranks OK < NOTE < WARNING < ERROR < FAIL (check did not complete).

options(warn = 1)
WORK <- Sys.getenv("REVDEP_WORK", unset = file.path("/scratch", Sys.getenv("USER"), "pROC_revdeps"))
req  <- readRDS(file.path(WORK, "requirements.rds"))
run_dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else
  readLines(file.path(WORK, "last_run_dir.txt"))[1]
message("run directory: ", run_dir)

RANK <- c(OK = 0L, NOTE = 1L, WARNING = 2L, ERROR = 3L, FAIL = 4L)
rank_of <- function(x) { r <- RANK[x]; r[is.na(r)] <- 0L; r }

read_side <- function(side) {
  logs <- Sys.glob(file.path(run_dir, side, "*.Rcheck", "00check.log"))
  out <- lapply(logs, function(f) {
    pkg <- sub("\\.Rcheck$", "", sub("^rdepends_", "", basename(dirname(f))))
    txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
    it   <- grep("^\\* checking ", txt, value = TRUE)
    name <- trimws(sub("^\\* checking ", "", sub("\\.\\.\\..*$", "", it)))
    res  <- toupper(trimws(sub("^.*\\.\\.\\.", "", it)))
    res[!res %in% names(RANK)] <- "OK"     # e.g. "... Package", "... OK"
    ## An item can appear more than once; keep the worst result for each.
    items <- if (length(name)) tapply(rank_of(res), name, max) else integer()
    st <- grep("^Status:", txt, value = TRUE)
    overall <- if (!length(st)) "FAIL" else {
      s <- tail(st, 1)
      if (grepl("ERROR", s)) "ERROR" else if (grepl("WARNING", s)) "WARNING" else
      if (grepl("NOTE", s)) "NOTE" else "OK"
    }
    list(pkg = pkg, status = overall, items = items, log = f, text = txt)
  })
  setNames(out, vapply(out, `[[`, "", "pkg"))
}

old <- read_side("old"); new <- read_side("new")
## The package under test is checked on each side too, but it is not one of its
## own reverse dependencies: report it separately rather than in the table.
self_old <- old[[req$pkg]]; self_new <- new[[req$pkg]]
old[[req$pkg]] <- NULL; new[[req$pkg]] <- NULL
pkgs <- sort(union(names(old), names(new)))
message(sprintf("old side: %d packages, new side: %d", length(old), length(new)))
if (!length(pkgs)) stop("no check logs found under ", run_dir)

status_of <- function(side, p) if (is.null(side[[p]])) "MISSING" else side[[p]]$status
orank_of  <- function(side, p) if (is.null(side[[p]])) RANK[["FAIL"]] else RANK[[side[[p]]$status]]

## Items that degraded, per package.
degraded <- lapply(pkgs, function(p) {
  oi <- if (is.null(old[[p]])) integer() else old[[p]]$items
  ni <- if (is.null(new[[p]])) integer() else new[[p]]$items
  if (!length(ni)) return(character())
  common <- intersect(names(oi), names(ni))
  worse  <- common[ni[common] > oi[common]]
  ## an item present only on the new side, and not clean there
  newbad <- setdiff(names(ni)[ni > 0L], names(oi))
  c(worse, newbad)
})
names(degraded) <- pkgs

summary_df <- data.frame(
  package  = pkgs,
  old      = vapply(pkgs, function(p) status_of(old, p), ""),
  new      = vapply(pkgs, function(p) status_of(new, p), ""),
  old_rank = vapply(pkgs, function(p) orank_of(old, p), 0L),
  new_rank = vapply(pkgs, function(p) orank_of(new, p), 0L),
  n_degraded = vapply(degraded, length, 0L),
  stringsAsFactors = FALSE)
summary_df$verdict <- with(summary_df,
  ifelse(new_rank > old_rank | n_degraded > 0, "WORSE",
  ifelse(new_rank < old_rank, "better", "same")))
summary_df$degraded_items <- vapply(pkgs, function(p)
  paste(degraded[[p]], collapse = "; "), "")
write.csv(summary_df[, c("package","old","new","verdict","degraded_items")],
          file.path(run_dir, "revdep_summary.csv"), row.names = FALSE)

worse <- summary_df[summary_df$verdict == "WORSE", ]
worse <- worse[order(-worse$new_rank, -worse$n_degraded, worse$package), ]

sev <- names(RANK)
md <- c("# pROC reverse-dependency check: regressions", "",
  sprintf("- **%s**: CRAN **%s** (old) vs working tree **%s** (new)",
          req$pkg, req$cran_version, req$local_version),
  sprintf("- reverse dependencies checked: **%d**", length(req$tocheck)),
  sprintf("- run: `%s`", run_dir),
  sprintf("- generated: %s", format(Sys.time())),
  sprintf("- %s itself: old %s, new %s", req$pkg,
          if (is.null(self_old)) "not checked" else self_old$status,
          if (is.null(self_new)) "not checked" else self_new$status), "",
  "A package is listed when any individual check item got worse, or when its",
  "overall status got worse. Item-level comparison matters: a package already",
  "failing for an unrelated reason would otherwise mask a real regression.", "",
  sprintf("## Verdict: %d regression%s out of %d checked",
          nrow(worse), if (nrow(worse) == 1) "" else "s", length(pkgs)), "")

if (!nrow(worse)) {
  md <- c(md, "**No reverse dependency got worse.**", "")
} else {
  md <- c(md, "| package | old | new | what degraded |", "|---|---|---|---|",
          sprintf("| %s | %s | %s | %s |", worse$package, worse$old, worse$new,
                  ifelse(nzchar(worse$degraded_items), worse$degraded_items, "overall status")), "")
  for (p in worse$package) {
    md <- c(md, sprintf("### %s — %s to %s", p,
                        summary_df$old[summary_df$package == p],
                        summary_df$new[summary_df$package == p]), "")
    d <- degraded[[p]]
    if (length(d)) {
      oi <- if (is.null(old[[p]])) integer() else old[[p]]$items
      ni <- new[[p]]$items
      md <- c(md, sprintf("- `%s`: %s to %s", d,
                          ifelse(d %in% names(oi), sev[oi[d] + 1L], "absent"),
                          sev[ni[d] + 1L]), "")
    }
    if (!is.null(new[[p]])) {
      txt <- new[[p]]$text
      idx <- grep("^\\* checking .*(ERROR|WARNING)$", txt)
      if (length(idx)) {
        ex <- unlist(lapply(head(idx, 3), function(i) txt[i:min(i + 20L, length(txt))]))
        md <- c(md, "```", head(ex, 60), "```", "")
      }
      md <- c(md, sprintf("Full log: `%s`", new[[p]]$log), "")
    }
  }
}
counts <- table(factor(summary_df$verdict, levels = c("WORSE","same","better")))
md <- c(md, "## All packages", "",
        sprintf("- worse: %d", counts[["WORSE"]]),
        sprintf("- unchanged: %d", counts[["same"]]),
        sprintf("- improved: %d", counts[["better"]]), "",
        sprintf("Full table: `%s`", file.path(run_dir, "revdep_summary.csv")), "")
if (length(req$unavailable))
  md <- c(md, sprintf("Suggested but not installable on either side: %s",
                      paste(req$unavailable, collapse = ", ")), "")
md <- c(md, "",
  "Environment note: LaTeX comes from TinyTeX under $REVDEP_WORK, so vignettes",
  "build normally. PDF manual checks are still skipped (--no-manual): they test",
  "the formatting of each package's own Rd files and cannot be affected by this",
  "one. If a vignette fails on a missing .sty or font, run",
  "tools/revdep/install_missing_tex.R and check again.", "")

writeLines(md, file.path(run_dir, "revdep_report.md"))
message(sprintf("\n%d worse, %d same, %d better",
                counts[["WORSE"]], counts[["same"]], counts[["better"]]))
message("report: ", file.path(run_dir, "revdep_report.md"))
if (nrow(worse)) print(worse[, c("package","old","new","degraded_items")], row.names = FALSE)
