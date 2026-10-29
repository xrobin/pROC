# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison.
# Copyright (C) 2010-2014 Xavier Robin, Alexandre Hainard, Natacha Turck,
# Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez
# and Markus Müller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Helper functions for the ROC curves. These functions should not be called directly as they perform very specific tasks and do nearly no argument validity checks. Not documented in RD and not exported.

# Internal consistency check failed: save the data needed to reproduce it, and stop
roc_utils_bug <- function(msg, ...) {
  data <- list(..., sessionInfo = utils::sessionInfo())
  file <- tempfile("pROC_bug", fileext = ".RData")
  if (interactive()) {
    wd.file <- file.path(getwd(), basename(file))
    if (isTRUE(utils::askYesNo(sprintf("pROC hit an internal error. Save the data to reproduce it in %s, for a bug report?", wd.file)))) {
      file <- wd.file
    }
  }
  save(list = names(data), envir = list2env(data), file = file)
  stop(simpleError(sprintf(
    "pROC: %s. Data to reproduce it saved in %s. Please report this bug to <%s> (check the file for sensitive data before posting it).",
    msg, file, utils::packageDescription("pROC")$BugReports
  ), sys.call(-1)))
}

# Whether a roc (or smooth.roc) object was built from an ordered predictor.
# Checks every field that may carry the class so the answer is stable during
# construction (before original.predictor is assigned) and on finished objects.
roc_utils_is_ordered_roc <- function(roc) {
  is.ordered(roc$predictor) || is.ordered(roc$thresholds) ||
    is.ordered(roc$cases) || is.ordered(roc$controls) ||
    is.ordered(roc$original.predictor)
}

# A curve of an ordered predictor built by pROC 1.19 or earlier stores numeric
# thresholds (the level codes): its thresholds can't be read as levels
roc_utils_check_legacy_ordered <- function(roc) {
  if (is.ordered(roc$original.predictor) && is.numeric(roc$thresholds)) {
    stop(simpleError("This ROC curve of an ordered predictor was built by pROC 1.19 or earlier. Rebuild it with roc(), e.g. roc(x$original.response, x$original.predictor).", sys.call(-1)))
  }
}

# Concatenate ordered predictors. On R >= 4.1, c() of ordered inputs
# stays ordered only when the level sets are identical; mismatched
# levels become an unordered factor with the union of the labels.
# DESCRIPTION still supports R >= 2.14, where c() dropped to integer
# codes, so the helper remains required on both sides.
# levels(x) and levels(y) must be identical (same labels, same order).
roc_utils_c_ordered <- function(...) {
  args <- list(...)
  if (length(args) == 0) {
    stop("Need at least one ordered vector")
  }
  if (!all(vapply(args, is.ordered, logical(1)))) {
    stop("Arguments to roc_utils_c_ordered must be ordered.")
  }
  lvls <- levels(args[[1]])
  for (x in args) {
    if (!identical(levels(x), lvls)) {
      stop("Levels of cases and controls differ.")
    }
  }
  # exclude = NULL: an NA level is a level like the others, as in base R
  ordered(unlist(lapply(args, as.character), use.names = FALSE), levels = lvls, exclude = NULL)
}

# Concatenate predictors, preserving ordered class when needed.
roc_utils_combine_predictor <- function(...) {
  args <- list(...)
  if (any(vapply(args, is.ordered, logical(1)))) {
    return(do.call(roc_utils_c_ordered, args))
  }
  unlist(args, recursive = FALSE, use.names = FALSE)
}

# Auto direction from typical control vs case values.
# Ordered data are compared by as.integer() rank so median() is defined when n is even;
# cases/controls themselves are not converted to numeric.
roc_utils_auto_direction <- function(controls, cases) {
  if (is.ordered(controls)) {
    ctrl <- median(as.integer(controls))
    case <- median(as.integer(cases))
  } else {
    ctrl <- median(controls)
    case <- median(cases)
  }
  if (ctrl <= case) "<" else ">"
}

# P-scale (predictor/cases/controls) rank is as.integer(x).
# T-scale (thresholds) rank is as.integer(x) - (direction == ">"),
# because ">" prepends "-Inf" at position 1.
roc_utils_threshold_rank <- function(threshold, direction) {
  as.integer(threshold) - (direction == ">")
}

roc_utils_sesp_from_ranks <- function(thr_i, ctrl_i, case_i, direction, scale = 1) {
  if (direction == ">") {
    tp <- sum(case_i <= thr_i, na.rm = TRUE)
    tn <- sum(ctrl_i > thr_i, na.rm = TRUE)
  } else if (direction == "<") {
    tp <- sum(case_i >= thr_i, na.rm = TRUE)
    tn <- sum(ctrl_i < thr_i, na.rm = TRUE)
  }
  c(sp = scale * tn / length(ctrl_i), se = scale * tp / length(case_i))
}

roc_utils_c_thresholds <- function(...) {
  args <- list(...)
  if (any(vapply(args, is.ordered, logical(1)))) {
    lvls <- levels(args[[which(vapply(args, is.ordered, logical(1)))[1]]])
    ch <- unlist(lapply(args, as.character), use.names = FALSE)
    return(ordered(ch, levels = lvls))
  }
  unlist(args, recursive = FALSE, use.names = FALSE)
}

roc_utils_na_thresholds <- function(n, template) {
  if (is.ordered(template)) {
    ordered(rep(NA_character_, n), levels = levels(template))
  } else {
    rep(NA_real_, n)
  }
}

# The data.frames of coords(), built once with all their columns. data.frame()
# deparses every column and is most of the cost of a coords() call: list2DF()
# builds the same data.frame, unless a column (the threshold of a named or
# matrix 'x') gives row names or columns of its own, or must be recycled (the
# NA threshold of a smoothed curve).
roc_utils_coords_df <- function(columns, check.names = TRUE) {
  plain <- vapply(columns, function(column) is.null(names(column)) && is.null(dim(column)), logical(1))
  if (all(plain) && length(unique(lengths(columns))) == 1) {
    return(list2DF(columns))
  }
  do.call(data.frame, c(columns, check.names = check.names, stringsAsFactors = FALSE))
}

roc_utils_coords_basic <- function(threshold, specificity, sensitivity) {
  roc_utils_coords_df(list(
    threshold = threshold,
    specificity = as.vector(specificity),
    sensitivity = as.vector(sensitivity)
  ))
}

# returns a list of sensitivities (se) and specificities (sp) for the given data. Fast algorithm
# scale: 100 for percent curves. The rates are computed as scale * k / n, one
# rounding, so that a typed value such as 57.5 (23/40) matches the point it
# denotes, on percent as on fraction curves ((k / n) * 100 rounds twice).
roc_utils_perfs_all <- function(thresholds, controls, cases, direction, scale = 1) {
  if (is.ordered(thresholds) || is.ordered(cases) || is.ordered(controls)) {
    return(roc_utils_perfs_all_ordered(thresholds, controls, cases, direction, scale))
  }
  ncontrols <- length(controls)
  ncases <- length(cases)
  predictor <- c(controls, cases)
  response <- c(rep(0, length(controls)), rep(1, length(cases)))
  decr <- direction == "<"
  predictor.order <- order(predictor, decreasing = decr)
  predictor.sorted <- predictor[predictor.order]
  response.sorted <- response[predictor.order]

  tp <- cumsum(response.sorted == 1)
  fp <- cumsum(response.sorted == 0)
  se <- scale * tp / ncases
  sp <- scale * (ncontrols - fp) / ncontrols
  # filter duplicate thresholds
  dups.pred <- rev(duplicated(rev(predictor.sorted)))
  dups.sesp <- duplicated(se) & duplicated(sp)
  dups <- dups.pred | dups.sesp
  # Make sure we have the right length
  if (sum(!dups) != length(thresholds) - 1) {
    roc_utils_bug("fast algorithm computed an incorrect number of sensitivities and specificities", thresholds = thresholds, controls = controls, cases = cases, direction = direction)
  }
  if (direction == "<") {
    se <- rev(c(0, se[!dups]))
    sp <- rev(c(scale, sp[!dups]))
  } else {
    se <- c(0, se[!dups])
    sp <- c(scale, sp[!dups])
  }
  return(list(se = se, sp = sp))
}

roc_utils_perfs_all_ordered <- function(thresholds, controls, cases, direction, scale = 1) {
  user_levels <- if (is.ordered(cases)) levels(cases) else levels(controls)
  expected <- if (direction == "<") c(user_levels, "Inf") else c("-Inf", user_levels)
  stopifnot(is.ordered(thresholds), identical(levels(thresholds), expected))
  case_i <- as.integer(cases)
  ctrl_i <- as.integer(controls)
  offset <- as.integer(direction == ">")
  n <- length(thresholds)
  se <- numeric(n)
  sp <- numeric(n)
  for (i in seq_len(n)) {
    p <- roc_utils_sesp_from_ranks(as.integer(thresholds[i]) - offset, ctrl_i, case_i, direction, scale)
    sp[i] <- unname(p["sp"])
    se[i] <- unname(p["se"])
  }
  return(list(se = se, sp = sp))
}

roc_utils_fun_sesp <- function(...) {
  warning("pROC::roc$fun.sesp is deprecated")
  roc_utils_perfs_all(...)
}

# returns a vector with two elements, sensitivity and specificity, given the threshold at which to evaluate the performance, the values of controls and cases and the direction of the comparison, a character '>' or '<' as controls CMP cases
# sp <- roc_utils_perfs(...)[1,]
# se <- roc_utils_perfs(...)[2,]
roc_utils_perfs <- function(threshold, controls, cases, direction) {
  if (is.ordered(cases) || is.ordered(controls) || is.ordered(threshold)) {
    return(roc_utils_perfs_ordered(threshold, controls, cases, direction))
  }
  if (direction == ">") {
    tp <- sum(cases <= threshold)
    tn <- sum(controls > threshold)
  } else if (direction == "<") {
    tp <- sum(cases >= threshold)
    tn <- sum(controls < threshold)
  }
  # return(c(sp, se))
  return(c(sp = tn / length(controls), se = tp / length(cases)))
}

roc_utils_perfs_ordered <- function(threshold, controls, cases, direction) {
  user_levels <- if (is.ordered(cases)) levels(cases) else levels(controls)
  expected <- if (direction == "<") c(user_levels, "Inf") else c("-Inf", user_levels)
  if (!is.ordered(threshold) || !identical(levels(threshold), expected) || length(threshold) != 1L) {
    stop("Invalid ordered threshold.")
  }
  roc_utils_sesp_from_ranks(
    roc_utils_threshold_rank(threshold, direction),
    as.integer(controls),
    as.integer(cases),
    direction
  )
}

# Evaluate roc_utils_perfs at each threshold without dropping ordered class.
roc_utils_perfs_each <- function(thresholds, controls, cases, direction) {
  nthr <- length(thresholds)
  out <- matrix(NA_real_, nrow = 2L, ncol = nthr)
  rownames(out) <- c("sp", "se")
  if (is.ordered(cases) || is.ordered(controls) || is.ordered(thresholds)) {
    user_levels <- if (is.ordered(cases)) levels(cases) else levels(controls)
    expected <- if (direction == "<") c(user_levels, "Inf") else c("-Inf", user_levels)
    stopifnot(is.ordered(thresholds), identical(levels(thresholds), expected))
    case_i <- as.integer(cases)
    ctrl_i <- as.integer(controls)
    offset <- as.integer(direction == ">")
    for (i in which(!is.na(thresholds))) { # a missing threshold stays NA
      out[, i] <- roc_utils_sesp_from_ranks(as.integer(thresholds[i]) - offset, ctrl_i, case_i, direction)
    }
  } else {
    for (i in seq_len(nthr)) {
      out[, i] <- roc_utils_perfs(thresholds[i], controls = controls, cases = cases, direction = direction)
    }
  }
  out
}

# as roc_utils_perfs, but for densities
roc_utils_perfs_dens <- function(threshold, x, dens.controls, dens.cases, direction) {
  if (direction == ">") {
    tp <- sum(dens.cases[x <= threshold])
    tn <- sum(dens.controls[x > threshold])
  } else if (direction == "<") {
    tp <- sum(dens.cases[x >= threshold])
    tn <- sum(dens.controls[x < threshold])
  }
  # return(c(sp, se))
  return(c(sp = tn / sum(dens.controls), se = tp / sum(dens.cases)))
}

# return the thresholds to evaluate in the ROC curve, given the 'predictor' values. Returns all unique values of 'predictor' plus 2 extreme values
roc_utils_thresholds <- function(predictor, direction) {
  if (is.ordered(predictor)) {
    return(roc_utils_thresholds_ordered(predictor, direction))
  }
  unique.candidates <- sort(unique(predictor))
  thresholds1 <- (c(-Inf, unique.candidates) + c(unique.candidates, +Inf)) / 2
  thresholds2 <- (c(-Inf, unique.candidates) / 2 + c(unique.candidates, +Inf) / 2)
  thresholds <- ifelse(abs(thresholds1) > 1e100, thresholds2, thresholds1)
  if (any(ties <- thresholds %in% predictor)) {
    # If we get here, some thresholds are identical to the predictor
    # This is caused by near numeric ties that caused the mean to equal
    # one of the candidate
    # We need to make sure we select the right threshold more carefully
    if (direction == ">") {
      # We have:
      # tp <- sum(cases <= threshold)
      # tn <- sum(controls > threshold)
      # We need to make sure the selected threshold
      # Corresponds to the lowest observation of the predictor
      # Identify problematic thresholds
      # rows <- which(ties)
      for (tie.idx in which(ties)) {
        if (thresholds[tie.idx] == unique.candidates[tie.idx - 1]) {
          # We're already good, nothing to do
        } else if (thresholds[tie.idx] == unique.candidates[tie.idx]) {
          thresholds[tie.idx] <- unique.candidates[tie.idx - 1]
        } else {
          roc_utils_bug(sprintf("Couldn't fix near ties in thresholds: %s, %s, %s, %s", thresholds[tie.idx], unique.candidates[tie.idx - 1], unique.candidates[tie.idx], direction), predictor = predictor, direction = direction)
        }
      }
    } else if (direction == "<") {
      # We have:
      # tp <- sum(cases >= threshold)
      # tn <- sum(controls < threshold)
      # We need to make sure the selected threshold
      # Corresponds to the highest observation of the predictor
      # Identify the problematic thresholds:
      # rows <- which(apply(o, 1, any))
      for (tie.idx in which(ties)) {
        if (thresholds[tie.idx] == unique.candidates[tie.idx - 1]) {
          # Easy to fix: should be unique.candidates[tie.idx]
          thresholds[tie.idx] <- unique.candidates[tie.idx]
        } else if (thresholds[tie.idx] == unique.candidates[tie.idx]) {
          # We're already good, nothing to do
        } else {
          roc_utils_bug(sprintf("Couldn't fix near ties in thresholds: %s, %s, %s, %s", thresholds[tie.idx], unique.candidates[tie.idx - 1], unique.candidates[tie.idx], direction), predictor = predictor, direction = direction)
        }
      }
    }
  }
  return(thresholds)
}

# Ordered thresholds are the defined levels plus one Inf sentinel.
# direction "<": c(levels, Inf); direction ">": c(-Inf, levels).
# An observation equal to the threshold is classified as positive
# (>= when direction is "<", <= when direction is ">").
roc_utils_thresholds_ordered <- function(predictor, direction) {
  user_levels <- levels(predictor)
  if (direction == "<") {
    if ("Inf" %in% user_levels || "+Inf" %in% user_levels) {
      stop("Predictor level 'Inf' collides with the Inf sentinel used for direction \"<\". Rename the level.")
    }
    thr_levels <- c(user_levels, "Inf")
  } else {
    if ("-Inf" %in% user_levels) {
      stop("Predictor level '-Inf' collides with the -Inf sentinel used for direction \">\". Rename the level.")
    }
    thr_levels <- c("-Inf", user_levels)
  }
  ordered(thr_levels, levels = thr_levels, exclude = NULL)
}

# Find the local maximas of the ROC curve, as a logical vector: the
# thresholds that, for some weighting of sensitivity against specificity, are
# strictly better than each of their neighbours. Exact, on the counts: sorted
# by tn increasing and tp decreasing, a point is a local maximum iff the cross
# product of its incoming and outgoing directions is negative. (1, 0) before
# the first point and (0, -1) after the last compare each end with its single
# neighbour. Duplicated points count once, keeping 1.19.1's label.
roc_utils_max_thresholds_idx <- function(thresholds, sp, se, ncontrols, ncases, scale = 1) {
  tn <- round(sp / scale * ncontrols)
  tp <- round(se / scale * ncases)
  idx <- if (is.unsorted(sp)) rev(seq_along(sp)) else seq_along(sp)
  o <- idx[order(tn[idx], -tp[idx])]
  o <- o[!duplicated(cbind(tn[o], tp[o]))]
  local.maximas <- rep(FALSE, length(thresholds))
  n <- length(o)
  if (n == 1) {
    local.maximas[o] <- TRUE
    return(local.maximas)
  }
  dtn <- c(1, diff(tn[o]), 0)
  dtp <- c(0, diff(tp[o]), -1)
  local.maximas[o[dtn[1:n] * dtp[2:(n + 1)] - dtp[1:n] * dtn[2:(n + 1)] < 0]] <- TRUE
  return(local.maximas)
}

detect.env.true <- function(x) {
  value <- Sys.getenv(x)
  if (value == "T" || value == "True" || value == "TRUE" || value == "true") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# sort roc curve. Make sure specificities are increasing.
sort_roc <- function(roc) {
  if (is.unsorted(roc$specificities, na.rm = TRUE)) {
    roc$sensitivities <- rev(roc$sensitivities)
    roc$specificities <- rev(roc$specificities)
    roc$thresholds <- rev(roc$thresholds)
  }
  return(roc)
}

# sort smoothed roc curve. Make sure specificities are increasing.
sort_smooth_roc <- function(roc) {
  if (is.unsorted(roc$specificities, na.rm = TRUE)) {
    roc$sensitivities <- rev(roc$sensitivities)
    roc$specificities <- rev(roc$specificities)
  }
  return(roc)
}

# The list of valid coordinate arguments, without 'thresholds'
roc.utils.valid.coords <- c(
  "specificity", "sensitivity", "accuracy",
  "tn", "tp", "fn", "fp",
  "npv", "ppv", "fdr",
  "fpr", "tpr", "tnr", "fnr",
  "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv",
  "precision", "recall", "lr_pos", "lr_neg",
  "youden", "closest.topleft"
)

# Arguments which can be returned by coords
# @param threshold: FALSE for smooth.roc where threshold isn't valid
roc_utils_match_coords_ret_args <- function(x, threshold = TRUE) {
  valid.ret.args <- roc.utils.valid.coords
  if (threshold) {
    valid.ret.args <- c("threshold", valid.ret.args)
  }

  if ("all" %in% x) {
    if (length(x) > 1) {
      stop("ret='all' can't be used with other 'ret' options.")
    }
    x <- valid.ret.args
  }
  x <- replace(x, x == "topleft", "closest.topleft")
  x <- replace(x, x == "t", "threshold")
  x <- replace(x, x == "npe", "1-npv")
  x <- replace(x, x == "ppe", "1-ppv")
  # match.arg(several.ok = TRUE) silently drops the elements that don't match
  unknown <- x[is.na(pmatch(x, valid.ret.args, duplicates.ok = TRUE))]
  if (length(unknown) > 0) {
    quoted <- paste0("'", unknown, "'", collapse = ", ")
    if (length(unknown) == length(x)) {
      stop(simpleError(sprintf("Invalid coordinate(s) in 'ret': %s.", quoted), sys.call(-1)))
    }
    warning(simpleWarning(sprintf("Invalid coordinate(s) in 'ret' ignored: %s.", quoted), sys.call(-1)))
  }
  matched <- match.arg(x, valid.ret.args, several.ok = TRUE)
  return(matched)
}

# Arguments which can be used as input for coords
# @param threshold: FALSE for smooth.roc where threshold isn't valid
roc_utils_match_coords_input_args <- function(x, threshold = TRUE) {
  valid.args <- roc.utils.valid.coords
  if (threshold) {
    valid.args <- c("threshold", valid.args)
  }
  x <- replace(x, x == "topleft", "closest.topleft")
  x <- replace(x, x == "t", "threshold")
  x <- replace(x, x == "npe", "1-npv")
  x <- replace(x, x == "ppe", "1-ppv")
  matched <- match.arg(x, valid.args, several.ok = FALSE)
  # We only handle monotone coords
  if (!coord.is.monotone[matched]) {
    stop(sprintf("Coordinate '%s' is not monotone and cannot be used as input.", matched))
  }
  return(matched)
}

# Compute the min/max for partial AUC
# ... with an auc
roc_utils_min_partial_auc_auc <- function(auc) {
  roc_utils_min_partial_auc(attr(auc, "partial.auc"), attr(auc, "percent"))
}

roc_utils_max_partial_auc_auc <- function(roc) {
  roc_utils_max_partial_auc(attr(auc, "partial.auc"), attr(auc, "percent"))
}

# ... with partial.auc/percent
roc_utils_min_partial_auc <- function(partial.auc, percent) {
  if (!identical(partial.auc, FALSE)) {
    min <- sum(ifelse(percent, 100, 1) - partial.auc) * abs(diff(partial.auc)) / 2 / ifelse(percent, 100, 1)
  } else {
    min <- 0.5 * ifelse(percent, 100, 1)
  }
  return(min)
}

roc_utils_max_partial_auc <- function(partial.auc, percent) {
  if (!identical(partial.auc, FALSE)) {
    max <- abs(diff(partial.auc))
  } else {
    max <- 1 * ifelse(percent, 100, 1)
  }
  return(max)
}

# Checks if the
# Input: roc object
# Output: boolean, true the curve reaches 100%/100%, false otherwise
roc_utils_is_perfect_curve <- function(roc) {
  best.point <- max(roc$sensitivities + roc$specificities) / ifelse(roc$percent, 100, 1)
  return(abs(best.point - 2) < .Machine$double.eps^0.5) # or best.point == 2, with numerical tolerance
}

# roc() returns NA (NA values with na.rm = FALSE) or NaN (infinite predictor
# values) instead of a curve. The interfaces that build a curve with roc()
# check its result with this helper and return the NA/NaN unchanged.
# Input: the result of roc(); plot: whether the caller was asked to draw it
# Output: TRUE for a curve, FALSE (with a warning if plot = TRUE) otherwise
roc_utils_is_curve <- function(roc, plot = FALSE) {
  if (inherits(roc, c("roc", "smooth.roc"))) {
    return(TRUE)
  }
  if (plot) {
    warning(simpleWarning("Nothing plotted: roc() returned NA or NaN (NA values with na.rm = FALSE, or infinite predictor values).", call = sys.call(-1)))
  }
  FALSE
}

# Load package namespace 'pkg'.
# Input: package name
# Returns: TRUE upon success (or if the package was already loaded)
# Stops if the package can't be loaded
load.suggested.package <- function(pkg) {
  if (requireNamespace(pkg)) {
    return(TRUE)
  } else if (interactive()) {
    auto.install <- utils::askYesNo(sprintf("Package %s not available, do you want to install it now?", pkg))
    if (isTRUE(auto.install)) {
      utils::install.packages(pkg)
      if (requireNamespace(pkg)) {
        return(TRUE)
      } else {
        stop(sprintf("Installation of %s failed!", pkg))
      }
    }
  }
  stop(sprintf("Package '%s' not available.", pkg))
}


# Calculate coordinates
# @param roc: the roc curve, used to guess if data is in percent and number of cases and controls.
# @param thr, se, sp
# @param best.weights: see coords
# @return data.frame
# exact: TRUE for the rows that are points of an empirical curve (recycled).
# Their se and sp were built as scale * k / n (roc_utils_perfs_all), so
# round() recovers the counts k exactly: it inverts that construction, it is
# not a tolerance. Never set it for interpolated or smoothed values, whose
# counts are not integers. Every other coordinate is then computed from the
# counts, with one rounding.
roc_utils_calc_coords <- function(roc, thr, se, sp, best.weights, exact = FALSE) {
  ncases <- ifelse(methods::is(roc, "smooth.roc"), length(attr(roc, "roc")$cases), length(roc$cases))
  ncontrols <- ifelse(methods::is(roc, "smooth.roc"), length(attr(roc, "roc")$controls), length(roc$controls))
  substr.percent <- ifelse(roc$percent, 100, 1)

  exact <- rep_len(exact, length(se))
  tp <- se * ncases / substr.percent
  tp[exact] <- round(tp[exact])
  fn <- ncases - tp
  tn <- sp * ncontrols / substr.percent
  tn[exact] <- round(tn[exact])
  fp <- ncontrols - tn
  npv <- substr.percent * tn / (tn + fn)
  ppv <- substr.percent * tp / (tp + fp)
  # res <- matrix(NA, nrow = length(ret), ncol = length(se))
  # if ("tp" %in% ret) {}
  accuracy <- substr.percent * (tp + tn) / (ncases + ncontrols)
  precision <- ppv
  recall <- tpr <- se
  fpr <- substr.percent * fp / ncontrols
  tnr <- sp
  fnr <- substr.percent * fn / ncases
  fdr <- substr.percent * fp / (tp + fp)
  youden <- roc_utils_optim_crit(se, sp, substr.percent, best.weights, "youden")
  closest.topleft <- -roc_utils_optim_crit(se, sp, substr.percent, best.weights, "closest.topleft") / substr.percent

  return(roc_utils_coords_df(list(
    threshold = thr,
    sensitivity = se,
    specificity = sp,
    accuracy = accuracy,
    tn = tn,
    tp = tp,
    fn = fn,
    fp = fp,
    npv = npv,
    ppv = ppv,
    tpr = tpr,
    tnr = tnr,
    fpr = fpr,
    fnr = fnr,
    fdr = fdr,
    "1-specificity" = fpr,
    "1-sensitivity" = fnr,
    "1-accuracy" = substr.percent * (fp + fn) / (ncases + ncontrols),
    "1-npv" = substr.percent * fn / (tn + fn),
    "1-ppv" = fdr,
    precision = precision,
    recall = recall,
    "lr_pos" = tp * ncontrols / (fp * ncases),
    "lr_neg" = fn * ncontrols / (tn * ncases),
    youden = youden,
    closest.topleft = closest.topleft
  ), check.names = FALSE))
}

# Match arbitrary user-supplied thresholds to the threshold of the ROC curve.
# We need to be careful to assign x to the right thresholds around exact data point
# values. This means this function cannot look at the ROC thresholds themselves
# but must instead use the predictor values to assess the thresholds exactly.
# Returns the indices of the thresholds x, NA for a missing x.
# @param roc: the roc curve
# @param x: the threshold to determine indices
# @return integer vector of indices along roc$thresholds/roc$se/roc$sp.
roc_utils_thr_idx <- function(roc, x) {
  missing.x <- is.na(x)
  if (any(missing.x)) {
    # A missing threshold gives a NA row, as a missing value of any other input
    thr_idx <- rep(NA_integer_, length(x))
    thr_idx[!missing.x] <- roc_utils_thr_idx(roc, x[!missing.x])
    return(thr_idx)
  }
  if (roc_utils_is_ordered_roc(roc)) {
    return(roc_utils_thr_idx_ordered(roc, x))
  }
  cut_points <- sort(unique(roc$predictor))
  thr_idx <- rep(NA_integer_, length(x))
  if (roc$direction == "<") {
    cut_points <- c(cut_points, Inf)
    j <- 1
    o <- order(x)
    for (i in seq_along(x)) {
      t <- x[o[i]]
      while (cut_points[j] < t) {
        j <- j + 1
      }
      thr_idx[o[i]] <- j
    }
  } else {
    cut_points <- c(rev(cut_points), -Inf)
    j <- 1
    o <- order(x, decreasing = TRUE)
    for (i in seq_along(x)) {
      t <- x[o[i]]
      while (cut_points[j] > t) {
        j <- j + 1
      }
      thr_idx[o[i]] <- j
    }
  }
  return(thr_idx)
}

roc_utils_thr_idx_ordered <- function(roc, x) {
  ch <- as.character(x)
  x_pos <- match(ch, levels(roc$thresholds))
  if (anyNA(x_pos)) {
    unknown <- unique(ch[is.na(x_pos)])
    extra <- ""
    if (identical(unknown, "-Inf") && roc$direction == "<") {
      extra <- " Direction \"<\" uses the \"Inf\" sentinel."
    } else if (identical(unknown, "Inf") && roc$direction == ">") {
      extra <- " Direction \">\" uses the \"-Inf\" sentinel."
    }
    stop(sprintf("Unknown threshold level(s): %s.%s", paste(unknown, collapse = ", "), extra))
  }
  x_ranks <- x_pos - (roc$direction == ">")
  thr_ranks <- roc_utils_threshold_rank(roc$thresholds, roc$direction)
  matched <- match(x_ranks, thr_ranks)
  if (anyNA(matched)) {
    unknown <- unique(ch[is.na(matched)])
    stop(sprintf("Unknown threshold level(s): %s.", paste(unknown, collapse = ", ")))
  }
  return(matched)
}


# Whether x is a single non-missing number, to validate arguments.
# A logical counts as a number, as in base R arithmetic (and t.test(conf.level = TRUE))
roc_utils_is_single_number <- function(x) {
  (is.numeric(x) || is.logical(x)) && length(x) == 1 && !is.na(x)
}

# Check the best.weights argument
# @param weights: see coords(best.weights=)
# @param arg: the name of the argument, for the error messages
roc_utils_check_best_weights <- function(weights, arg = "best.weights") {
  if (!is.numeric(weights) || length(weights) != 2) {
    stop(sprintf("'%s' must be a numeric vector of length 2", arg))
  }
  # A NA or zero cost gives a NA criterion
  if (!is.na(weights[1]) && weights[1] < 0) {
    stop(sprintf("cost ('%s[1]') must not be negative.", arg))
  }
  if (is.na(weights[2]) || weights[2] <= 0 || weights[2] >= 1) {
    stop(sprintf("prevalence ('%s[2]') must be in the interval ]0,1[.", arg))
  }
}

# Get optimal criteria Youden or Closest Topleft
# @param se, sp: the roc curve
# @param max: the maximum value, 1 or 100, based on percent. Namely ifelse(percent, 100, 1)
# @param weights: see coords(best.weights=)
# @param method: youden or closest.topleft coords(best.method=)
# @return numeric vector along roc$thresholds/roc$se/roc$sp.
roc_utils_optim_crit <- function(se, sp, max, weights, method) {
  roc_utils_check_best_weights(weights)
  r <- (1 - weights[2]) / (weights[1] * weights[2]) # r should be 1 by default

  if (method == "youden") {
    optim.crit <- se + r * sp
  } else if (method == "closest.topleft" || method == "topleft") {
    optim.crit <- -((max - se)^2 + r * (max - sp)^2)
  }
  return(optim.crit)
}

# The points of a curve in the range of its partial AUC, where coords()
# searches "best": all of them without a partial AUC, or with ignore = TRUE
# (coords(ignore.partial.auc = TRUE))
roc_utils_in_partial_auc <- function(roc, ignore = FALSE) {
  partial.auc <- if (is.null(roc$auc)) FALSE else attr(roc$auc, "partial.auc")
  if (ignore || identical(partial.auc, FALSE)) {
    return(rep(TRUE, length(roc$sensitivities)))
  }
  values <- if (attr(roc$auc, "partial.auc.focus") == "sensitivity") roc$sensitivities else roc$specificities
  values <= partial.auc[1] & values >= partial.auc[2]
}

# The criterion of roc_utils_optim_crit computed from the counts of an
# empirical curve, to select the "best" thresholds with an exact ==:
# youden is tp * n0 + tn * n1 and closest.topleft is
# -(fn^2 * n0^2 + fp^2 * n1^2), roc_utils_optim_crit scaled by n0 * n1 and
# n0^2 * n1^2, so they order the thresholds the same way. They are exact
# integers only with the default weights (r = 1) and while they stay below
# 2^53: NULL otherwise (and for smoothed curves), to select on
# roc_utils_optim_crit as before.
roc_utils_optim_crit_counts <- function(roc, weights, method) {
  r <- (1 - weights[2]) / (weights[1] * weights[2])
  if (!isTRUE(r == 1) || methods::is(roc, "smooth.roc")) {
    return(NULL)
  }
  n1 <- length(roc$cases)
  n0 <- length(roc$controls)
  scale <- ifelse(roc$percent, 100, 1)
  # round() inverts scale * k / n (see roc_utils_calc_coords)
  tp <- round(roc$sensitivities * n1 / scale)
  tn <- round(roc$specificities * n0 / scale)
  if (method == "youden") {
    if (2 * n0 * n1 > 2^53) {
      return(NULL)
    }
    tp * n0 + tn * n1
  } else {
    if (2 * n0^2 * n1^2 > 2^53) {
      return(NULL)
    }
    -((n1 - tp)^2 * n0^2 + (n0 - tn)^2 * n1^2)
  }
}

coord.is.monotone <- c(
  "threshold" = TRUE,
  "sensitivity" = TRUE,
  "specificity" = TRUE,
  "accuracy" = FALSE,
  "tn" = TRUE,
  "tp" = TRUE,
  "fn" = TRUE,
  "fp" = TRUE,
  "npv" = FALSE,
  "ppv" = FALSE,
  "tpr" = TRUE,
  "tnr" = TRUE,
  "fpr" = TRUE,
  "fnr" = TRUE,
  "fdr" = FALSE,
  "1-specificity" = TRUE,
  "1-sensitivity" = TRUE,
  "1-accuracy" = FALSE,
  "1-npv" = FALSE,
  "1-ppv" = FALSE,
  "precision" = FALSE,
  "recall" = TRUE,
  "lr_pos" = FALSE,
  "lr_neg" = FALSE,
  "youden" = FALSE,
  "closest.topleft" = FALSE
)

coord.is.decreasing <- c(
  "threshold" = NA, # Depends on direction
  "sensitivity" = TRUE,
  "specificity" = FALSE,
  "accuracy" = NA,
  "tn" = FALSE,
  "tp" = TRUE,
  "fn" = FALSE,
  "fp" = TRUE,
  "npv" = NA,
  "ppv" = NA,
  "tpr" = TRUE,
  "tnr" = FALSE,
  "fpr" = TRUE,
  "fnr" = FALSE,
  "fdr" = NA,
  "1-specificity" = TRUE,
  "1-sensitivity" = FALSE,
  "1-accuracy" = NA,
  "1-npv" = NA,
  "1-ppv" = NA,
  "precision" = NA,
  "recall" = TRUE,
  "lr_pos" = NA,
  "lr_neg" = NA,
  "youden" = NA,
  "closest.topleft" = NA
)

# Get response and predictor(s) from a formula.
# This function takes care of all the logic to handle
# weights, subset, na.action etc. It handles formulas with
# and without data. It rejects weights and certain na.actions.
# @param formula
# @param data
# @param data.missing
# @param call the call from the parent function
# @param ... the ... from the parent function
# @return a list with 3 elements: response (vector), predictor.names (character),
#         predictors (data.frame).
roc_utils_extract_formula <- function(formula, data, data.missing, call, ...) {
  # Get predictors (easy)
  if (data.missing) {
    formula.terms <- terms(formula)
  } else {
    formula.terms <- terms(formula, data = data)
  }
  predictors <- attr(formula.terms, "term.labels")
  if (attr(formula.terms, "response") == 0) {
    stop("Error in the formula: a response is required in a formula of type response~predictor.")
  }

  indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(call), nomatch = 0)
  if (indx[1] == 0) {
    stop("A formula argument is required")
  }
  # Keep the standard arguments and run them in model.frame
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  # Only na.pass and na.fail should be used
  if (indx[5] != 0) {
    na.action.value <- as.character(call[indx[5]])
    if (!na.action.value %in% c("na.pass", "na.fail")) {
      warning(paste0(
        sprintf("Value %s of na.action is not supported ", na.action.value),
        "and will break pairing in roc.test and are.paired. ",
        "Please use 'na.rm = TRUE' instead."
      ))
    }
  } else {
    temp$na.action <- "na.pass"
  }

  # If data is missing we have to "find" it
  if (data.missing) {
    m <- tryCatch(
      {
        # First try normally. This works if the formula refers to variables
        # on the caller's environment or an attached database.
        eval.parent(temp, n = 2)
      },
      error = function(cond) {
        # If the function was called in a "with" call, the data is
        # a few frames up. 7 and 5 look like magic numbers here...
        # but they will change depending on the context and if the
        # roc function was wrapped. So we start with 7 and walk down
        # the frames until we find the data we want
        frames <- sys.frames()
        for (frame in frames[(length(frames) - 7):1]) {
          temp$data <- frame
          for (n in 5:9) {
            # n controls where the given formula will be found.
            # Here as well, 5-9 seem to be magic numbers
            try(
              {
                return(eval.parent(temp, n = n))
              },
              silent = TRUE
            )
          }
        }
        # If we're here it means we didn't find anything.
        # Stop with the initial message
        stop(cond)
      }
    )
  } else {
    # Run model.frame in the parent
    m <- eval.parent(temp, n = 2)
  }

  if (!is.null(model.weights(m))) stop("weights are not supported")

  # Sanity checks
  if (length(predictors) != ncol(m) - 1) {
    stop("Invalid formula: only formulas of type response~predictor or response~predictor1+predictor2+... are supported.")
  }
  stopifnot(all.equal(model.response(m), m[[1]], check.attributes = FALSE))

  return(list(
    response.name = names(m)[1],
    response = model.response(m),
    predictor.names = names(m)[-1],
    predictors = m[-1]
  ))
}

# Shared preconditions for the base-graphics add-on functions in
# R/plot.auc.R and R/plot.ci.R (polygon_auc, polygon_max_auc, polygon_ci,
# text.auc and friends, plus the pre-existing plot.ci.se/plot.ci.sp/
# plot.ci.thresholds), so they fail with a clear message instead of a
# low-level graphics error.

roc_utils_stop_if_no_device <- function(fun.name) {
  if (dev.cur() == 1) {
    stop(sprintf("'%s' needs an open plot. Call plot() first, or pass it to plot.roc(panel.first = ...).", fun.name))
  }
}

roc_utils_stop_if_no_auc <- function(x, arg = "x") {
  if (is.null(x$auc)) {
    stop(sprintf("'%s' has no 'auc'; call auc() on it (or roc(..., auc = TRUE)) first.", arg))
  }
}

# 'parallel' was deprecated in 1.19 when the bootstrap became sequential. The
# argument is still accepted so that old scripts keep running, but it has no
# effect: warn rather than silently ignoring a request to go parallel.
roc_utils_warn_deprecated_parallel <- function(parallel) {
  if (!identical(parallel, FALSE)) {
    warning("Parallel processing is deprecated in pROC 1.19. Ignoring 'parallel' argument")
  }
}

# Returns the empirical ROC curve a smooth.roc was built from, with the AUC
# specification of the smoothed curve (partial.auc etc.) rather than the one
# of the empirical curve, so that the smoothed curve can be rebuilt as it was.
roc_utils_unsmooth <- function(smooth.roc) {
  roc <- attr(smooth.roc, "roc")
  if (!is.null(smooth.roc$auc)) {
    # Only the specification matters: the (corrected) partial AUC of the
    # empirical curve may not be defined, and is never reported
    roc$auc <- auc(roc,
      partial.auc = attr(smooth.roc$auc, "partial.auc"),
      partial.auc.focus = attr(smooth.roc$auc, "partial.auc.focus"),
      partial.auc.correct = attr(smooth.roc$auc, "partial.auc.correct"),
      allow.invalid.partial.auc.correct = TRUE
    )
  }
  return(roc)
}

# Smooths 'roc' again as a curve smoothed with 'smoothing.args' (the
# smoothing.args of a smooth.roc). Every replay of a smoothing goes through
# here: roc.test, cov, are.paired and the bootstrap of smoothed curves.
# - The stored arguments are passed as values (quote = TRUE): a language
#   object such as roc(..., plot = TRUE, main = quote(italic(x))) stored is
#   not evaluated again.
# - The AUC specification is not a smoothing argument (smooth(r, partial.auc =
#   ...) and older versions of roc(smooth = TRUE) stored it): it is kept in
#   the auc of the curve.
# - The CI of 'roc' is dropped: smooth.roc(reuse.ci = TRUE) would bootstrap it
#   again. So is its AUC, unless reuse.auc is TRUE: the callers that only need
#   the curve (bootstrap replicates) should not compute it.
roc_utils_resmooth <- function(roc, smoothing.args, reuse.auc = FALSE) {
  roc$ci <- NULL
  if (!reuse.auc) {
    roc$auc <- NULL
  }
  smoothing.args[roc_utils_auc_arg_names] <- NULL
  do.call(smooth.roc, c(list(roc = roc), smoothing.args), quote = TRUE)
}

# plot.roc() for roc(..., plot = TRUE). The dots may carry ci.coords()' 'x'
# (roc(..., ci = TRUE, of = "coords", x = ...)), which would bind to the 'x'
# of plot.roc() and replace the curve: capture it here so it is not forwarded.
roc_utils_plot_roc <- function(curve, ..., x) {
  plot.roc(curve, ...)
}

# axis() for plot.roc() without explicit labels. The 'lab' graphical parameter
# (number of ticks, used by plot()) would partially match the 'labels' argument
# of axis(): capture it here so it is not forwarded.
roc_utils_axis <- function(..., lab) {
  axis(...)
}

# smooth.roc() for roc(..., smooth = TRUE). roc() forwards the same ... to
# smooth.roc, auc, ci and plot.roc. Capture here the few names of the other
# functions that the smoothers would take (as plot.default does for Axis):
# - x (ci.coords) would become the data of density() (BUG-076, BUG-133);
# - panel.first and panel.last (plot.roc) would be forced before the plot
#   exists when smooth.roc stores the dots (BUG-161);
# - partial.auc, partial.auc.focus, partial.auc.correct and
#   allow.invalid.partial.auc.correct (auc) would be passed to fitdistr and
#   the density function (BUG-042);
# - adj (graphical parameter) would partially match the adjust of density
#   smoothing (BUG-169).
# Everything else is forwarded as is.
roc_utils_smooth_roc <- function(roc, smooth.method, smooth.n, ...,
                                 x, panel.first, panel.last,
                                 partial.auc, partial.auc.focus, partial.auc.correct,
                                 allow.invalid.partial.auc.correct, adj) {
  smooth.roc(roc, method = smooth.method, n = smooth.n, ...)
}

# roc() for the default and formula methods of ci and ci.auc. Their ... also
# holds the method of ci.auc ("delong" or "bootstrap"), which roc(smooth = TRUE)
# would pass to smooth.roc as a second 'method' (BUG-166): capture it here.
roc_utils_roc_ci <- function(response, predictor, ..., method) {
  roc.default(response, predictor, ...)
}

# The AUC specification, which smooth.roc stores in smoothing.args when it is
# passed to it (see roc_utils_resmooth)
roc_utils_auc_arg_names <- c("partial.auc", "partial.auc.focus", "partial.auc.correct", "allow.invalid.partial.auc.correct")

# Curves smoothed from numeric densities cannot be resampled
roc_utils_stop_if_numeric_density <- function(smooth.roc, what) {
  args <- smooth.roc$smoothing.args
  if (is.null(args) || is.numeric(args$density.cases) || is.numeric(args$density.controls)) {
    stop(sprintf("Cannot compute %s ROC curves smoothed with numeric density.controls and density.cases.", what))
  }
}

# Set the call of a curve built by a wrapper (and the variable names, for the
# formula methods), also on the curve a smoothed curve keeps in attr(x, "roc"),
# from which print() takes the data names.
roc_utils_set_call_names <- function(roc, call, predictor.name = NULL, response.name = NULL) {
  roc$call <- call
  if (!is.null(predictor.name)) {
    roc$predictor.name <- predictor.name
    roc$response.name <- response.name
  }
  if (!is.null(attr(roc, "roc"))) {
    attr(roc, "roc") <- roc_utils_set_call_names(attr(roc, "roc"), call, predictor.name, response.name)
  }
  roc
}


# Labels of the thresholds marked by plot.roc(print.thres = ...) and
# geom_roc_threshold: threshold (the level for ordered predictors),
# specificity, sensitivity. With the default pattern (NULL), the threshold is
# shown as R prints a number, each one on its own (format() of the vector
# would pad them to common decimals); a custom pattern gets the number.
roc_utils_threshold_labels <- function(pattern, threshold, specificity, sensitivity, percent) {
  if (is.ordered(threshold) || is.character(threshold)) {
    threshold <- as.character(threshold)
  }
  if (is.null(pattern)) {
    pattern <- if (percent) "%s (%.1f%%, %.1f%%)" else "%s (%.3f, %.3f)"
    if (is.numeric(threshold)) {
      threshold <- vapply(threshold, format, "")
    }
  }
  sprintf(pattern, threshold, specificity, sensitivity)
}
