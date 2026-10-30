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

##########  AUC of two ROC curves (roc.test, cov)  ##########

bootstrap.cov <- function(roc1, roc2, boot.n, boot.stratified, boot.return, smoothing.args,
                          progress = FALSE, cl = NULL) {
  # One column per replicate, holding the statistic of each curve.
  resampled.values <- bootstrap.replicates(boot.n, bootstrap.test.replicate,
    roc1 = roc1, roc2 = roc2, stratified = boot.stratified, test = "boot",
    x = NULL, paired = TRUE, smoothing.args = smoothing.args,
    FUN.VALUE = c(NA_real_, NA_real_), progress = progress, cl = cl
  )
  resampled.values <- roc_utils_drop_na_replicates(resampled.values, margin = 2L)

  cov <- stats::cov(resampled.values[1, ], resampled.values[2, ])
  if (boot.return) {
    attr(cov, "resampled.values") <- resampled.values
  }
  return(cov)
}

# Bootstrap test, used by roc.test.roc
bootstrap.test <- function(roc1, roc2, test, x, paired, boot.n, boot.stratified, smoothing.args,
                           progress = FALSE, cl = NULL) {
  # One column per replicate, holding the statistic of each curve.
  resampled.values <- bootstrap.replicates(boot.n, bootstrap.test.replicate,
    roc1 = roc1, roc2 = roc2, stratified = boot.stratified, test = test,
    x = x, paired = paired, smoothing.args = smoothing.args,
    FUN.VALUE = c(NA_real_, NA_real_), progress = progress, cl = cl
  )

  # compute the statistics
  diffs <- roc_utils_drop_na_replicates(
    resampled.values[1, ] - resampled.values[2, ]
  )

  # Restore smoothing if necessary
  if (!is.null(smoothing.args$roc1)) {
    roc1 <- roc_utils_resmooth(roc1, smoothing.args$roc1, reuse.auc = TRUE)
  }
  if (!is.null(smoothing.args$roc2)) {
    roc2 <- roc_utils_resmooth(roc2, smoothing.args$roc2, reuse.auc = TRUE)
  }

  if (test == "sp") {
    coord1 <- coords(roc1, x = x, input = c("specificity"), ret = c("sensitivity"))[1, 1]
    coord2 <- coords(roc2, x = x, input = c("specificity"), ret = c("sensitivity"))[1, 1]
    D <- (coord1 - coord2) / sd(diffs)
  } else if (test == "se") {
    coord1 <- coords(roc1, x = x, input = c("sensitivity"), ret = c("specificity"))[1, 1]
    coord2 <- coords(roc2, x = x, input = c("sensitivity"), ret = c("specificity"))[1, 1]
    D <- (coord1 - coord2) / sd(diffs)
  } else {
    D <- (roc1$auc - roc2$auc) / sd(diffs)
  }
  estimates.equal <- if (test %in% c("sp", "se")) coord1 == coord2 else roc1$auc == roc2$auc
  if (is.nan(D) && all(diffs == 0) && estimates.equal) {
    D <- 0
  } # special case: no difference between the tested values produces a NaN
  if (test %in% c("sp", "se")) {
    # the tested values, reported as estimates by roc.test
    attr(D, "estimate") <- c(coord1, coord2)
  }

  return(D)
}

##########  Running the replicates  ##########

# The one place bootstrap replicates are iterated.
#
# Calls FUN(i, ...) for i in 1:boot.n and assembles the results. Every
# bootstrap in pROC goes through here, so what applies to all of them --
# progress reporting and parallel execution -- has a single home instead of
# being repeated at two dozen call sites.
#
# FUN.VALUE is the NA replicate of the caller: a vector, or a matrix (a
# replicate that is a data frame, such as ci.coords', is taken as a matrix).
# Every replicate must have its length and type: vapply() makes a wrong shape
# an error instead of values shifted or recycled between replicates. The
# replicates are assembled as vapply() does, one column (or one slice of the
# array) per replicate, or with simplify = "rows" one row per replicate
# (ci.se, ci.sp).
bootstrap.replicates <- function(boot.n, FUN, ..., FUN.VALUE,
                                 simplify = c("columns", "rows"),
                                 progress = FALSE, cl = NULL) {
  boot.n <- roc_utils_check_boot_n(boot.n)
  simplify <- match.arg(simplify)
  force(FUN)
  cluster <- roc_utils_resolve_cluster(cl)
  if (!is.null(cluster$cluster)) {
    if (cluster$owned) {
      on.exit(stopCluster(cluster$cluster))
    }
    if (isTRUE(progress)) {
      message(sprintf(
        "Bootstrapping %i replicates on %i workers...",
        boot.n, length(cluster$cluster)
      ))
    }
    replicates <- bootstrap.replicates.parallel(cluster$cluster, boot.n, FUN, ...)
  } else if (isTRUE(progress)) {
    replicates <- bootstrap.replicates.with.progress(boot.n, FUN, ...)
  } else {
    replicates <- lapply(seq_len(boot.n), FUN, ...)
  }
  values <- vapply(replicates, function(value) {
    if (is.data.frame(value)) as.matrix(value) else value
  }, FUN.VALUE)
  if (simplify == "rows") {
    t(matrix(values, nrow = length(FUN.VALUE)))
  } else {
    values
  }
}

# boot.n is a replicate count, read as seq_len() does: a fraction is
# truncated, a vector uses its first element with a warning, and boot.n = 0
# runs no replicate. Values that seq_len() rejects (NA, negative, infinite or
# not a number) stop with a clear error. Values that seq_len() coerces (such
# as "10") are accepted. Returns the number of replicates.
roc_utils_check_boot_n <- function(boot.n) {
  n <- if (length(boot.n) >= 1 && is.atomic(boot.n)) as.numeric(boot.n[1]) else NA
  if (!is.finite(n) || n < 0) {
    stop("'boot.n' must be a non-negative number.")
  }
  if (length(boot.n) > 1) {
    warning("first element used of 'boot.n' argument")
  }
  trunc(n)
}

# As above, with a text progress bar.
#
# The bar is updated about a hundred times rather than once per replicate: a
# setTxtProgressBar() call costs around 6 microseconds against the 0.7-5
# milliseconds a replicate takes, which is a percent or so of the total if it
# is called every time, and nothing at all when throttled.
bootstrap.replicates.with.progress <- function(boot.n, FUN, ...) {
  every <- max(1L, boot.n %/% 100L)
  pb <- utils::txtProgressBar(min = 0, max = boot.n, style = 3)
  on.exit({
    utils::setTxtProgressBar(pb, boot.n)
    close(pb)
  })
  replicates <- vector("list", boot.n)
  for (i in seq_len(boot.n)) {
    replicates[[i]] <- FUN(i, ...)
    if (i %% every == 0L) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  replicates
}

# Run the replicates on a cluster.
#
# Each replicate gets its own L'Ecuyer-CMRG stream, so which worker happens to
# run it, and how many workers there are, make no difference to the result.
# The worker closure is built in an empty environment and given only what it
# needs: a closure defined here would otherwise drag this frame -- including
# the cluster object itself -- to every worker.
bootstrap.replicates.parallel <- function(cluster, boot.n, FUN, ...) {
  force(FUN)
  streams <- roc_utils_rng_streams(boot.n)
  # The extra arguments ride in the worker's environment rather than through
  # parLapply's '...': parLapply forwards those into clusterApply(x = , fun = )
  # and ci.coords() passes an argument called 'x', which would collide.
  #
  # quote = TRUE matters. The smoothing.args of a smoothed curve can hold
  # unevaluated expressions (such as a plotmath title), and do.call() would
  # otherwise evaluate them here instead of passing them on.
  #
  # The worker is made by a factory whose enclosure is the global environment,
  # so its own environment holds the factory's three arguments and nothing
  # else. They are forced so that no promise pointing back here is shipped.
  make.worker <- function(fun, fun.args, rng.streams) {
    force(fun)
    force(fun.args)
    force(rng.streams)
    function(i) {
      # The cluster may be the user's, its generator seeded for their own work
      # (clusterSetRNGStream()): leave it as it was, as on the master.
      had.seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      if (had.seed) {
        worker.seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      }
      worker.kind <- RNGkind()
      on.exit({
        # Set only the kinds the replicate changed: setting an unchanged
        # sample.kind = "Rounding" again would repeat its warning
        kinds <- list(kind = worker.kind[1], normal.kind = worker.kind[2], sample.kind = worker.kind[3])
        changed <- RNGkind() != worker.kind
        if (any(changed)) {
          do.call(RNGkind, kinds[changed])
        }
        if (had.seed) {
          assign(".Random.seed", worker.seed, envir = globalenv())
        } else {
          rm(".Random.seed", envir = globalenv())
        }
      })
      assign(".Random.seed", rng.streams[[i]], envir = globalenv())
      do.call(fun, c(list(i), fun.args), quote = TRUE)
    }
  }
  environment(make.worker) <- globalenv()
  worker <- make.worker(FUN, list(...), streams)
  parLapply(cluster, seq_len(boot.n), worker)
}

# One independent random number stream per replicate.
#
# L'Ecuyer-CMRG lets a single seed be split into sub-sequences far enough apart
# that they never overlap. Giving one to each replicate -- rather than one to
# each worker, as parallel::clusterSetRNGStream does -- means a replicate's
# numbers depend only on the seed and its index, so the same seed gives the
# same answer whatever the number of workers.
#
# The seed is drawn from the caller's own stream, so set.seed() still governs
# the result; the caller's generator is then restored exactly as it was.
# Restoring .Random.seed is enough: its first element encodes the three
# RNGkind() kinds, which R switches back to the next time it reads the seed
# (any draw, or RNGkind()). Unlike setting RNGkind() again, this does not
# repeat the warning about the "Rounding" sampler.
roc_utils_rng_streams <- function(boot.n) {
  if (!exists(".Random.seed", envir = globalenv())) {
    set.seed(NULL)
  }
  seed <- sample.int(.Machine$integer.max, 1L)
  caller.seed <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", caller.seed, envir = globalenv()))
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  streams <- vector("list", boot.n)
  stream <- get(".Random.seed", envir = globalenv())
  for (i in seq_len(boot.n)) {
    streams[[i]] <- stream
    stream <- nextRNGStream(stream)
  }
  streams
}

# What to run the bootstrap on.
#
# Returns the cluster and whether pROC created it, since a cluster pROC made
# must also be stopped by pROC. See ?ci for the accepted values.
roc_utils_resolve_cluster <- function(cl) {
  none <- list(cluster = NULL, owned = FALSE)
  if (is.null(cl) || isFALSE(cl)) {
    return(none)
  }
  if (isTRUE(cl)) {
    cluster <- getDefaultCluster()
    if (is.null(cluster)) {
      stop("'cl = TRUE' needs a default cluster: register one with parallel::setDefaultCluster(), pass a cluster object, or give the number of workers.")
    }
    return(list(cluster = cluster, owned = FALSE))
  }
  # Covers PSOCK, FORK, MPI and mirai clusters alike: they all extend
  # "cluster", so parLapply() works on any of them.
  if (inherits(cl, "cluster")) {
    return(list(cluster = cl, owned = FALSE))
  }
  if (is.numeric(cl) && length(cl) == 1L && !is.na(cl) && cl >= 1) {
    # Forking is nearly free to start, so creating the cluster for the
    # duration of one call costs little; Windows has no fork and pays the
    # socket cluster's startup instead.
    cluster <- if (.Platform$OS.type == "unix") {
      makeForkCluster(as.integer(cl))
    } else {
      makePSOCKcluster(as.integer(cl))
    }
    return(list(cluster = cluster, owned = TRUE))
  }
  stop("'cl' must be NULL or FALSE for a sequential bootstrap, a cluster from parallel::makeCluster(), TRUE to use parallel::getDefaultCluster(), or the number of workers.")
}

# Which progress bar, if any, the user asked for.
#
# TRUE/FALSE is the supported form. The pre-1.19 plyr-backed API took the name
# of a plyr progress bar ("none", "text", "win", "tk"), or the list stored in
# the pROCProgress option; those spellings are still understood so that old
# scripts and .Rprofile settings keep working, but only a text bar is drawn.
roc_utils_normalise_progress <- function(progress) {
  if (is.null(progress)) {
    return(FALSE)
  }
  if (is.list(progress)) {
    progress <- progress$name
  }
  if (is.character(progress)) {
    return(!identical(progress, "none"))
  }
  isTRUE(progress)
}

##########  Resampling  ##########

# The indices of one bootstrap resample of a ROC curve's observations.
#
# A stratified resample draws controls and cases separately, so each class
# keeps its original size; a non-stratified resample draws observations, so
# the class sizes vary between replicates. See the "stratified bootstrap"
# discussion in ?ci for which to prefer. Paired curves share these indices.
roc_utils_resample_idx <- function(roc, stratified) {
  if (stratified) {
    list(
      controls = sample.int(length(roc$controls), replace = TRUE),
      cases = sample.int(length(roc$cases), replace = TRUE)
    )
  } else {
    list(observations = sample.int(length(roc$predictor), replace = TRUE))
  }
}

# One bootstrap resample of a ROC curve's observations.
roc_utils_resample <- function(roc, stratified, idx = roc_utils_resample_idx(roc, stratified)) {
  if (stratified) {
    controls <- roc$controls[idx$controls]
    cases <- roc$cases[idx$cases]
    predictor <- roc_utils_combine_predictor(controls, cases)
    response <- c(
      rep(roc$levels[1], length(controls)),
      rep(roc$levels[2], length(cases))
    )
  } else {
    predictor <- roc$predictor[idx$observations]
    response <- roc$response[idx$observations]
    # Subset rather than split(): a class missing from the resample must give
    # an empty vector of the predictor's type (split() gives NULL), so that
    # ordered predictors can still be combined and the replicate dropped
    controls <- predictor[response %in% roc$levels[1]]
    cases <- predictor[response %in% roc$levels[2]]
  }
  list(controls = controls, cases = cases, predictor = predictor, response = response)
}

# The ROC curve of a resample (roc_utils_resample()), with every field
# consistent with the resampled data. Returns a 'roc' ready for auc(),
# coords() or smooth().
roc_utils_resampled_roc <- function(roc, resampled) {
  thresholds <- roc_utils_thresholds(
    roc_utils_combine_predictor(resampled$controls, resampled$cases),
    roc$direction
  )
  perfs <- roc_utils_perfs_all(
    thresholds = thresholds,
    controls = resampled$controls,
    cases = resampled$cases,
    direction = roc$direction,
    scale = ifelse(roc$percent, 100, 1)
  )
  roc$controls <- resampled$controls
  roc$cases <- resampled$cases
  roc$predictor <- resampled$predictor
  roc$response <- resampled$response
  roc$thresholds <- thresholds
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp
  roc
}

##########  One replicate  ##########

# One bootstrap replicate of one curve, as boot() calls its statistic: the
# only place where a replicate can fail. Resamples 'roc', smooths the
# resample again with 'smoothing.args' if given (the smoothing.args of a
# smoothed curve), and returns statistic(curve, ...). Returns 'na', NA
# values with the shape of a successful replicate given by the caller, when:
#  - a non-stratified resample lost a class (no control or no case);
#  - the resample could not be smoothed;
#  - statistic() returned NULL (no value defined on this resample, such as
#    coords "best" with no point in the partial AUC range).
# Errors in statistic() itself propagate. The caller then drops the NA
# replicates once, with one warning (roc_utils_drop_na_replicates()).
# With curve = FALSE, statistic() gets the resampled observations
# (roc_utils_resample()) instead of a curve, for ci.thresholds.
bootstrap.replicate <- function(n, roc, stratified, statistic, ..., na, smoothing.args = NULL, curve = TRUE) {
  resampled <- bootstrap.resampled.curve(roc, stratified, smoothing.args, curve = curve)
  if (is.null(resampled)) {
    return(na)
  }
  value <- statistic(resampled, ...)
  if (is.null(value)) na else value
}

# The resampled curve of a replicate, or NULL if it cannot be built.
# 'idx' lets paired curves share their resample (bootstrap.test.replicate).
bootstrap.resampled.curve <- function(roc, stratified, smoothing.args = NULL,
                                      idx = roc_utils_resample_idx(roc, stratified),
                                      curve = TRUE) {
  resampled <- roc_utils_resample(roc, stratified, idx)
  if (length(resampled$controls) == 0 || length(resampled$cases) == 0) {
    return(NULL)
  }
  if (!curve) {
    return(resampled)
  }
  resampled <- roc_utils_resampled_roc(roc, resampled)
  if (is.null(smoothing.args)) {
    return(resampled)
  }
  # Smoothers fail on some resamples (too few distinct values, fits that do
  # not converge): only the smoothing errors make the replicate NA
  tryCatch(roc_utils_resmooth(resampled, smoothing.args), error = function(e) NULL)
}

# Drop the replicates that produced NA and warn once.
#
# A resampled curve is occasionally degenerate -- a class can vanish from a
# non-stratified resample, and a smoother can fail to converge -- and the
# convention throughout pROC is to ignore those replicates rather than fail.
#
# 'x' is either a vector with one value per replicate, or a matrix. For a
# matrix, 'margin' says which dimension indexes the replicates: 1 for one row
# per replicate (ci.se, ci.sp), 2 for one column per replicate (the paired
# statistics of cov and roc.test). A replicate is always dropped whole, so
# that paired statistics stay aligned -- getting this margin wrong silently
# discarded a whole statistic instead of a replicate.
# 'ignore' flags, along the other dimension, the statistics that are NA in
# every replicate because they were requested for a missing value (NA 'x' in
# ci.se, ci.sp or ci.thresholds): they do not make a replicate NA.
roc_utils_drop_na_replicates <- function(x, margin = 1L, ignore = FALSE) {
  if (is.matrix(x)) {
    bad <- apply(x, margin, function(replicate) anyNA(replicate[!ignore]))
  } else {
    bad <- is.na(x)
  }
  if (!any(bad)) {
    return(x)
  }
  warning(sprintf(
    "%i NA value(s) produced during bootstrap were ignored.", sum(bad)
  ), call. = FALSE)
  if (!is.matrix(x)) {
    x[!bad]
  } else if (margin == 1L) {
    x[!bad, , drop = FALSE]
  } else {
    x[, !bad, drop = FALSE]
  }
}

##########  Statistics computed on each replicate  ##########

# AUC (ci.auc, var, roc.test, cov). 'auc.args' are the partial AUC settings
# of the original curve (bootstrap.auc.args()). auc.roc() also serves
# smoothed curves.
bootstrap.auc <- function(curve, auc.args) {
  # as.numeric() drops the 'roc' attribute auc.roc() attaches: it is the whole
  # resampled curve, which no caller of a bootstrap replicate ever reads.
  as.numeric(do.call(auc.roc, c(list(curve), auc.args, allow.invalid.partial.auc.correct = TRUE)))
}

# The partial AUC settings of an auc object, as bootstrap.auc() takes them
bootstrap.auc.args <- function(auc) {
  attributes(auc)[intersect(names(attributes(auc)), c("partial.auc", "partial.auc.focus", "partial.auc.correct"))]
}

# SE at given SP (ci.se) and SP at given SE (ci.sp)
bootstrap.se <- function(curve, sp) {
  coords(curve, sp, input = "specificity", ret = "sensitivity")[, 1]
}
bootstrap.sp <- function(curve, se) {
  coords(curve, se, input = "sensitivity", ret = "specificity")[, 1]
}

# SE and SP at given thresholds (ci.thresholds): on the resampled
# observations, without building the whole curve
bootstrap.thresholds <- function(resampled, thresholds, direction) {
  roc_utils_perfs_each(thresholds, resampled$controls, resampled$cases, direction)
}

# Coordinates (ci.coords). 'window' is the partial AUC range of a smoothed
# curve, where "best" is searched: the smoothing of the resample does not set
# it (roc_utils_resmooth() drops the AUC). An empirical resample keeps the auc
# of the original curve, and its range.
bootstrap.coords <- function(curve, x, input, ret, best.method, best.weights,
                             best.policy, window = NULL) {
  if (!is.null(window)) {
    # uncorrected: only the range is used
    curve$auc <- do.call(auc.roc, c(list(curve), window))
  }
  if (identical(x, "best") && !any(roc_utils_in_partial_auc(curve))) {
    # No point in the partial AUC range of this resample: NA replicate
    return(NULL)
  }
  if (is.numeric(x) && input %in% c("tn", "fp", "tp", "fn")) {
    # A count above the class size of this resample: NA, as coords() gives for an NA x
    empirical <- if (methods::is(curve, "smooth.roc")) attr(curve, "roc") else curve
    x[x > length(if (input %in% c("tn", "fp")) empirical$controls else empirical$cases)] <- NA
  }
  res <- coords(curve,
    x = x, input = input, ret = ret,
    best.method = best.method, best.weights = best.weights
  )
  ci_coords_numeric_columns(enforce.best.policy.if.needed(res, x, best.policy))
}

# The NA replicate of ci.coords
bootstrap.coords.na <- function(x, ret) {
  as.data.frame(matrix(NA_real_, length(x), length(ret), dimnames = list(NULL, ret)))
}

##########  Two curves (roc.test, cov)  ##########

# One bootstrap replicate of a paired or unpaired comparison of two curves.
# Returns the statistic for each curve, or c(NA, NA) if either resampled curve
# could not be built. 'smoothing.args' holds the smoothing.args of each
# smoothed curve (roc1, roc2), NULL for an empirical curve.
bootstrap.test.replicate <- function(n, roc1, roc2, stratified, test, x, paired, smoothing.args) {
  # Draw every index before building anything, so that the random stream
  # does not depend on which replicates fail
  idx1 <- roc_utils_resample_idx(roc1, stratified)
  # Paired curves must see the same patients.
  idx2 <- if (paired) idx1 else roc_utils_resample_idx(roc2, stratified)
  boot.roc1 <- bootstrap.resampled.curve(roc1, stratified, smoothing.args$roc1, idx1)
  boot.roc2 <- bootstrap.resampled.curve(roc2, stratified, smoothing.args$roc2, idx2)
  if (is.null(boot.roc1) || is.null(boot.roc2)) {
    return(c(NA_real_, NA_real_))
  }
  switch(test,
    sp = c(bootstrap.se(boot.roc1, x)[1], bootstrap.se(boot.roc2, x)[1]),
    se = c(bootstrap.sp(boot.roc1, x)[1], bootstrap.sp(boot.roc2, x)[1]),
    c(bootstrap.auc(boot.roc1, bootstrap.auc.args(roc1$auc)), bootstrap.auc(boot.roc2, bootstrap.auc.args(roc2$auc)))
  )
}

##########  AUC of one ROC curve (ci.auc)  ##########

ci_auc_bootstrap <- function(roc, conf.level, boot.n, boot.stratified, progress = FALSE,
                             cl = NULL, ...) {
  aucs <- roc_utils_drop_na_replicates(
    bootstrap.replicates(boot.n, bootstrap.replicate,
      roc = roc, stratified = boot.stratified,
      statistic = bootstrap.auc, auc.args = bootstrap.auc.args(roc$auc), na = NA_real_,
      FUN.VALUE = NA_real_, progress = progress, cl = cl
    )
  )
  # TODO: Maybe apply a correction (it's in the Tibshirani?) What do Carpenter-Bithell say about that?
  quantile(aucs, c(0 + (1 - conf.level) / 2, .5, 1 - (1 - conf.level) / 2))
}

# "best" can match several points; the policy decides which one the replicate
# contributes. Only reachable when x is the single value "best".
enforce.best.policy.if.needed <- function(res, x, best.policy) {
  if (length(x) == 1 && x == "best" && nrow(res) != 1) {
    enforce.best.policy(res, best.policy)
  } else {
    res
  }
}
