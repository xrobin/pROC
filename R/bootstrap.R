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
  # rename method into smooth.method for roc
  smoothing.args$roc1$smooth.method <- smoothing.args$roc1$method
  smoothing.args$roc1$method <- NULL
  smoothing.args$roc2$smooth.method <- smoothing.args$roc2$method
  smoothing.args$roc2$method <- NULL

  # Prepare arguments for later calls to roc
  auc1skeleton <- attributes(roc1$auc)
  auc1skeleton$roc <- NULL
  auc1skeleton$direction <- roc1$direction
  auc1skeleton$class <- NULL
  auc1skeleton$allow.invalid.partial.auc.correct <- TRUE
  auc1skeleton <- c(auc1skeleton, smoothing.args$roc1)
  names(auc1skeleton)[which(names(auc1skeleton) == "n")] <- "smooth.n"
  auc2skeleton <- attributes(roc2$auc)
  auc2skeleton$roc <- NULL
  auc2skeleton$direction <- roc2$direction
  auc2skeleton$class <- NULL
  auc2skeleton$allow.invalid.partial.auc.correct <- TRUE
  auc2skeleton <- c(auc2skeleton, smoothing.args$roc2)
  names(auc2skeleton)[which(names(auc2skeleton) == "n")] <- "smooth.n"

  auc1skeleton$auc <- auc2skeleton$auc <- TRUE

  # Some attributes may be duplicated in AUC skeletons and will mess the boostrap later on when we do.call().
  # If this condition happen, it probably means we have a bug elsewhere.
  # Rather than making a complicated processing to remove the duplicates,
  # just throw an error and let us solve the bug when a user reports it.
  duplicated.auc1skeleton <- duplicated(names(auc1skeleton))
  duplicated.auc2skeleton <- duplicated(names(auc2skeleton))
  if (any(duplicated.auc1skeleton)) {
    sessionInfo <- sessionInfo()
    save(roc1, roc2, boot.n, boot.stratified, boot.return, smoothing.args, sessionInfo, file = "pROC_bug.RData")
    stop(sprintf("pROC: duplicated argument(s) in AUC1 skeleton: \"%s\". Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", paste(names(auc1skeleton)[duplicated(names(auc1skeleton))], collapse = ", "), utils::packageDescription("pROC")$BugReports))
  }
  if (any(duplicated.auc2skeleton)) {
    sessionInfo <- sessionInfo()
    save(roc1, roc2, boot.n, boot.stratified, boot.return, smoothing.args, sessionInfo, file = "pROC_bug.RData")
    stop(sprintf("duplicated argument(s) in AUC2 skeleton: \"%s\". Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", paste(names(auc2skeleton)[duplicated(names(auc2skeleton))], collapse = ", "), utils::packageDescription("pROC")$BugReports))
  }
  if (!boot.stratified) {
    # The non-stratified resample rebuilds each curve from response/predictor,
    # so the skeletons need the levels and direction to do it without checks.
    auc1skeleton$levels <- roc1$levels
    auc1skeleton$direction <- roc1$direction
    auc2skeleton$levels <- roc2$levels
    auc2skeleton$direction <- roc2$direction
  }
  # One column per replicate, holding the statistic of each curve.
  resampled.values <- bootstrap.replicates(boot.n, bootstrap.test.replicate,
    roc1 = roc1, roc2 = roc2, stratified = boot.stratified, test = "boot",
    x = NULL, paired = TRUE,
    auc1skeleton = auc1skeleton, auc2skeleton = auc2skeleton,
    simplify = "columns", progress = progress, cl = cl
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
  # rename method into smooth.method for roc
  smoothing.args$roc1$smooth.method <- smoothing.args$roc1$method
  smoothing.args$roc1$method <- NULL
  smoothing.args$roc2$smooth.method <- smoothing.args$roc2$method
  smoothing.args$roc2$method <- NULL

  # Prepare arguments for later calls to roc
  auc1skeleton <- attributes(roc1$auc)
  auc1skeleton$roc <- NULL
  auc1skeleton$direction <- roc1$direction
  auc1skeleton$class <- NULL
  auc1skeleton$allow.invalid.partial.auc.correct <- TRUE
  auc1skeleton <- c(auc1skeleton, smoothing.args$roc1)
  names(auc1skeleton)[which(names(auc1skeleton) == "n")] <- "smooth.n"
  auc2skeleton <- attributes(roc2$auc)
  auc2skeleton$roc <- NULL
  auc2skeleton$direction <- roc2$direction
  auc2skeleton$class <- NULL
  auc2skeleton$allow.invalid.partial.auc.correct <- TRUE
  auc2skeleton <- c(auc2skeleton, smoothing.args$roc2)
  names(auc2skeleton)[which(names(auc2skeleton) == "n")] <- "smooth.n"

  auc1skeleton$auc <- auc2skeleton$auc <- test == "boot"

  # Some attributes may be duplicated in AUC skeletons and will mess the boostrap later on when we do.call().
  # If this condition happen, it probably means we have a bug elsewhere.
  # Rather than making a complicated processing to remove the duplicates,
  # just throw an error and let us solve the bug when a user reports it.
  duplicated.auc1skeleton <- duplicated(names(auc1skeleton))
  duplicated.auc2skeleton <- duplicated(names(auc2skeleton))
  if (any(duplicated.auc1skeleton)) {
    sessionInfo <- sessionInfo()
    save(roc1, roc2, test, x, paired, boot.n, boot.stratified, smoothing.args, sessionInfo, file = "pROC_bug.RData")
    stop(sprintf("pROC: duplicated argument(s) in AUC1 skeleton: \"%s\". Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", paste(names(auc1skeleton)[duplicated(names(auc1skeleton))], collapse = ", "), utils::packageDescription("pROC")$BugReports))
  }
  if (any(duplicated.auc2skeleton)) {
    sessionInfo <- sessionInfo()
    save(roc1, roc2, test, x, paired, boot.n, boot.stratified, smoothing.args, sessionInfo, file = "pROC_bug.RData")
    stop(sprintf("duplicated argument(s) in AUC2 skeleton: \"%s\". Diagnostic data saved in pROC_bug.RData. Please report this bug to <%s>.", paste(names(auc2skeleton)[duplicated(names(auc2skeleton))], collapse = ", "), utils::packageDescription("pROC")$BugReports))
  }

  if (!boot.stratified) {
    # The non-stratified resample rebuilds each curve from response/predictor,
    # so the skeletons need the levels and direction to do it without checks.
    auc1skeleton$levels <- roc1$levels
    auc1skeleton$direction <- roc1$direction
    auc2skeleton$levels <- roc2$levels
    auc2skeleton$direction <- roc2$direction
  }
  # One column per replicate, holding the statistic of each curve.
  resampled.values <- bootstrap.replicates(boot.n, bootstrap.test.replicate,
    roc1 = roc1, roc2 = roc2, stratified = boot.stratified, test = test,
    x = x, paired = paired,
    auc1skeleton = auc1skeleton, auc2skeleton = auc2skeleton,
    simplify = "columns", progress = progress, cl = cl
  )

  # compute the statistics
  diffs <- roc_utils_drop_na_replicates(
    resampled.values[1, ] - resampled.values[2, ]
  )

  # Restore smoothing if necessary
  if (smoothing.args$roc1$smooth) {
    smoothing.args$roc1$method <- smoothing.args$roc1$smooth.method
    roc1 <- do.call("smooth.roc", c(list(roc = roc1), smoothing.args$roc1))
  }
  if (smoothing.args$roc2$smooth) {
    smoothing.args$roc2$method <- smoothing.args$roc2$smooth.method
    roc2 <- do.call("smooth.roc", c(list(roc = roc2), smoothing.args$roc2))
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
  if (is.nan(D) && all(diffs == 0) && roc1$auc == roc2$auc) {
    D <- 0
  } # special case: no difference between AUCs produces a NaN

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
# 'simplify' says how the replicates are assembled, because the callers
# genuinely need different shapes:
#   "list"    one element per replicate, untouched
#   "vector"  one value per replicate
#   "columns" a matrix with one COLUMN per replicate (paired statistics)
#   "rows"    a matrix with one ROW per replicate (ci.se, ci.sp)
bootstrap.replicates <- function(boot.n, FUN, ...,
                                 simplify = c("list", "vector", "columns", "rows"),
                                 progress = FALSE, cl = NULL) {
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
  switch(simplify,
    list = replicates,
    vector = unlist(replicates),
    columns = simplify2array(replicates),
    rows = do.call(rbind, replicates)
  )
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
  # quote = TRUE matters. Some of these arguments are unevaluated calls -- the
  # smoothing call each replicate has to make -- and do.call() would otherwise
  # evaluate them here instead of passing them on.
  worker <- function(i) {
    assign(".Random.seed", rng.streams[[i]], envir = globalenv())
    do.call(fun, c(list(i), fun.args), quote = TRUE)
  }
  environment(worker) <- list2env(
    list(fun = FUN, fun.args = list(...), rng.streams = streams),
    parent = globalenv()
  )
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
roc_utils_rng_streams <- function(boot.n) {
  if (!exists(".Random.seed", envir = globalenv())) {
    set.seed(NULL)
  }
  seed <- sample.int(.Machine$integer.max, 1L)
  caller.seed <- get(".Random.seed", envir = globalenv())
  caller.kind <- RNGkind()
  on.exit({
    RNGkind(caller.kind[1], caller.kind[2], caller.kind[3])
    assign(".Random.seed", caller.seed, envir = globalenv())
  })
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
    if (cl == 1) {
      return(none)
    }
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

# One bootstrap resample of a ROC curve's observations.
#
# A stratified resample draws controls and cases separately, so each class
# keeps its original size; a non-stratified resample draws observations, so
# the class sizes vary between replicates. See the "stratified bootstrap"
# discussion in ?ci for which to prefer.
roc_utils_resample <- function(roc, stratified) {
  if (stratified) {
    controls <- roc$controls[sample.int(length(roc$controls), replace = TRUE)]
    cases <- roc$cases[sample.int(length(roc$cases), replace = TRUE)]
    predictor <- roc_utils_combine_predictor(controls, cases)
    response <- c(
      rep(roc$levels[1], length(controls)),
      rep(roc$levels[2], length(cases))
    )
  } else {
    idx <- sample.int(length(roc$predictor), replace = TRUE)
    predictor <- roc$predictor[idx]
    response <- roc$response[idx]
    splitted <- split(predictor, response)
    controls <- splitted[[as.character(roc$levels[1])]]
    cases <- splitted[[as.character(roc$levels[2])]]
  }
  list(controls = controls, cases = cases, predictor = predictor, response = response)
}

# One bootstrap replicate of a ROC curve: resample the observations and
# recompute the curve on them. Returns a 'roc' ready for auc(), coords() or
# smooth(), with every field consistent with the resampled data.
roc_utils_resampled_roc <- function(roc, stratified) {
  resampled <- roc_utils_resample(roc, stratified)
  thresholds <- roc_utils_thresholds(
    roc_utils_combine_predictor(resampled$controls, resampled$cases),
    roc$direction
  )
  perfs <- roc_utils_perfs_all(
    thresholds = thresholds,
    controls = resampled$controls,
    cases = resampled$cases,
    direction = roc$direction
  )
  scale <- ifelse(roc$percent, 100, 1)
  roc$controls <- resampled$controls
  roc$cases <- resampled$cases
  roc$predictor <- resampled$predictor
  roc$response <- resampled$response
  roc$thresholds <- thresholds
  roc$sensitivities <- perfs$se * scale
  roc$specificities <- perfs$sp * scale
  roc
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
roc_utils_drop_na_replicates <- function(x, margin = 1L) {
  if (is.matrix(x)) {
    bad <- apply(x, margin, anyNA)
  } else {
    bad <- is.na(x)
  }
  if (!any(bad)) {
    return(x)
  }
  warning(sprintf(
    "%i NA value(s) produced during bootstrap were ignored.", sum(bad)
  ))
  if (!is.matrix(x)) {
    x[!bad]
  } else if (margin == 1L) {
    x[!bad, , drop = FALSE]
  } else {
    x[, !bad, drop = FALSE]
  }
}

##########  AUC of two ROC curves (roc.test, cov)  ##########

# One bootstrap replicate of a paired or unpaired comparison of two curves.
# Returns the statistic for each curve, or c(NA, NA) if either resampled curve
# could not be built.
bootstrap.test.replicate <- function(n, roc1, roc2, stratified, test, x, paired,
                                     auc1skeleton, auc2skeleton) {
  if (stratified) {
    # Resample controls and cases separately, keeping each class size.
    idx.controls <- sample.int(length(roc1$controls), replace = TRUE)
    idx.cases <- sample.int(length(roc1$cases), replace = TRUE)
    auc1skeleton$controls <- roc1$controls[idx.controls]
    auc1skeleton$cases <- roc1$cases[idx.cases]
    if (paired) {
      # Paired curves must see the same patients.
      auc2skeleton$controls <- roc2$controls[idx.controls]
      auc2skeleton$cases <- roc2$cases[idx.cases]
    } else {
      auc2skeleton$controls <- roc2$controls[sample.int(length(roc2$controls), replace = TRUE)]
      auc2skeleton$cases <- roc2$cases[sample.int(length(roc2$cases), replace = TRUE)]
    }
    builder <- "roc_cc_nochecks"
  } else {
    idx <- sample.int(length(roc1$response), replace = TRUE)
    auc1skeleton$response <- roc1$response[idx]
    auc1skeleton$predictor <- roc1$predictor[idx]
    if (paired) {
      auc2skeleton$response <- roc2$response[idx]
      auc2skeleton$predictor <- roc2$predictor[idx]
    } else {
      idx2 <- sample.int(length(roc2$response), replace = TRUE)
      auc2skeleton$response <- roc2$response[idx2]
      auc2skeleton$predictor <- roc2$predictor[idx2]
    }
    builder <- "roc_rp_nochecks"
  }

  boot.roc1 <- try(do.call(builder, auc1skeleton), silent = TRUE)
  boot.roc2 <- try(do.call(builder, auc2skeleton), silent = TRUE)
  # A resampled curve may not be smoothable: report the replicate as missing.
  if (methods::is(boot.roc1, "try-error") || methods::is(boot.roc2, "try-error")) {
    return(c(NA_real_, NA_real_))
  }

  switch(test,
    sp = c(
      coords(boot.roc1, x = x, input = "specificity", ret = "sensitivity")[1, 1],
      coords(boot.roc2, x = x, input = "specificity", ret = "sensitivity")[1, 1]
    ),
    se = c(
      coords(boot.roc1, x = x, input = "sensitivity", ret = "specificity")[1, 1],
      coords(boot.roc2, x = x, input = "sensitivity", ret = "specificity")[1, 1]
    ),
    c(as.numeric(boot.roc1$auc), as.numeric(boot.roc2$auc))
  )
}

##########  AUC of one ROC curve (ci.auc, var)  ##########

ci_auc_bootstrap <- function(roc, conf.level, boot.n, boot.stratified, progress = FALSE,
                             cl = NULL, ...) {
  aucs <- roc_utils_drop_na_replicates(
    bootstrap.replicates(boot.n, bootstrap.auc,
      roc = roc, stratified = boot.stratified, simplify = "vector",
      progress = progress, cl = cl
    )
  )
  # TODO: Maybe apply a correction (it's in the Tibshirani?) What do Carpenter-Bithell say about that?
  quantile(aucs, c(0 + (1 - conf.level) / 2, .5, 1 - (1 - conf.level) / 2))
}

bootstrap.auc <- function(n, roc, stratified) {
  resampled <- roc_utils_resampled_roc(roc, stratified)
  # as.numeric() drops the 'roc' attribute auc.roc() attaches: it is the whole
  # resampled curve, which no caller of a bootstrap replicate ever reads.
  as.numeric(auc.roc(resampled,
    partial.auc = attr(roc$auc, "partial.auc"),
    partial.auc.focus = attr(roc$auc, "partial.auc.focus"),
    partial.auc.correct = attr(roc$auc, "partial.auc.correct"),
    allow.invalid.partial.auc.correct = TRUE
  ))
}

##########  AUC of a smooth ROC curve (ci.auc, var)  ##########

bootstrap.smooth.auc <- function(n, roc, stratified, smooth.roc.call, auc.call) {
  smooth.roc.call$roc <- roc_utils_resampled_roc(roc, stratified)
  auc.call$smooth.roc <- try(eval(smooth.roc.call), silent = TRUE)
  if (methods::is(auc.call$smooth.roc, "try-error")) {
    return(NA_real_)
  }
  as.numeric(eval(auc.call))
}

##########  SE and SP of a ROC curve (ci.se, ci.sp)  ##########

bootstrap.se <- function(n, roc, stratified, sp) {
  resampled <- roc_utils_resampled_roc(roc, stratified)
  coords.roc(resampled, sp, input = "specificity", ret = "sensitivity")[, 1]
}

bootstrap.sp <- function(n, roc, stratified, se) {
  resampled <- roc_utils_resampled_roc(roc, stratified)
  coords.roc(resampled, se, input = "sensitivity", ret = "specificity")[, 1]
}

##########  SE and SP of a smooth ROC curve (ci.se, ci.sp)  ##########

bootstrap.smooth.se <- function(n, roc, stratified, sp, smooth.roc.call) {
  smooth.roc.call$roc <- roc_utils_resampled_roc(roc, stratified)
  smooth.roc <- try(eval(smooth.roc.call), silent = TRUE)
  if (methods::is(smooth.roc, "try-error")) {
    return(NA_real_)
  }
  coords.smooth.roc(smooth.roc, sp, input = "specificity", ret = "sensitivity")[, 1]
}

bootstrap.smooth.sp <- function(n, roc, stratified, se, smooth.roc.call) {
  smooth.roc.call$roc <- roc_utils_resampled_roc(roc, stratified)
  smooth.roc <- try(eval(smooth.roc.call), silent = TRUE)
  if (methods::is(smooth.roc, "try-error")) {
    return(NA_real_)
  }
  coords.smooth.roc(smooth.roc, se, input = "sensitivity", ret = "specificity")[, 1]
}

##########  Thresholds of a ROC curve (ci.thresholds)  ##########

bootstrap.thresholds <- function(n, roc, stratified, thresholds) {
  resampled <- roc_utils_resample(roc, stratified)
  roc_utils_perfs_each(thresholds, resampled$controls, resampled$cases, roc$direction)
}

##########  Coords of a ROC curve (ci.coords)  ##########

bootstrap.coords <- function(n, roc, stratified, x, input, ret,
                             best.method, best.weights, best.policy) {
  resampled <- roc_utils_resampled_roc(roc, stratified)
  res <- coords.roc(resampled,
    x = x, input = input, ret = ret,
    best.method = best.method, best.weights = best.weights
  )
  enforce.best.policy.if.needed(res, x, best.policy)
}

##########  Coords of a smooth ROC curve (ci.coords)  ##########

bootstrap.smooth.coords <- function(n, roc, stratified, x, input, ret,
                                    best.method, best.weights, smooth.roc.call,
                                    best.policy) {
  smooth.roc.call$roc <- roc_utils_resampled_roc(roc, stratified)
  smooth.roc <- try(eval(smooth.roc.call), silent = TRUE)
  if (methods::is(smooth.roc, "try-error")) {
    return(NA)
  }
  # coords.smooth.roc(), not coords.roc(): a smoothed curve has no thresholds,
  # and the smooth method is what fills them with NA and resolves x = "best"
  # before delegating. Calling coords.roc() directly left the "best" search
  # comparing against absent thresholds.
  res <- coords.smooth.roc(smooth.roc,
    x = x, input = input, ret = ret,
    best.method = best.method, best.weights = best.weights
  )
  enforce.best.policy.if.needed(res, x, best.policy)
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
