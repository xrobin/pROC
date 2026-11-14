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

roc.test <- function(...) {
  UseMethod("roc.test")
}

roc.test.formula <- function(formula, data, ..., subset, na.action) {
  data.missing <- missing(data)
  call <- match.call()
  roc.data <- roc_utils_extract_formula(formula, data, ...,
    data.missing = data.missing,
    call = call
  )
  if (length(roc.data$predictor.names) != 2) {
    stop("Invalid formula: exactly 2 predictors are required in a formula of type response~predictor1+predictor2.")
  }
  response <- roc.data$response
  predictors <- roc.data$predictors

  testres <- roc.test.default(response, predictors, ...)
  testres$call <- call
  # data.names for pretty print()ing
  if (data.missing) {
    testres$data.names <- sprintf("%s and %s by %s (%s, %s)", roc.data$predictor.names[1], roc.data$predictor.names[2], roc.data$response.name, testres$roc1$levels[1], testres$roc1$levels[2])
    testres$data.name <- testres$data.names # the standard htest element
  } else {
    testres$data.names <- sprintf("%s and %s in %s by %s (%s, %s)", roc.data$predictor.names[1], roc.data$predictor.names[2], deparse1(substitute(data)), roc.data$response.name, testres$roc1$levels[1], testres$roc1$levels[2])
    testres$data.name <- testres$data.names # the standard htest element
  }

  return(testres)
}

roc.test.default <- function(response, predictor1, predictor2 = NULL, na.rm = TRUE, method = NULL, ...) {
  if (is.matrix(predictor1) | is.data.frame(predictor1)) {
    if (!is.null(predictor2)) {
      stop("Predictor2 must not be specified if predictor1 is a matrix or a data.frame.")
    }
    if (dim(predictor1)[2] == 2 & length(response) == dim(predictor1)[1]) {
      roc1 <- roc(response, predictor1[, 1], ...)
      roc2 <- roc(response, predictor1[, 2], ...)
      if (!roc_utils_is_curve(roc1) || !roc_utils_is_curve(roc2)) {
        stop("Cannot compare ROC curves: roc() returned NA or NaN (NA values with na.rm = FALSE, or infinite predictor values).")
      }
      if (!is.null(names(predictor1))) {
        data.names <- sprintf("%s and %s in %s by %s (%s, %s)", names(predictor1)[1], names(predictor1)[2], deparse1(substitute(predictor1)), deparse1(substitute(response)), roc1$levels[1], roc1$levels[2])
      } else if (!is.null(colnames(predictor1))) {
        data.names <- sprintf("%s and %s in %s by %s (%s, %s)", colnames(predictor1)[1], colnames(predictor1)[2], deparse1(substitute(predictor1)), deparse1(substitute(response)), roc1$levels[1], roc1$levels[2])
      } else {
        data.names <- sprintf("%s by %s (%s, %s)", deparse1(substitute(predictor1)), deparse1(substitute(response)), roc1$levels[1], roc1$levels[2])
      }
    } else {
      stop("Wrong dimension for predictor1 as a matrix or a data.frame.")
    }
  } else {
    if (is.null(predictor2)) {
      stop("Missing argument predictor2 with predictor1 as a vector.")
    }
    # Need to remove NAs
    if (na.rm) {
      nas <- is.na(response) | is.na(predictor1) | is.na(predictor2)
      response <- response[!nas]
      predictor1 <- predictor1[!nas]
      predictor2 <- predictor2[!nas]
    }
    roc1 <- roc(response, predictor1, ...)
    roc2 <- roc(response, predictor2, ...)
    if (!roc_utils_is_curve(roc1) || !roc_utils_is_curve(roc2)) {
      stop("Cannot compare ROC curves: roc() returned NA or NaN (NA values with na.rm = FALSE, or infinite predictor values).")
    }
    call <- match.call()
    data.names <- sprintf("%s and %s by %s (%s, %s)", deparse1(call$predictor1), deparse1(call$predictor2), deparse1(call$response), roc1$levels[1], roc1$levels[2])
  }
  test <- roc.test.roc(roc1, roc2, method = method, ...)
  test$data.names <- data.names
  test$data.name <- test$data.names # the standard htest element
  return(test)
}

roc.test.auc <- function(roc1, roc2, ...) {
  # First save the names
  data.names <- paste(deparse1(substitute(roc1)), "and", deparse1(substitute(roc2)))
  # Change roc1 from an auc to a roc object but keep the auc specifications
  auc1 <- roc1
  attr(auc1, "roc") <- NULL
  roc1 <- attr(roc1, "roc")
  roc1$auc <- auc1
  # Pass to roc.test.roc
  testres <- roc.test.roc(roc1, roc2, ...)
  testres$call <- match.call()
  testres$data.names <- data.names
  testres$data.name <- testres$data.names # the standard htest element
  return(testres)
}

roc.test.smooth.roc <- function(roc1, roc2, ...) {
  testres <- roc.test.roc(roc1, roc2, ...)
  testres$call <- match.call()
  testres$data.names <- paste(deparse1(substitute(roc1)), "and", deparse1(substitute(roc2)))
  testres$data.name <- testres$data.names # the standard htest element
  return(testres)
}

# DeLong's test for k paired ROC curves: DeLong et al. (1988) p. 841
roc.test.list <- function(rocs, contrast = NULL,
                          alternative = c("two.sided", "less", "greater"),
                          conf.level = 0.95, ...) {
  chkDots(...)
  alternative <- match.arg(alternative)
  data.name <- deparse1(substitute(rocs))
  rocs <- delong_list_rocs(rocs)
  k <- length(rocs)

  if (is.null(contrast)) {
    # the k - 1 successive differences: all the AUCs are equal
    contrast <- cbind(diag(k - 1), 0) - cbind(0, diag(k - 1))
  } else if (!is.matrix(contrast)) {
    contrast <- matrix(contrast, nrow = 1)
  }
  if (!is.numeric(contrast) || ncol(contrast) != k || anyNA(contrast)) {
    stop("'contrast' must be a numeric vector with one coefficient per ROC curve, or a matrix with one such row per contrast.")
  }

  estimate <- vapply(rocs, function(roc) as.numeric(roc$auc), numeric(1))
  names(estimate) <- paste("AUC of", if (is.null(names(rocs))) paste0("roc", seq_len(k)) else names(rocs))
  htest <- list(
    alternative = alternative,
    data.name = data.name,
    estimate = estimate,
    method = sprintf("DeLong's test for %d correlated ROC curves", k)
  )
  class(htest) <- "htest"

  calcs <- delong.cov(rocs)
  d <- drop(contrast %*% calcs$theta) # L theta
  M <- contrast %*% calcs$S %*% t(contrast) # L S L'

  if (nrow(contrast) == 1) {
    # One contrast: L theta / sqrt(L S L') is standard normal
    roc_test_check_conf_level(conf.level)
    sig <- sqrt(M[1, 1])
    stat <- delong.test(list(d = d, sig = sig))
    htest$statistic <- c("Z" = stat)
    htest$null.value <- c("contrast of AUCs" = 0)
    htest$conf.int <- roc_test_ci_normal(d, sig, conf.level, alternative)
    if (rocs[[1]]$percent) {
      # the placements are on the 0-1 scale, the estimates on 0-100
      htest$conf.int <- htest$conf.int * 100
    }
    htest$p.value <- switch(alternative,
      two.sided = 2 * pnorm(-abs(stat)),
      greater = pnorm(-stat),
      less = pnorm(stat)
    )
  } else {
    # Several contrasts, Eq. 5: (L theta)' (L S L')^-1 (L theta) is chi-square
    # with the rank of L S L' as degrees of freedom
    if (alternative != "two.sided") {
      stop("Only two-sided tests are available with several contrasts.")
    }
    if (anyNA(M)) {
      # A single case or control: the variance is undefined
      stat <- df <- NA_real_
    } else {
      # Redundant contrasts do not change the statistic: keep a full-rank subset
      qrM <- qr(M)
      df <- qrM$rank
      keep <- qrM$pivot[seq_len(df)]
      # rank 0: no difference between the curves, as for delong.test
      stat <- if (df > 0) sum(d[keep] * solve(M[keep, keep, drop = FALSE], d[keep])) else 0
    }
    htest$statistic <- c("X-squared" = stat)
    htest$parameter <- c("df" = df)
    htest$p.value <- pchisq(stat, df, lower.tail = FALSE)
  }
  htest$contrast <- contrast
  htest$rocs <- rocs
  return(htest)
}

roc.test.roc <- function(roc1, roc2,
                         method = c("delong", "bootstrap", "venkatraman", "sensitivity", "specificity"),
                         sensitivity = NULL, specificity = NULL,
                         alternative = c("two.sided", "less", "greater"),
                         paired = NULL,
                         reuse.auc = TRUE,
                         boot.n = 2000, boot.stratified = TRUE,
                         ties.method = "random",
                         progress = getOption("pROCProgress", interactive()),
                         parallel = FALSE,
                         conf.level = 0.95,
                         cl = NULL,
                         ...) {
  alternative <- match.arg(alternative)
  data.names <- paste(deparse1(substitute(roc1)), "and", deparse1(substitute(roc2)))
  # If roc2 is an auc, take the roc but keep the auc specifications
  if (methods::is(roc2, "auc")) {
    auc2 <- roc2
    attr(auc2, "roc") <- NULL
    roc2 <- attr(roc2, "roc")
    roc2$auc <- auc2
  }

  if (roc_utils_is_perfect_curve(roc1) && roc_utils_is_perfect_curve(roc2)) {
    warning("roc.test() of two ROC curves with AUC == 1 has always p.value = 1 and can be misleading.")
  }
  progress <- roc_utils_normalise_progress(progress)
  roc_utils_warn_deprecated_parallel(parallel)

  # store which objects are smoothed, and how
  smoothing.args <- list()
  smooth1 <- methods::is(roc1, "smooth.roc")
  if (smooth1) {
    roc_utils_stop_if_numeric_density(roc1, "the statistic on")
    smoothing.args$roc1 <- roc1$smoothing.args
    smooth.call1 <- roc1$call
    roc1 <- roc_utils_unsmooth(roc1)
  }
  smooth2 <- methods::is(roc2, "smooth.roc")
  if (smooth2) {
    roc_utils_stop_if_numeric_density(roc2, "the statistic on")
    smoothing.args$roc2 <- roc2$smoothing.args
    smooth.call2 <- roc2$call
    roc2 <- roc_utils_unsmooth(roc2)
  }

  # Check if we do a paired or unpaired roc.test
  rocs.are.paired <- roc_utils_paired(roc1, roc2, paired)
  paired <- isTRUE(rocs.are.paired)
  if (paired) {
    roc1 <- attr(rocs.are.paired, "roc1")
    roc2 <- attr(rocs.are.paired, "roc2")
  }

  # check that the AUC was computed, or do it now
  if (is.null(roc1$auc) | !reuse.auc) {
    if (smooth1) {
      roc1$auc <- auc(smooth.roc = roc_utils_resmooth(roc1, smoothing.args$roc1), ...)
    } else {
      roc1$auc <- auc(roc1, ...)
    }
  }
  if (is.null(roc2$auc) | !reuse.auc) {
    if (smooth2) {
      roc2$auc <- auc(smooth.roc = roc_utils_resmooth(roc2, smoothing.args$roc2), ...)
    } else {
      roc2$auc <- auc(roc2, ...)
    }
  }

  # check that the same region was requested in auc. Otherwise, issue a warning
  if (!identical(attributes(roc1$auc)[names(attributes(roc1$auc)) != "roc"], attributes(roc2$auc)[names(attributes(roc2$auc)) != "roc"])) {
    warning("Different AUC specifications in the ROC curves. Enforcing the inconsistency, but unexpected results may be produced.")
  }
  # check that the same smoothing params were requested in auc. Otherwise, issue a warning
  if (!identical(smoothing.args$roc1, smoothing.args$roc2)) {
    warning("Different smoothing parameters in the ROC curves. Enforcing the inconsistency, but unexpected results may be produced.")
  }

  # Check the method
  if (missing(method) | is.null(method)) {
    # determine method if missing
    if (has.partial.auc(roc1) || has.partial.auc(roc2)) {
      # partial auc: go for bootstrap
      method <- "bootstrap"
    } else if (smooth1 || smooth2) {
      # smoothing in one or both: bootstrap
      method <- "bootstrap"
    } else if (roc1$direction != roc2$direction) {
      # delong doesn't work well with opposite directions (will report high significance if roc1$auc and roc2$auc are similar and high)
      method <- "bootstrap"
    } else {
      method <- "delong"
    }
  } else {
    method <- match.arg(method)
    if (method == "delong") {
      # delong NA to pAUC: warn + change
      if (has.partial.auc(roc1) || has.partial.auc(roc2)) {
        stop("DeLong's test is not supported for partial AUC. Use method=\"bootstrap\" instead.")
      }
      if (smooth1 || smooth2) {
        stop("DeLong's test is not supported for smoothed ROCs. Use method=\"bootstrap\" instead.")
      }
      if (roc1$direction != roc2$direction) {
        warning("DeLong's test should not be applied to ROC curves with a different direction.")
      }

    } else if (method == "venkatraman") {
      if (has.partial.auc(roc1) || has.partial.auc(roc2)) {
        stop("Partial AUC is not supported for Venkatraman's test.")
      }
      if (smooth1 || smooth2) {
        stop("Venkatraman's test is not supported for smoothed ROCs")
      }
      if (roc1$direction != roc2$direction) {
        warning("Venkatraman's test should not be applied to ROC curves with different directions.")
      }
      if (alternative != "two.sided") {
        stop("Only two-sided tests are available for Venkatraman.")
      }
    }
  }

  # Prepare the return value htest
  if (smooth1) {
    estimate <- roc_utils_resmooth(roc1, smoothing.args$roc1, reuse.auc = TRUE)$auc
  } else {
    estimate <- roc1$auc
  }
  if (smooth2) {
    estimate <- c(estimate, roc_utils_resmooth(roc2, smoothing.args$roc2, reuse.auc = TRUE)$auc)
  } else {
    estimate <- c(estimate, roc2$auc)
  }
  if (identical(attr(roc1$auc, "partial.auc"), FALSE)) {
    nest <- paste(ifelse(smooth1, "Smoothed ", ""), "AUC of roc1", sep = "")
  } else {
    nest <- paste(ifelse(attr(roc1$auc, "partial.auc.correct"), "Corrected ", ""),
      ifelse(smooth1, "Smoothed ", ""),
      "pAUC (", attr(roc1$auc, "partial.auc")[1], "-", attr(roc1$auc, "partial.auc")[2], " ", attr(roc1$auc, "partial.auc.focus"),
      ") of roc1",
      sep = ""
    )
  }
  if (identical(attr(roc2$auc, "partial.auc"), FALSE)) {
    nest <- c(nest, paste(ifelse(smooth2, "Smoothed ", ""), "AUC of roc2", sep = ""))
  } else {
    nest <- c(nest, paste(ifelse(attr(roc2$auc, "partial.auc.correct"), "Corrected ", ""),
      ifelse(smooth2, "Smoothed ", ""),
      "pAUC (", attr(roc2$auc, "partial.auc")[1], "-", attr(roc2$auc, "partial.auc")[2], " ", attr(roc2$auc, "partial.auc.focus"),
      ") of roc2",
      sep = ""
    ))
  }
  nest <- sub("Corrected Smoothed", "Corrected smoothed", nest) # no upper on smoothed if corrected.
  names(estimate) <- nest
  null.value <- 0
  names(null.value) <- "difference in AUC"
  htest <- list(
    alternative = alternative,
    data.names = data.names,
    data.name = data.names, # the standard htest element
    estimate = estimate,
    null.value = null.value
  )
  class(htest) <- "htest"

  if (method == "delong") {
    roc_test_check_conf_level(conf.level)
    if (paired) {
      delong.calcs <- delong.paired.calculations(roc1, roc2)
      stat <- c("Z" = delong.test(delong.calcs))
      htest$method <- "DeLong's test for two correlated ROC curves"
    } else {
      delong.calcs <- delong.unpaired.calculations(roc1, roc2)
      stat <- c("D" = delong.test(delong.calcs))
      htest$method <- "DeLong's test for two ROC curves"
    }
    htest$statistic <- stat
    htest$conf.int <- roc_test_ci_normal(delong.calcs$d, delong.calcs$sig, conf.level, alternative)
    if (roc1$percent) {
      # the placements are on the 0-1 scale, the estimates on 0-100
      htest$conf.int <- htest$conf.int * 100
    }

    if (alternative == "two.sided") {
      pval <- 2 * pnorm(-abs(stat))
    } else if (alternative == "greater") {
      pval <- pnorm(-stat)
    } else {
      pval <- pnorm(stat)
    }
    htest$p.value <- pval
  } else if (method == "venkatraman") {
    boot.n <- roc_utils_check_boot_n(boot.n)
    if (paired) {
      stats <- venkatraman.paired.test(roc1, roc2, boot.n, ties.method)
      htest$method <- "Venkatraman's test for two paired ROC curves"
    } else {
      stats <- venkatraman.unpaired.test(roc1, roc2, boot.n, ties.method)
      htest$method <- "Venkatraman's test for two unpaired ROC curves"
    }
    stat <- stats[[1]]
    names(stat) <- "E"
    htest$statistic <- stat
    parameter <- c(boot.n)
    names(parameter) <- "boot.n"
    htest$parameter <- parameter
    pval <- sum(stats[[2]] >= stats[[1]]) / boot.n
    htest$p.value <- pval
    names(htest$null.value) <- "difference in at least one ROC operating point"
    htest$estimate <- NULL # AUC not relevant in venkatraman
  } else { # method == "bootstrap" or "sensitivity" or "specificity"
    roc_test_check_conf_level(conf.level)
    boot.n <- roc_utils_check_boot_n(boot.n)
    if (method == "specificity") {
      if (!is.numeric(specificity) || length(specificity) != 1) {
        stop("Argument 'specificity' must be numeric of length 1 for a specificity test.")
      }
      boot <- bootstrap.test(roc1, roc2, "sp", specificity, paired, boot.n, boot.stratified, smoothing.args,
        progress = progress, cl = cl
      )
      if (paired) {
        htest$method <- "Specificity test for two correlated ROC curves"
      } else {
        htest$method <- "Specificity test for two ROC curves"
      }
      names(htest$null.value) <- sprintf(
        "difference in sensitivity at %s specificity",
        specificity
      )
      htest$estimate <- boot$estimate
      names(htest$estimate) <- sprintf("sensitivity of roc%d at %s specificity", 1:2, specificity)
    } else if (method == "sensitivity") {
      if (!is.numeric(sensitivity) || length(sensitivity) != 1) {
        stop("Argument 'sensitivity' must be numeric of length 1 for a sensitivity test.")
      }
      boot <- bootstrap.test(roc1, roc2, "se", sensitivity, paired, boot.n, boot.stratified, smoothing.args,
        progress = progress, cl = cl
      )
      if (paired) {
        htest$method <- "Sensitivity test for two correlated ROC curves"
      } else {
        htest$method <- "Sensitivity test for two ROC curves"
      }

      names(htest$null.value) <- sprintf(
        "difference in specificity at %s sensitivity",
        sensitivity
      )
      htest$estimate <- boot$estimate
      names(htest$estimate) <- sprintf("specificity of roc%d at %s sensitivity", 1:2, sensitivity)
    } else {
      boot <- bootstrap.test(roc1, roc2, "boot", NULL, paired, boot.n, boot.stratified, smoothing.args,
        progress = progress, cl = cl
      )
      if (paired) {
        htest$method <- "Bootstrap test for two correlated ROC curves"
      } else {
        htest$method <- "Bootstrap test for two ROC curves"
      }
    }
    stat <- c("D" = boot$statistic)
    htest$statistic <- stat
    htest$conf.int <- roc_test_ci_normal(boot$d, boot$sig, conf.level, alternative)
    parameter <- c(boot.n, boot.stratified)
    names(parameter) <- c("boot.n", "boot.stratified")
    htest$parameter <- parameter

    if (alternative == "two.sided") {
      pval <- 2 * pnorm(-abs(stat))
    } else if (alternative == "greater") {
      pval <- pnorm(-stat)
    } else {
      pval <- pnorm(stat)
    }
    htest$p.value <- pval
  }

  htest$roc1 <- roc1
  htest$roc2 <- roc2
  # Remove name from p value
  htest$p.value <- unname(htest$p.value)
  # Restore smoothing if necessary
  if (smooth1) {
    htest$roc1 <- roc_utils_resmooth(roc1, smoothing.args$roc1, reuse.auc = TRUE)
    # do.call() put the whole curve in the call: print the call of the input curve
    htest$roc1$call <- smooth.call1
  }
  if (smooth2) {
    htest$roc2 <- roc_utils_resmooth(roc2, smoothing.args$roc2, reuse.auc = TRUE)
    # do.call() put the whole curve in the call: print the call of the input curve
    htest$roc2$call <- smooth.call2
  }
  return(htest)
}

# Checks conf.level for the tests that return a CI
roc_test_check_conf_level <- function(conf.level) {
  # A logical counts as a number, as in t.test()
  if (!is.numeric(conf.level) && !is.logical(conf.level)) {
    stop("conf.level must be numeric between 0 and 1.")
  } else if (!roc_utils_is_single_number(conf.level)) {
    stop("conf.level must be a single number between 0 and 1.")
  } else if (0 > conf.level || 1 < conf.level) {
    stop("conf.level must be between 0 and 1.")
  }
}

# Normal CI of the difference d with standard error sig, the inverse of the
# Z test. One-sided alternatives give a one-sided interval, as in t.test.
roc_test_ci_normal <- function(d, sig, conf.level, alternative) {
  crit_z <- qnorm(if (alternative == "two.sided") 1 - (1 - conf.level) / 2 else conf.level)
  conf.int <- c(
    if (alternative == "less") -Inf else d - crit_z * sig,
    if (alternative == "greater") Inf else d + crit_z * sig
  )
  attr(conf.int, "conf.level") <- conf.level
  conf.int
}
