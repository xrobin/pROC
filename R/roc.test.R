# pROC: Tools Receiver operating characteristic (ROC curves) with
# (partial) area under the curve, confidence intervals and comparison. 
# Copyright (C) 2010, 2011 Xavier Robin, Alexandre Hainard, Natacha Turck,
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

roc.test.formula <- function (formula, data, ...){
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  term.labels <- attr(attr(m, "terms"), "term.labels")
  response <- model.extract(m, "response")
  
  if (length(term.labels) != 2) {
    stop("Invalid formula: exactly 2 predictors are required in a formula of type response~predictor1+predictor2.")
  }
  if (length(response) == 0) {
    stop("Error in the formula: a response is required in a formula of type response~predictor1+predictor2.")
  }

  testres <- roc.test.default(response, m[[term.labels[1]]], m[[term.labels[2]]], ...)
  testres$call <- cl
  # data.names for pretty print()ing
  if (missing(data))
    testres$data.names <- sprintf("%s and %s by %s (%s, %s)", term.labels[1], term.labels[2], names(m)[1], testres$roc1$levels[1], testres$roc1$levels[2])
  else
    testres$data.names <- sprintf("%s and %s in %s by %s (%s, %s)", term.labels[1], term.labels[2], deparse(substitute(data)), names(m)[1], testres$roc1$levels[1], testres$roc1$levels[2])

  return(testres)
}

roc.test.default <- function(response, predictor1, predictor2=NULL, na.rm=TRUE, method=NULL, ...) {
  if (is.matrix(predictor1) | is.data.frame(predictor1)) {
    if (!is.null(predictor2))
      stop("Predictor2 must not be specified if predictor1 is a matrix or a data.frame.")
    if (dim(predictor1)[2] == 2 & length(response) == dim(predictor1)[1]) {
      roc1 <- roc(response, predictor1[,1], ...)
      roc2 <- roc(response, predictor1[,2], ...)
      if (!is.null(names(predictor1)))
        data.names <- sprintf("%s and %s in %s by %s (%s, %s)", names(predictor1)[1], names(predictor1)[2], deparse(substitute(predictor1)), deparse(substitute(response)), roc1$levels[1], roc1$levels[2])
      else if (!is.null(colnames(predictor1)))
        data.names <- sprintf("%s and %s in %s by %s (%s, %s)", colnames(predictor1)[1], colnames(predictor1)[2], deparse(substitute(predictor1)), deparse(substitute(response)), roc1$levels[1], roc1$levels[2])
      else
        data.names <- sprintf("%s by %s (%s, %s)", deparse(substitute(predictor1)), deparse(substitute(response)), roc1$levels[1], roc1$levels[2])
    }
    else {
      stop("Wrong dimension for predictor1 as a matrix or a data.frame.")
    }
  }
  else {
    if (missing(predictor2))
      stop("Missing argument predictor2 with predictor1 as a vector.")
    # Need to remove NAs
    if (na.rm) {
      nas <- is.na(response) | is.na(predictor1) | is.na(predictor2)
      response <- response[!nas]
      predictor1 <- predictor1[!nas]
      predictor2 <- predictor2[!nas]
    }
    roc1 <- roc(response, predictor1, ...)
    roc2 <- roc(response, predictor2, ...)
    call <- match.call()
    data.names <- sprintf("%s and %s by %s (%s, %s)", deparse(call$predictor1), deparse(call$predictor2), deparse(call$response), roc1$levels[1], roc1$levels[2])
  }
  test <- roc.test.roc(roc1, roc2, method=method, ...)
  test$data.names <- data.names
  return(test)
}

roc.test.auc <- function(roc1, roc2, ...) {
  testres <- roc.test.roc(attr(roc1, "roc"), roc2, ...)
  testres$call <- match.call()
  testres$data.names <- paste(deparse(substitute(roc1)), "and", deparse(substitute(roc2)))
  return(testres)
}

roc.test.smooth.roc <- function(roc1, roc2, ...) {
  testres <- roc.test.roc(roc1, roc2, ...)
  testres$call <- match.call()
  testres$data.names <- paste(deparse(substitute(roc1)), "and", deparse(substitute(roc2)))
  return(testres)
}

roc.test.roc <- function(roc1, roc2,
                         method=c("delong", "bootstrap", "venkatraman", "sensitivity", "specificity"),
                         sensitivity=NULL, specificity=NULL,
                         alternative = c("two.sided", "less", "greater"),
                         paired=NULL,
                         reuse.auc=TRUE,
                         boot.n=2000, boot.stratified=TRUE,
                         ties.method="first",
                         progress=getOption("pROCProgress")$name,
                         ...) {
  alternative <- match.arg(alternative)
  data.names <- paste(deparse(substitute(roc1)), "and", deparse(substitute(roc2)))
  if (class(roc2) == "auc")
    roc2 <- attr(roc2, "roc")

  # store which objects are smoothed, and how
  smoothing.args <- list()
  if (class(roc1) == "smooth.roc") {
    smoothing.args$roc1 <- roc1$smoothing.args
    smoothing.args$roc1$smooth <- TRUE
    roc1 <- attr(roc1, "roc")
  }
  else {
    smoothing.args$roc1 <- list(smooth=FALSE)
  }
  if (class(roc2) == "smooth.roc") {
    smoothing.args$roc2 <- roc2$smoothing.args
    smoothing.args$roc2$smooth <- TRUE
    roc2 <- attr(roc2, "roc")
  }
  else {
    smoothing.args$roc2 <- list(smooth=FALSE)
  }

  # Check if we do a paired or unpaired roc.test
  if (is.null(paired)) {
    # then determine whether the rocs are paired or not
    rocs.are.paired <- are.paired(roc1, roc2, return.paired.rocs=TRUE, reuse.auc=TRUE, reuse.ci=FALSE, reuse.smooth=TRUE)
    if (rocs.are.paired) {
      paired <- TRUE
      roc1 <- attr(rocs.are.paired, "roc1")
      roc2 <- attr(rocs.are.paired, "roc2")
    }
    else {
      paired <- FALSE
      roc1 <- roc1
      roc2 <- roc2
    }
  }
  else if (paired) {
    # make sure the rocs are really paired
    rocs.are.paired <- rocs.are.paired <- are.paired(roc1, roc2, return.paired.rocs=TRUE, reuse.auc=TRUE, reuse.ci=FALSE, reuse.smooth=TRUE)
    if (! rocs.are.paired) 
      stop("The paired ROC test cannot be applied to unpaired curves.")
    roc1 <- attr(rocs.are.paired, "roc1")
    roc2 <- attr(rocs.are.paired, "roc2")
  }
  else { # assume unpaired
    rocs.are.paired <- are.paired(roc1, roc2, return.paired.rocs=FALSE)
    if (rocs.are.paired) 
      warning("The ROC curves seem to be paired. Consider performing a paired roc.test.")
    roc1 <- roc1
    roc2 <- roc2
  }

  # check that the AUC was computed, or do it now
  if (is.null(roc1$auc) | !reuse.auc) {
    if (smoothing.args$roc1$smooth) {
      roc1$auc <- auc(smooth.roc=do.call("smooth.roc", c(list(roc=roc1), smoothing.args$roc1)), ...)
      # remove partial.auc.* arguments that are now in roc1$auc and that will mess later processing
      # (formal argument "partial.auc(.*)" matched by multiple actual arguments)
      # This removal should be safe because we always use smoothing.args with roc1 in the following processing,
      # however it is a potential source of bugs.
      smoothing.args$roc1$partial.auc <- NULL
      smoothing.args$roc1$partial.auc.correct <- NULL
      smoothing.args$roc1$partial.auc.focus <- NULL
    }
    else
      roc1$auc <- auc(roc1, ...)
  }
  if (is.null(roc2$auc) | !reuse.auc) {
    if (smoothing.args$roc2$smooth) {
      roc2$auc <- auc(smooth.roc=do.call("smooth.roc", c(list(roc=roc2), smoothing.args$roc2)), ...)
      # remove partial.auc.* arguments that are now in roc1$auc and that will mess later processing
      # (formal argument "partial.auc(.*)" matched by multiple actual arguments)
      # This removal should be safe because we always use smoothing.args with roc2 in the following processing,
      # however it is a potential source of bugs.
      smoothing.args$roc2$partial.auc <- NULL
      smoothing.args$roc2$partial.auc.correct <- NULL
      smoothing.args$roc2$partial.auc.focus <- NULL
    }
    else
      roc2$auc <- auc(roc2, ...)
  }
    
  # check that the same region was requested in auc. Otherwise, issue a warning
  if (!identical(attributes(roc1$auc)[names(attributes(roc1$auc))!="roc"], attributes(roc2$auc)[names(attributes(roc2$auc))!="roc"]))
    warning("Different AUC specifications in the ROC curves. Enforcing the inconsistency, but unexpected results may be produced.")
  # check that the same smoothing params were requested in auc. Otherwise, issue a warning
  if (!identical(smoothing.args$roc1, smoothing.args$roc2))
    warning("Different smoothing parameters in the ROC curves. Enforcing the inconsistency, but unexpected results may be produced.")

  # Check the method
  if (missing(method) | is.null(method)) {
    # determine method if missing
    if (is.numeric(attr(roc1$auc, "partial.auc")) && length(attr(roc1$auc, "partial.auc") == 2)) {
      # partial auc: go for bootstrap
      method <- "bootstrap"
    }
    else if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth) {
      # smoothing in one or both: bootstrap
      method <- "bootstrap"
    }
    else if (roc1$direction != roc2$direction) {
      # delong doesn't work well with opposite directions (will report high significance if roc1$auc and roc2$auc are similar and high)
      method <- "bootstrap"
    }
    else {
      method <- "delong"
    }
  }
  else {
    method <- match.arg(method)
    if (method == "delong") {
      # delong NA to pAUC: warn + change
      if ((is.numeric(attr(roc1$auc, "partial.auc")) && length(attr(roc1$auc, "partial.auc") == 2)) || (is.numeric(attr(roc2$auc, "partial.auc")) && length(attr(roc2$auc, "partial.auc") == 2))) {
        warning("Using DeLong's test for partial AUC is not supported. Using bootstrap test instead.")
        method <- "bootstrap"
      }
      if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth) {
        warning("Using DeLong's test for smoothed ROCs is not supported. Using bootstrap test instead.")
        method <- "bootstrap"
      }
      if (roc1$direction != roc2$direction)
        warning("DeLong's test should not be applied to ROC curves with a different direction.")
    }
    else if (method == "venkatraman") {
      if (is.numeric(attr(roc1$auc, "partial.auc")) && length(attr(roc1$auc, "partial.auc") == 2))
        warning("Partial AUC is ignored in Venkatraman's test.")
      if (smoothing.args$roc1$smooth || smoothing.args$roc2$smooth)
        stop("Using Venkatraman's test for smoothed ROCs is not supported.")
      if (roc1$direction != roc2$direction)
        warning("Venkatraman's test should not be applied to ROC curves with different directions.")
      if (alternative != "two.sided") {
        warning("Only two-sided tests are available for Venkatraman. Performing two-sided test instead.")
        alternative <- "two.sided"
      }
    }
  }

  # Prepare the return value htest
  if (smoothing.args$roc1$smooth)
    estimate <- do.call("smooth.roc", c(list(roc=roc1), smoothing.args$roc1))$auc
  else
    estimate <- roc1$auc
  if (smoothing.args$roc2$smooth)
    estimate <- c(estimate, do.call("smooth.roc", c(list(roc=roc2), smoothing.args$roc2))$auc)
  else
    estimate <- c(estimate, roc2$auc)
  if (identical(attr(roc1$auc, "partial.auc"), FALSE)) {
    nest <- paste(ifelse(smoothing.args$roc1$smooth, "Smoothed ", ""), "AUC of roc1", sep="")
  }
  else {
    nest <- paste(ifelse (attr(roc1$auc, "partial.auc.correct"), "Corrected ", ""),
                  ifelse (smoothing.args$roc1$smooth, "Smoothed ", ""),
                  "pAUC (", attr(roc1$auc, "partial.auc")[1], "-", attr(roc1$auc, "partial.auc")[2], " ", attr(roc1$auc, "partial.auc.focus"),
                        ") of roc1", sep="")
  }
  if (identical(attr(roc2$auc, "partial.auc"), FALSE)) {
    nest <- c(nest, paste(ifelse(smoothing.args$roc2$smooth, "Smoothed ", ""), "AUC of roc2", sep=""))
  }
  else {
    nest <- c(nest, paste(ifelse (attr(roc2$auc, "partial.auc.correct"), "Corrected ", ""),
                          ifelse (smoothing.args$roc2$smooth, "Smoothed ", ""),
                          "pAUC (", attr(roc2$auc, "partial.auc")[1], "-", attr(roc2$auc, "partial.auc")[2], " ", attr(roc2$auc, "partial.auc.focus"),
                          ") of roc2", sep=""))
  }
  nest <- sub("Corrected Smoothed", "Corrected smoothed", nest) # no upper on smoothed if corrected.
  names(estimate) <- nest
  null.value <- 0
  names(null.value) <- "difference in AUC"
  htest <- list(
                alternative = alternative,
                data.names = data.names,
                estimate = estimate,
                null.value = null.value
                )
  class(htest) <- "htest"
  
  if (method == "delong") {
    if (paired) {
      stat <- delong.paired.test(roc1, roc2)
      names(stat) <- "Z"
      htest$statistic <- stat
      htest$method <- "DeLong's test for two correlated ROC curves"
      
      if (alternative == "two.sided")
        pval <- 2*pnorm(-abs(stat))
      else if (alternative == "greater")
        pval <- pnorm(-stat)
      else
        pval <- pnorm(stat)
      htest$p.value <- pval
    }
    else {
      stats <- delong.unpaired.test(roc1, roc2)
      stat <- stats[1]
      df <- stats[2]
      htest$statistic <- c("D"=stat)
      htest$parameter <- c("df"=df)
      htest$method <- "DeLong's test for two ROC curves"

      if (alternative == "two.sided")
        pval <- 2*pt(-abs(stat), df=df)
      else if (alternative == "greater")
        pval <- pt(-stat, df=df)
      else
        pval <- pt(stat, df=df)
      htest$p.value <- pval
    }
  }
  else if (method == "venkatraman") {
    if(class(progress) != "list")
      progress <- roc.utils.get.progress.bar(progress, title="Venkatraman ROC test", label="Permutations in progress...", ...)
    if (paired) {
      stats <- venkatraman.paired.test(roc1, roc2, boot.n, ties.method, progress)
      htest$method <- "Venkatraman's test for two paired ROC curves"
    }
    else {
      stats <- venkatraman.unpaired.test(roc1, roc2, boot.n, ties.method, progress)
      htest$method <- "Venkatraman's test for two unpaired ROC curves"
    }
    stat <- stats[[1]]
    names(stat) <- "E"
    htest$statistic <- stat
    parameter <- c(boot.n)
    names(parameter) <- "boot.n"
    htest$parameter <- parameter
    pval <- sum(stats[[2]]>=stats[[1]])/boot.n
    htest$p.value <- pval
    names(null.value) <- "difference in ROC operating points"
    htest$estimate <- NULL # AUC not relevant in venkatraman
  }
  else { # method == "bootstrap" or "sensitivity" or "specificity"
    # Check if called with density.cases or density.controls
    if (is.null(smoothing.args) || is.numeric(smoothing.args$density.cases) || is.numeric(smoothing.args$density.controls))
      stop("Cannot compute the statistic on ROC curves smoothed with numeric density.controls and density.cases.")

    if(class(progress) != "list")
      progress <- roc.utils.get.progress.bar(progress, title="Bootstrap ROC test", label="Bootstrap in progress...", ...)

    if (method == "specificity") {
      if (! is.numeric(specificity) || length(specificity) != 1) {
        stop("Argument 'specificity' must be numeric of length 1 for a specificity test.")
      }
      stat <- bootstrap.test(roc1, roc2, "sp", specificity, paired, boot.n, boot.stratified, smoothing.args, progress)
      if (paired)
        htest$method <- "Specificity test for two correlated ROC curves"
      else
        htest$method <- "Specificity test for two ROC curves"
    }
    else if (method == "sensitivity") {
      if (! is.numeric(sensitivity) || length(sensitivity) != 1) {
        stop("Argument 'sensitivity' must be numeric of length 1 for a sensitivity test.")
      }
      stat <- bootstrap.test(roc1, roc2, "se", sensitivity, paired, boot.n, boot.stratified, smoothing.args, progress)
      if (paired)
        htest$method <- "Sensitivity test for two correlated ROC curves"
      else
        htest$method <- "Sensitivity test for two ROC curves"
    }
    else {
      stat <- bootstrap.test(roc1, roc2, "boot", NULL, paired, boot.n, boot.stratified, smoothing.args, progress)
      if (paired)
        htest$method <- "Bootstrap test for two correlated ROC curves"
      else
        htest$method <- "Bootstrap test for two ROC curves"
    }
    stat <- as.vector(stat) # remove auc attributes
    names(stat) <- "D"
    htest$statistic <- stat
    parameter <- c(boot.n, boot.stratified)
    names(parameter) <- c("boot.n", "boot.stratified")
    htest$parameter <- parameter
    
    if (alternative == "two.sided")
      pval <- 2*pnorm(-abs(stat))
    else if (alternative == "greater")
      pval <- pnorm(-stat)
    else
      pval <- pnorm(stat)
    htest$p.value <- pval
  }

  htest$roc1 <- roc1
  htest$roc2 <- roc2
  # Restore smoothing if necessary
  if (smoothing.args$roc1$smooth)
    htest$roc1 <- do.call("smooth.roc", c(list(roc=roc1), smoothing.args$roc1))
  if (smoothing.args$roc2$smooth)
    htest$roc2 <- do.call("smooth.roc", c(list(roc=roc2), smoothing.args$roc2))
  return(htest)
}

venkatraman.paired.test <- function(roc1, roc2, boot.n, ties.method="first", progress) {
  X <- roc1$predictor
  Y <- roc2$predictor
  R <- rank(X, ties.method = ties.method)
  S <- rank(Y, ties.method = ties.method)
  D <- roc1$response # because roc1&roc2 are paired

  E <- venkatraman.paired.stat(R, S, D, roc1$levels)
  EP <- raply(boot.n, venkatraman.paired.permutation(R, S, D, roc1$levels, ties.method), .progress=progress)
  return(list(E, EP))
}

venkatraman.unpaired.test <- function(roc1, roc2, boot.n, ties.method="first", progress) {
  X <- roc1$predictor
  Y <- roc2$predictor
  R <- rank(X, ties.method = ties.method)
  S <- rank(Y, ties.method = ties.method)
  D1<- roc1$response
  D2 <- roc2$response
  mp <- (sum(D1 == roc1$levels[2]) + sum(D2 == roc2$levels[2]))/(length(D1) + length(D1)) # mixing proportion, kappa

  E <- venkatraman.unpaired.stat(R, S, D1, D2, roc1$levels, roc2$levels, mp)
  EP <- raply(boot.n, venkatraman.unpaired.permutation(R, S, D1, D2, roc1$levels, roc2$levels, mp, ties.method), .progress=progress)
  return(list(E, EP))
}


venkatraman.paired.permutation <- function(R, S, D, levels, ties.method) {
  # Break ties
  R2 <- R + runif(length(D)) - 0.5 # Add small amount of random but keep same mean
  S2 <- S + runif(length(D)) - 0.5

  # Permutation
  q <- 1 - round(runif(length(D)))
  R3 <- R2 * q + (1 - q) * S
  S3 <- S2 * q + (1 - q) * R

  return(venkatraman.paired.stat(rank(R3, ties.method=ties.method), rank(S3, ties.method=ties.method), D, levels))
}


venkatraman.unpaired.permutation <- function(R, S, D1, D2, levels1, levels2, mp, ties.method) {
  # Break ties
  R <- R + runif(length(D1)) - 0.5 # Add small amount of random but keep same mean
  S <- S + runif(length(D2)) - 0.5

  R.controls <- R[D1==levels1[1]]
  R.cases <- R[D1==levels1[2]]
  S.controls <- S[D2==levels2[1]]
  S.cases <- S[D2==levels2[2]]

  # Permutation
  controls <- sample(c(R.controls, S.controls))
  cases <- sample(c(R.cases, S.cases))
  R[D1==levels1[1]] <- controls[1:length(R.controls)]
  S[D2==levels2[1]] <- controls[(length(R.controls)+1):length(controls)]
  R[D1==levels1[2]] <- cases[1:length(R.cases)]
  S[D2==levels2[2]] <- cases[(length(R.cases)+1):length(cases)]

  return(venkatraman.unpaired.stat(rank(R, ties.method=ties.method), rank(S, ties.method=ties.method), D1, D2, levels1, levels2, mp))
}

venkatraman.paired.stat <- function(R, S, D, levels) {
  R.controls <- R[D==levels[1]]
  R.cases <- R[D==levels[2]]
  S.controls <- S[D==levels[1]]
  S.cases <- S[D==levels[2]]
  n <- length(D)

  R.fn <- sapply(1:n, function(x) sum(R.cases <= x))
  R.fp <- sapply(1:n, function(x) sum(R.controls > x))
  S.fn <- sapply(1:n, function(x) sum(S.cases <= x))
  S.fp <- sapply(1:n, function(x) sum(S.controls > x))

  return(sum(abs((S.fn + S.fp) - (R.fn + R.fp))))
}

venkatraman.unpaired.stat <- function(R, S, D1, D2, levels1, levels2, mp) {
  R.controls <- R[D1==levels1[1]]
  R.cases <- R[D1==levels1[2]]
  S.controls <- S[D2==levels2[1]]
  S.cases <- S[D2==levels2[2]]
  n <- length(D1)
  m <- length(D2)

  R.fx <- sapply(1:n, function(x) sum(R.cases <= x)) / length(R.cases)
  R.gx <- sapply(1:n, function(x) sum(R.controls <= x)) / length(R.controls)
  S.fx <- sapply(1:m, function(x) sum(S.cases <= x)) / length(S.cases)
  S.gx <- sapply(1:m, function(x) sum(S.controls <= x)) / length(S.controls)
  R.p <- mp*R.fx + (1 - mp)*R.gx
  S.p <- mp*S.fx + (1 - mp)*S.gx
  R.exp <- mp*R.fx + (1 - mp)*(1-R.gx)
  S.exp <- mp*S.fx + (1 - mp)*(1-S.gx)

  # Do the integration
  x <- sort(c(R.p, S.p))
  R.f <- approxfun(R.p, R.exp)
  S.f <- approxfun(S.p, S.exp)
  f  <- function(x) abs(R.f(x)-S.f(x))
  y <- f(x)
  #trapezoid integration:
  idx <- 2:length(x)
  integral <- sum(((y[idx] + y[idx-1]) * (x[idx] - x[idx-1])) / 2, na.rm=TRUE) # remove NA that can appear in the borders
  return(integral)
}

# Delong's test paired, used by roc.test.roc
delong.paired.test <- function(roc1, roc2) {

  n <- length(roc1$controls)
  m <- length(roc1$cases)

  VR <- roc.utils.delong.placements(roc1)
  VS <- roc.utils.delong.placements(roc2)

  SX <- matrix(NA, ncol=2, nrow=2)
  SX[1,1] <- sum((VR$X - VR$theta) * (VR$X - VR$theta))/(m-1)
  SX[1,2] <- sum((VR$X - VR$theta) * (VS$X - VS$theta))/(m-1)
  SX[2,1] <- sum((VS$X - VS$theta) * (VR$X - VR$theta))/(m-1)
  SX[2,2] <- sum((VS$X - VS$theta) * (VS$X - VS$theta))/(m-1)

  
  SY <- matrix(NA, ncol=2, nrow=2)
  SY[1,1] <- sum((VR$Y - VR$theta) * (VR$Y - VR$theta))/(n-1)
  SY[1,2] <- sum((VR$Y - VR$theta) * (VS$Y - VS$theta))/(n-1)
  SY[2,1] <- sum((VS$Y - VS$theta) * (VR$Y - VR$theta))/(n-1)
  SY[2,2] <- sum((VS$Y - VS$theta) * (VS$Y - VS$theta))/(n-1)

  S <- SX/m + SY/n
  L <- c(1,-1)
  sig <- sqrt(L%*%S%*%L)
  zscore <- (VR$theta-VS$theta)/sig[1]
  if (is.nan(zscore) && VR$theta == VR$theta && sig[1] == 0)
    zscore <- 0 # special case: no difference between theta's produces a NaN
  return(zscore)
}

# Delong's test unpaired, used by roc.test.roc
delong.unpaired.test <- function(roc1, roc2) {

  nR <- length(roc1$controls)
  mR <- length(roc1$cases)

  nS <- length(roc2$controls)
  mS <- length(roc2$cases)

  VR <- roc.utils.delong.placements(roc1)
  VS <- roc.utils.delong.placements(roc2)

  SRX <- sum((VR$X - VR$theta) * (VR$X - VR$theta))/(mR-1)
  SSX <- sum((VS$X - VS$theta) * (VS$X - VS$theta))/(mS-1)

  SRY <- sum((VR$Y - VR$theta) * (VR$Y - VR$theta))/(nR-1)
  SSY <- sum((VS$Y - VS$theta) * (VS$Y - VS$theta))/(nS-1)

  SR <- SRX/mR + SRY/nR
  SS <- SSX/mS + SSY/nS

  ntotR <- nR + mR
  ntotS <- nS + mS
  SSR <- sqrt((SR) + (SS))
  t <- (VR$theta - VS$theta) / SSR
  df <- ((SR) + (SS))^2 /
    (((SR)^2 / (ntotR-1)) + ((SS)^2 / (ntotS -1 )))

  return(c(t, df))
}

# Paired bootstrap test, used by roc.test.roc
bootstrap.test <- function(roc1, roc2, test, x, paired, boot.n, boot.stratified, smoothing.args, progress) {

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
  auc1skeleton <- c(auc1skeleton, smoothing.args$roc1)
  auc2skeleton <- attributes(roc2$auc)
  auc2skeleton$roc <- NULL
  auc2skeleton$direction <- roc2$direction
  auc2skeleton$class <- NULL
  auc2skeleton <- c(auc2skeleton, smoothing.args$roc2)

  # Some attributes may be duplicated in AUC skeletons and will mess the boostrap later on when we do.call().
  # If this condition happen, it probably means we have a bug elsewhere.
  # Rather than making a complicated processing to remove the duplicates,
  # just throw an error and let us solve the bug when a user reports it.
  duplicated.auc1skeleton <- duplicated(names(auc1skeleton))
  duplicated.auc2skeleton <- duplicated(names(auc2skeleton))
  if (any(duplicated.auc1skeleton))
    stop(sprintf("duplicated argument(s) in AUC1 skeleton: \"%s\". Please report this bug to the package maintainer %s", paste(names(auc1skeleton)[duplicated(names(auc1skeleton))], collapse=", "), packageDescription("pROC")$Maintainer))
  if (any(duplicated.auc2skeleton))
    stop(sprintf("duplicated argument(s) in AUC2 skeleton: \"%s\". Please report this bug to the package maintainer %s", paste(names(auc2skeleton)[duplicated(names(auc2skeleton))], collapse=", "), packageDescription("pROC")$Maintainer))
#  # This was an attempt to remove duplicated arguments and check their equality. Now we die just above if it happens
#  if (any(duplicated.auc1skeleton)) {
#    #store and remove duplicates
#    auc1skeleton.removed <- auc1skeleton[names(auc1skeleton)[duplicated.auc1skeleton]] # store
#    auc1skeleton[names(auc1skeleton)[duplicated.auc1skeleton]] <- NULL # remove
#    auc1skeleton.kept <- auc1skeleton[names(auc1skeleton.removed)] # check what is left
#    if (!all.equal(auc1skeleton.kept, auc1skeleton.removed))
#      stop("conflicting arguments in AUC1 skeleton")
#  }
#  if (any(duplicated.auc2skeleton)) {
#    #store and remove duplicates
#    auc2skeleton.removed <- auc2skeleton[names(auc2skeleton)[duplicated.auc2skeleton]] # store
#    auc2skeleton[names(auc2skeleton)[duplicated.auc2skeleton]] <- NULL # remove
#    auc2skeleton.kept <- auc2skeleton[names(auc2skeleton.removed)] # check what is left
#    if (!all.equal(auc2skeleton.kept, auc2skeleton.removed))
#      stop("conflicting arguments in AUC2 skeleton")
#  }

  if (boot.stratified) { # precompute sorted responses if stratified
    response.roc1 <- factor(c(rep(roc1$levels[1], length(roc1$controls)), rep(roc1$levels[2], length(roc1$cases))), levels=roc1$levels)
    response.roc2 <- factor(c(rep(roc2$levels[1], length(roc2$controls)), rep(roc2$levels[2], length(roc2$cases))), levels=roc2$levels)
    auc1skeleton$response <- response.roc1
    auc2skeleton$response <- response.roc2
    resampled.values <- raply(boot.n, stratified.bootstrap.test(roc1, roc2, test, x, paired, auc1skeleton, auc2skeleton), .progress=progress)
  }
  else {
    resampled.values <- raply(boot.n, nonstratified.bootstrap.test(roc1, roc2, test, x, paired, auc1skeleton, auc2skeleton), .progress=progress)
  }

  # compute the statistics
  diffs <- resampled.values[,1] - resampled.values[,2]

  # are there NA values?
  if ((num.NAs <- sum(is.na(diffs))) > 0) {
    warning(sprintf("%i NA value(s) produced during bootstrap were ignored.", num.NAs))
    diffs <- diffs[!is.na(diffs)]
  }

  # Restore smoothing if necessary
  if (smoothing.args$roc1$smooth) {
    smoothing.args$roc1$method <- smoothing.args$roc1$smooth.method
    roc1 <- do.call("smooth.roc", c(list(roc=roc1), smoothing.args$roc1))
  }
  if (smoothing.args$roc2$smooth) {
    smoothing.args$roc2$method <- smoothing.args$roc2$smooth.method
    roc2 <- do.call("smooth.roc", c(list(roc=roc2), smoothing.args$roc2))
  }

  if (test == "sp") {
    coord1 <- coords(roc1, x=x, input=c("specificity"), ret=c("sensitivity"), as.list=FALSE)
    coord2 <- coords(roc2, x=x, input=c("specificity"), ret=c("sensitivity"), as.list=FALSE)
    D <- (coord1 - coord2) / sd(diffs)
  }
  else if (test == "se") {
    coord1 <- coords(roc1, x=x, input=c("sensitivity"), ret=c("specificity"), as.list=FALSE)
    coord2 <- coords(roc2, x=x, input=c("sensitivity"), ret=c("specificity"), as.list=FALSE)
    D <- (coord1 - coord2) / sd(diffs)
  }
  else {
    D <- (roc1$auc - roc2$auc) / sd(diffs)
  }
  if (is.nan(D) && all(diffs == 0) && roc1$auc == roc2$auc)
    D <- 0 # special case: no difference between AUCs produces a NaN

  return(D)
}

stratified.bootstrap.test <- function(roc1, roc2, test, x, paired, auc1skeleton, auc2skeleton) {
  # sample control and cases separately for a stratified bootstrap
  idx.controls.roc1 <- sample(1:length(roc1$controls), replace=TRUE)
  idx.cases.roc1 <- sample(1:length(roc1$cases), replace=TRUE)
  # finish roc skeletons
  auc1skeleton$predictor <- c(roc1$controls[idx.controls.roc1], roc1$cases[idx.cases.roc1])
  if (paired) {
    auc2skeleton$predictor <- c(roc2$controls[idx.controls.roc1], roc2$cases[idx.cases.roc1])
  }
  else { # for unpaired, resample roc2 separately
    idx.controls.roc2 <- sample(1:length(roc2$controls), replace=TRUE)
    idx.cases.roc2 <- sample(1:length(roc2$cases), replace=TRUE)
    auc2skeleton$predictor <- c(roc2$controls[idx.controls.roc2], roc2$cases[idx.cases.roc2])
  }

  # re-compute the resampled ROC curves
  roc1 <- try(do.call("roc.default", auc1skeleton), silent=TRUE)
  roc2 <- try(do.call("roc.default", auc2skeleton), silent=TRUE)
  # resampled ROCs might not be smoothable: return NA
  if (class(roc1) == "try-error" || class(roc2) == "try-error") {
    return(c(NA, NA))
  }
  else {
    if (test == "sp") {
      coord1 <- coords(roc1, x=x, input=c("specificity"), ret=c("sensitivity"), as.list=FALSE)
      coord2 <- coords(roc1, x=x, input=c("specificity"), ret=c("sensitivity"), as.list=FALSE)
      return(c(coord1, coord2))
    }
    else if (test == "se") {
      coord1 <- coords(roc1, x=x, input=c("sensitivity"), ret=c("specificity"), as.list=FALSE)
      coord2 <- coords(roc1, x=x, input=c("sensitivity"), ret=c("specificity"), as.list=FALSE)
      return(c(coord1, coord2))
    }
    else {
      return(c(roc1$auc, roc2$auc))
    }
  }
}

nonstratified.bootstrap.test <- function(roc1, roc2, test, x, paired, auc1skeleton, auc2skeleton) {
  # sample all patients
  idx.all.roc1 <- sample(1:length(roc1$response), replace=TRUE)
  # finish roc skeletons
  auc1skeleton$response <- roc1$response[idx.all.roc1]
  auc1skeleton$predictor <- roc1$predictor[idx.all.roc1]
  if (paired) { # if paired, resample roc2 as roc1
    auc2skeleton$response <- roc2$response[idx.all.roc1]
    auc2skeleton$predictor <- roc2$predictor[idx.all.roc1]
  }
  else { # if unpaired, resample roc2 separately
    idx.all.roc2 <- sample(1:length(roc2$response), replace=TRUE)
    auc2skeleton$response <- roc2$response[idx.all.roc2]
    auc2skeleton$predictor <- roc2$predictor[idx.all.roc2]
  }

  # re-compute the resampled ROC curves
  roc1 <- try(do.call("roc.default", auc1skeleton), silent=TRUE)
  roc2 <- try(do.call("roc.default", auc2skeleton), silent=TRUE)
  # resampled ROCs might not be smoothable: return NA
  if (class(roc1) == "try-error" || class(roc2) == "try-error") {
    return(c(NA, NA))
  }
  else {
    if (test == "sp") {
      coord1 <- coords(roc1, x=x, input=c("specificity"), ret=c("sensitivity"), as.list=FALSE)
      coord2 <- coords(roc1, x=x, input=c("specificity"), ret=c("sensitivity"), as.list=FALSE)
      return(c(coord1, coord2))
    }
    else if (test == "se") {
      coord1 <- coords(roc1, x=x, input=c("sensitivity"), ret=c("specificity"), as.list=FALSE)
      coord2 <- coords(roc1, x=x, input=c("sensitivity"), ret=c("specificity"), as.list=FALSE)
      return(c(coord1, coord2))
    }
    else {
      return(c(roc1$auc, roc2$auc))
    }
  }
}
