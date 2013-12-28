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

##########  AUC of two ROC curves (roc.test, cov)  ##########

bootstrap.cov <- function(roc1, roc2, boot.n, boot.stratified, boot.return, smoothing.args, progress) {
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

  if (boot.stratified) { # precompute sorted responses if stratified
    response.roc1 <- factor(c(rep(roc1$levels[1], length(roc1$controls)), rep(roc1$levels[2], length(roc1$cases))), levels=roc1$levels)
    response.roc2 <- factor(c(rep(roc2$levels[1], length(roc2$controls)), rep(roc2$levels[2], length(roc2$cases))), levels=roc2$levels)
    auc1skeleton$response <- response.roc1
    auc2skeleton$response <- response.roc2
    resampled.values <- raply(boot.n, stratified.bootstrap.test(roc1, roc2, "boot", NULL, TRUE, auc1skeleton, auc2skeleton), .progress=progress)
  }
  else {
    resampled.values <- raply(boot.n, nonstratified.bootstrap.test(roc1, roc2, "boot", NULL, TRUE, auc1skeleton, auc2skeleton), .progress=progress)
  }

  # are there NA values?
  if ((num.NAs <- sum(apply(resampled.values, 1, is.na))) > 0) {
    warning(sprintf("%i NA value(s) produced during bootstrap were ignored.", num.NAs))
    resampled.values <- resampled.values[!apply(resampled.values, 1, function(x) any(is.na(x))),]
  }

  cov <- stats::cov(resampled.values[,1], resampled.values[,2])
  if (boot.return) {
    attr(cov, "resampled.values") <- resampled.values
  }
  return(cov)
}

# Bootstrap test, used by roc.test.roc
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
  if (is(roc1, "try-error") || is(roc2, "try-error")) {
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
  if (is(roc1, "try-error") || is(roc2, "try-error")) {
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

##########  AUC of one ROC curves (ci.auc, var)  ##########

ci.auc.bootstrap <- function(roc, conf.level, boot.n, boot.stratified, progress, ...) {
  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="AUC confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    aucs <- unlist(rlply(boot.n, stratified.ci.auc(roc), .progress=progress))
  }
  else {
    aucs <- unlist(rlply(boot.n, nonstratified.ci.auc(roc), .progress=progress))
  }

  if (sum(is.na(aucs)) > 0) {
    warning("NA value(s) produced during bootstrap were ignored.")
    aucs <- aucs[!is.na(aucs)]
  }
  # TODO: Maybe apply a correction (it's in the Tibshirani?) What do Carpenter-Bithell say about that?
  # Prepare the return value
  return(quantile(aucs, c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))
}

stratified.ci.auc <- function(roc) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp

  as.numeric(auc.roc(roc, partial.auc=attr(roc$auc, "partial.auc"), partial.auc.focus=attr(roc$auc, "partial.auc.focus"), partial.auc.correct=attr(roc$auc, "partial.auc.correct")))
}

nonstratified.ci.auc <- function(roc) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(controls, cases))

  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp
  
  as.numeric(auc.roc(roc, partial.auc=attr(roc$auc, "partial.auc"), partial.auc.focus=attr(roc$auc, "partial.auc.focus"), partial.auc.correct=attr(roc$auc, "partial.auc.correct")))
}

##########  AUC of a smooth ROC curve (ci.smooth.auc)  ##########

# Returns a smoothed auc in a stratified manner
stratified.ci.smooth.auc <- function(roc, smooth.roc.call, auc.call) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  # need to rebuild a ROC and smooth it
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- c(controls, cases)
  roc$response <- c(rep(roc$levels[1], length(controls)), rep(roc$levels[2], length(cases)))
  roc$thresholds <- thresholds

  # call smooth.roc and auc.smooth.roc
  smooth.roc.call$roc <- roc
  auc.call$smooth.roc <- try(eval(smooth.roc.call), silent=TRUE)
  if (is(auc.call$smooth.roc, "try-error")) {
    return(NA)
  }
  return(as.numeric(eval(auc.call)))
}

# Returns a smoothed auc in a non stratified manner
nonstratified.ci.smooth.auc <- function(roc, smooth.roc.call, auc.call) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(controls, cases))

  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- predictor
  roc$response <- response
  roc$thresholds <- thresholds

  # call smooth.roc and auc.smooth.roc
  smooth.roc.call$roc <- roc
  auc.call$smooth.roc <- try(eval(smooth.roc.call), silent=TRUE)
  if (is(auc.call$smooth.roc, "try-error")) {
    return(NA)
  }
  return(as.numeric(eval(auc.call)))
}

##########  AUC of a multiclass ROC (ci.multiclass.auc)  ##########
  
ci.multiclass.auc.bootstrap <- function(roc, conf.level, boot.n, boot.stratified, progress, ...) {
  if(class(progress) != "list")
    progress <- roc.utils.get.progress.bar(progress, title="Multi-class AUC confidence interval", label="Bootstrap in progress...", ...)

  if (boot.stratified) {
    aucs <- unlist(rlply(boot.n, stratified.ci.multiclass.auc(roc), .progress=progress))
  }
  else {
    aucs <- unlist(rlply(boot.n, nonstratified.ci.multiclass.auc(roc), .progress=progress))
  }

  if (sum(is.na(aucs)) > 0) {
    warning("NA value(s) produced during bootstrap were ignored.")
    aucs <- aucs[!is.na(aucs)]
  }
  # TODO: Maybe apply a correction (it's in the Tibshirani?) What do Carpenter-Bithell say about that?
  # Prepare the return value
  return(quantile(aucs, c(0+(1-conf.level)/2, .5, 1-(1-conf.level)/2)))

}

# Returns an auc in a stratified manner
stratified.ci.multiclass.auc <- function(roc) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp

  auc.roc(roc, partial.auc=attr(roc$auc, "partial.auc"), partial.auc.focus=attr(roc$auc, "partial.auc.focus"), partial.auc.correct=attr(roc$auc, "partial.auc.correct"))
}


# Returns an auc in a non stratified manner
nonstratified.ci.multiclass.auc <- function(roc) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(controls, cases))

  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp
  
  auc.roc(roc, partial.auc=attr(roc$auc, "partial.auc"), partial.auc.focus=attr(roc$auc, "partial.auc.focus"), partial.auc.correct=attr(roc$auc, "partial.auc.correct"))
}

##########  SE of a ROC curve (ci.se)  ##########

stratified.ci.se <- function(roc, sp) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$thresholds <- thresholds

  return(sapply(sp, function(x) coords.roc(roc, x, input="specificity", ret="sensitivity")))
}

nonstratified.ci.se <- function(roc, sp) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$thresholds <- thresholds

  return(sapply(sp, function(x) coords.roc(roc, x, input="specificity", ret="sensitivity")))
}

##########  SE of a smooth ROC curve (ci.se)  ##########

stratified.ci.smooth.se <- function(roc, sp, smooth.roc.call) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
    perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- c(controls, cases)
  roc$response <- c(rep(roc$levels[1], length(controls)), rep(roc$levels[2], length(cases)))
  roc$thresholds <- thresholds

  # call smooth.roc and auc.smooth.roc
  smooth.roc.call$roc <- roc
  smooth.roc <- try(eval(smooth.roc.call), silent=TRUE)
  if (is(smooth.roc, "try-error"))
    return(NA)
  return(sapply(sp, function(x) coords.smooth.roc(smooth.roc, x, input="specificity", ret="sensitivity")))
}

nonstratified.ci.smooth.se <- function(roc, sp, smooth.roc.call) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
    perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- predictor
  roc$response <- response
  roc$thresholds <- thresholds

  # call smooth.roc and auc.smooth.roc
  smooth.roc.call$roc <- roc
  smooth.roc <- try(eval(smooth.roc.call), silent=TRUE)
  if (is(smooth.roc, "try-error"))
    return(NA)
  return(sapply(sp, function(x) coords.smooth.roc(smooth.roc, x, input="specificity", ret="sensitivity")))
}

##########  SP of a ROC curve (ci.sp)  ##########

stratified.ci.sp <- function(roc, se) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$thresholds <- thresholds

  return(sapply(se, function(x) coords.roc(roc, x, input="sensitivity", ret="specificity")))
}

nonstratified.ci.sp <- function(roc, se) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$thresholds <- thresholds

  return(sapply(se, function(x) coords.roc(roc, x, input="sensitivity", ret="specificity")))
}

##########  SP of a smooth ROC curve (ci.sp)  ##########

stratified.ci.smooth.sp <- function(roc, se, smooth.roc.call) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
    perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- c(controls, cases)
  roc$response <- c(rep(roc$levels[1], length(controls)), rep(roc$levels[2], length(cases)))
  roc$thresholds <- thresholds

  # call smooth.roc and auc.smooth.roc
  smooth.roc.call$roc <- roc
  smooth.roc <- try(eval(smooth.roc.call), silent=TRUE)
  if (is(smooth.roc, "try-error"))
    return(NA)
  return(sapply(se, function(x) coords.smooth.roc(smooth.roc, x, input="sensitivity", ret="specificity")))
}

nonstratified.ci.smooth.sp <- function(roc, se, smooth.roc.call) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
    perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- predictor
  roc$response <- response
  roc$thresholds <- thresholds

  # call smooth.roc and auc.smooth.roc
  smooth.roc.call$roc <- roc
  smooth.roc <- try(eval(smooth.roc.call), silent=TRUE)
  if (is(smooth.roc, "try-error"))
    return(NA)
  return(sapply(se, function(x) coords.smooth.roc(smooth.roc, x, input="sensitivity", ret="specificity")))
}

##########  Threshold of a ROC curve (ci.thresholds)  ##########

stratified.ci.thresholds <- function(roc, thresholds) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  
  return(sapply(thresholds, roc.utils.perfs, controls=controls, cases=cases, direction=roc$direction))
}

# Returns an auc in a non stratified manner
nonstratified.ci.thresholds <- function(roc, thresholds) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]

  return(sapply(thresholds, roc.utils.perfs, controls=controls, cases=cases, direction=roc$direction))
}


##########  Coords of one ROC curves (ci.coords)  ##########
stratified.ci.coords <- function(roc, x, input, ret, best.method, best.weights) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- c(controls, cases)
  roc$response <- c(rep(roc$levels[1], length(controls)), rep(roc$levels[2], length(cases)))
  roc$thresholds <- thresholds

  as.numeric(coords.roc(roc, x=x, input=input, ret=ret, best.method=best.method, best.weights=best.weights))
}

nonstratified.ci.coords <- function(roc, x, input, ret, best.method, best.weights) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(controls, cases))


  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se
  roc$specificities <- perfs$sp
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- c(controls, cases)
  roc$response <- c(rep(roc$levels[1], length(controls)), rep(roc$levels[2], length(cases)))
  roc$thresholds <- thresholds
  
  as.numeric(coords.roc(roc, x=x, input=input, ret=ret, best.method=best.method, best.weights=best.weights))
}

##########  Coords of a smooth ROC curve (ci.coords)  ##########

stratified.ci.smooth.coords <- function(roc, x, input, ret, best.method, best.weights, smooth.roc.call) {
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  

  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- c(controls, cases)
  roc$response <- c(rep(roc$levels[1], length(controls)), rep(roc$levels[2], length(cases)))
  roc$thresholds <- thresholds

  # call smooth.roc and auc.smooth.roc
  smooth.roc.call$roc <- roc
  smooth.roc <- try(eval(smooth.roc.call), silent=TRUE)
  if (is(smooth.roc, "try-error"))
    return(NA)
  as.numeric(coords.roc(smooth.roc, x=x, input=input, ret=ret, best.method=best.method, best.weights=best.weights))
}

nonstratified.ci.smooth.coords <- function(roc, x, input, ret, best.method, best.weights, smooth.roc.call) {
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- roc$fun.sesp(thresholds=thresholds, controls=controls, cases=cases, direction=roc$direction)

  # update ROC
  roc$sensitivities <- perfs$se * ifelse(roc$percent, 100, 1)
  roc$specificities <- perfs$sp * ifelse(roc$percent, 100, 1)
  roc$cases <- cases
  roc$controls <- controls
  roc$predictor <- predictor
  roc$response <- response
  roc$thresholds <- thresholds

  # call smooth.roc and auc.smooth.roc
  smooth.roc.call$roc <- roc
  smooth.roc <- try(eval(smooth.roc.call), silent=TRUE)
  if (is(smooth.roc, "try-error"))
    return(NA)
  as.numeric(coords.roc(smooth.roc, x=x, input=input, ret=ret, best.method=best.method, best.weights=best.weights))
}
