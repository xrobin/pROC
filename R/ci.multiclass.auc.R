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

ci.multiclass.auc <- function(multiclass.auc, ...)
  ci.multiclass.roc(attr(multiclass.auc, "roc"), ...)

ci.multiclass.roc <- function(multiclass.roc,
                   conf.level = 0.95,
                   boot.n = 2000,
                   boot.stratified = TRUE,
                   reuse.auc=TRUE,
                   progress = getOption("pROCProgress")$name,
                   ...
                   ) {
  if (conf.level > 1 | conf.level < 0)
    stop("conf.level must be within the interval [0,1].")

  # We need an auc
  if (is.null(multiclass.roc$auc) | !reuse.auc)
    multiclass.roc$auc <- auc(multiclass.roc, ...)

  # do all the computations in fraction, re-transform in percent later if necessary
  percent <- multiclass.roc$percent
  multiclass.roc$percent <- FALSE
  oldauc <- multiclass.roc$auc
  if (percent) {
    attr(multiclass.roc$auc, "percent") <- FALSE
    if (! identical(attr(multiclass.roc$auc, "partial.auc"), FALSE))
      attr(multiclass.roc$auc, "partial.auc") <- attr(multiclass.roc$auc, "partial.auc") / 100
  }

  ci <- ci.multiclass.auc.bootstrap(multiclass.roc, conf.level, boot.n, boot.stratified, progress, ...)

  if (percent) {
    ci <- ci * 100
  }
  attr(ci, "conf.level") <- conf.level
  attr(ci, "boot.n") <- boot.n
  attr(ci, "boot.stratified") <- boot.stratified
  attr(ci, "multiclass.auc") <- oldauc
  class(ci) <- "ci.multiclass.auc"
  return(ci)
}
  
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
  browser()
  controls <- sample(roc$controls, replace=TRUE)
  cases <- sample(roc$cases, replace=TRUE)
  thresholds <- roc.utils.thresholds(c(cases, controls))
  
  perfs <- sapply(thresholds, roc.utils.perfs, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs[2,]
  roc$specificities <- perfs[1,]

  auc.roc(roc, partial.auc=attr(roc$auc, "partial.auc"), partial.auc.focus=attr(roc$auc, "partial.auc.focus"), partial.auc.correct=attr(roc$auc, "partial.auc.correct"))
}


# Returns an auc in a non stratified manner
nonstratified.ci.multiclass.auc <- function(roc) {
  browser()
  tmp.idx <- sample(1:length(roc$predictor), replace=TRUE)
  predictor <- roc$predictor[tmp.idx]
  response <- roc$response[tmp.idx]
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(roc$levels[1])]]
  cases <- splitted[[as.character(roc$levels[2])]]
  thresholds <- roc.utils.thresholds(c(controls, cases))

  perfs <- sapply(thresholds, roc.utils.perfs, controls=controls, cases=cases, direction=roc$direction)
  roc$sensitivities <- perfs[2,]
  roc$specificities <- perfs[1,]
  
  auc.roc(roc, partial.auc=attr(roc$auc, "partial.auc"), partial.auc.focus=attr(roc$auc, "partial.auc.focus"), partial.auc.correct=attr(roc$auc, "partial.auc.correct"))
}
