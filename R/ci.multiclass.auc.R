# # pROC: Tools Receiver operating characteristic (ROC curves) with
# # (partial) area under the curve, confidence intervals and comparison. 
# # Copyright (C) 2010-2014 Xavier Robin, Alexandre Hainard, Natacha Turck,
# # Natalia Tiberti, Frédérique Lisacek, Jean-Charles Sanchez
# # and Markus Müller
# #
# # This program is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.
# #
# # This program is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.
# #
# # You should have received a copy of the GNU General Public License
# # along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# ci.multiclass.auc <- function(multiclass.auc, ...) {
# 	stop("CI of a multiclass AUC not implemented")
# 	ci.multiclass.roc(attr(multiclass.auc, "roc"), ...)
# }
# 
# ci.multiclass.roc <- function(multiclass.roc,
#                    conf.level = 0.95,
#                    boot.n = 2000,
#                    boot.stratified = TRUE,
#                    reuse.auc=TRUE,
#                    progress = getOption("pROCProgress")$name,
#                    parallel = FALSE,
#                    ...
#                    ) {
# 	stop("ci of a multiclass ROC curve not implemented")
#   if (conf.level > 1 | conf.level < 0)
#     stop("conf.level must be within the interval [0,1].")
# 
#   # We need an auc
#   if (is.null(multiclass.roc$auc) | !reuse.auc)
#     multiclass.roc$auc <- auc(multiclass.roc, ...)
# 
#   # do all the computations in fraction, re-transform in percent later if necessary
#   percent <- multiclass.roc$percent
#   oldauc <- multiclass.roc$auc
#   if (percent) {
#   	multiclass.roc <- roc.utils.unpercent(multiclass.roc)
#   }
# 
#   ci <- ci.multiclass.auc.bootstrap(multiclass.roc, conf.level, boot.n, boot.stratified, progress, parallel, ...)
# 
#   if (percent) {
#     ci <- ci * 100
#   }
#   attr(ci, "conf.level") <- conf.level
#   attr(ci, "boot.n") <- boot.n
#   attr(ci, "boot.stratified") <- boot.stratified
#   attr(ci, "multiclass.auc") <- oldauc
#   class(ci) <- "ci.multiclass.auc"
#   return(ci)
# }
