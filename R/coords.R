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

coords <- function(...)
  UseMethod("coords")

coords.smooth.roc <- function(smooth.roc, x, input=c("specificity", "sensitivity"), ret=c("specificity", "sensitivity"), as.list=FALSE, drop=TRUE, best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5), ...) {
  # make sure x was provided
  if (missing(x))
    stop("'x' must be a numeric or character vector.")
  # match input 
  input <- match.arg(input)
  # match return
  ret <- roc.utils.match.coords.ret.args(ret, threshold = FALSE)

  if (is.character(x)) {
    x <- match.arg(x, c("best")) # no thresholds in smoothed roc: only best is possible
    partial.auc <- attr(smooth.roc$auc, "partial.auc")
    # What kind of "best" do we want?
    # Compute weights
    if (! is.numeric(best.weights) || length(best.weights) != 2)
      stop("'best.weights' must be a numeric vector of length 2.")
    if (best.weights[2] <= 0 || best.weights[2] >= 1)
      stop("prevalence ('best.weights[2]') must be in the interval ]0,1[.")
    r <- (1 - best.weights[2]) / (best.weights[1] * best.weights[2]) # r should be 1 by default
      
    # Compute optimality criterion and store it in the optim.crit vector
    best.method <- match.arg(best.method[1], c("youden", "closest.topleft", "topleft")) # cheat: allow the user to pass "topleft"
    if (best.method == "youden") {
      optim.crit <- smooth.roc$sensitivities + r * smooth.roc$specificities
    }
    else if (best.method == "closest.topleft" || best.method == "topleft") {
      fac.1 <- ifelse(smooth.roc$percent, 100, 1)
      optim.crit <- - ((fac.1 - smooth.roc$sensitivities)^2 + r * (fac.1 - smooth.roc$specificities)^2)
    }
    
    if (is.null(smooth.roc$auc) || identical(partial.auc, FALSE)) {
      se <- smooth.roc$sensitivities[optim.crit==max(optim.crit)]
      sp <- smooth.roc$specificities[optim.crit==max(optim.crit)]
    }
    else {
      if (attr(smooth.roc$auc, "partial.auc.focus") == "sensitivity") {
        optim.crit <- (optim.crit)[smooth.roc$se <= partial.auc[1] & smooth.roc$se >= partial.auc[2]]
        se <- smooth.roc$sensitivities[smooth.roc$sensitivities <= partial.auc[1] & smooth.roc$sensitivities >= partial.auc[2]][optim.crit==max(optim.crit)]
        sp <- smooth.roc$specificities[smooth.roc$sensitivities <= partial.auc[1] & smooth.roc$sensitivities >= partial.auc[2]][optim.crit==max(optim.crit)]
      }
      else {
        optim.crit <- (optim.crit)[smooth.roc$sp <= partial.auc[1] & smooth.roc$sp >= partial.auc[2]]
        se <- smooth.roc$sensitivities[smooth.roc$specificities <= partial.auc[1] & smooth.roc$specificities >= partial.auc[2]][optim.crit==max(optim.crit)]
        sp <- smooth.roc$specificities[smooth.roc$specificities <= partial.auc[1] & smooth.roc$specificities >= partial.auc[2]][optim.crit==max(optim.crit)]
      }
    }
    # Deduce additional tn, tp, fn, fp, npv, ppv
    ncases <- length(attr(smooth.roc, "roc")$cases)
    ncontrols <- length(attr(smooth.roc, "roc")$controls)
    substr.percent <- ifelse(smooth.roc$percent, 100, 1)
    co <- roc.utils.calc.coords(substr.percent,
    							se, sp, ncases, ncontrols)
    
    if (length(se) == 1) {
      if (as.list) {
        list <- as.list(co)[ret]
        if (drop == FALSE) {
          list <- list(list)
          names(list) <- x
        }
        return(list)
      }
      else {
        if (drop == FALSE) {
        	return(t(co))
        }
      	else {
      		return(t(co)[1,])
      	}
      }
    }
    else if (length(se) > 1) {
      if (as.list) {
        co <- apply(t(co), 2, as.list)
        names(co) <- rep("best", length(co))
        return(co)
      }
      else {
        co <- t(co)
        colnames(co) <- rep(x, dim(co)[2])
        return(co)
      }
    }
  }

  # use coords.roc
  smooth.roc$thresholds <- rep(NA, length(smooth.roc$specificities))
  coords.roc(smooth.roc, x, input, ret, as.list, drop, ...)
}

coords.roc <- function(roc, x, input=c("threshold", "specificity", "sensitivity"), ret=c("threshold", "specificity", "sensitivity"), as.list=FALSE, drop=TRUE, best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5), ...) {
  # make sure x was provided
  if (missing(x) || length(x) == 0)
    stop("'x' must be a numeric or character vector of positive length.")
  # match input 
  input <- match.arg(input)
  # match return
  ret <- roc.utils.match.coords.ret.args(ret)
  # make sure the sort of roc is correct
  roc <- sort(roc)

  if (is.character(x)) {
    x <- match.arg(x, c("all", "local maximas", "best"))
    partial.auc <- attr(roc$auc, "partial.auc")
    if (x == "all") {
      # Pre-filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        se <- roc$se
        sp <- roc$sp
        thres <- roc$thresholds
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- roc$se[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          sp <- roc$sp[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
        }
        else {
          se <- roc$se[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          sp <- roc$sp[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      co <- coords(roc, x=thres, input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (class(co) == "matrix")
        colnames(co) <- rep(x, dim(co)[2])
      else if (class(co) == "list" && class(co[[1]]) == "list")
        names(co) <- rep(x, length(co))
      return(co)
    }
    else if (x == "local maximas") {
      # Pre-filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        se <- roc$se
        sp <- roc$sp
        thres <- roc$thresholds
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- roc$se[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          sp <- roc$sp[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
        }
        else {
          se <- roc$se[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          sp <- roc$sp[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      lm.idx <- roc.utils.max.thresholds.idx(thres, sp=sp, se=se)
      co <- coords(roc, x=thres[lm.idx], input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (class(co) == "matrix")
        colnames(co) <- rep(x, dim(co)[2])
      else if (class(co) == "list" && class(co[[1]]) == "list")
        names(co) <- rep(x, length(co))
      return(co)
    }
    else { # x == "best"
      # What kind of "best" do we want?
      # Compute weights
      if (is.numeric(best.weights) && length(best.weights) == 2)
        r <- (1 - best.weights[2]) / (best.weights[1] * best.weights[2]) # r should be 1 by default
      else
        stop("'best.weights' must be a numeric vector of length 2")
      # Compute optimality criterion and store it in the optim.crit vector
      best.method <- match.arg(best.method[1], c("youden", "closest.topleft", "topleft")) # cheat: allow the user to pass "topleft"
      if (best.method == "youden") {
        optim.crit <- roc$sensitivities + r * roc$specificities
      }
      else if (best.method == "closest.topleft" || best.method == "topleft") {
        fac.1 <- ifelse(roc$percent, 100, 1)
        optim.crit <- - ((fac.1 - roc$sensitivities)^2 + r * (fac.1 - roc$specificities)^2)
      }

      # Filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        thres <- roc$thresholds[optim.crit==max(optim.crit)]
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          optim.crit <- (optim.crit)[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]][optim.crit==max(optim.crit)]
        }
        else {
          optim.crit <- (optim.crit)[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]][optim.crit==max(optim.crit)]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      co <- coords(roc, x=thres, input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (class(co) == "matrix")
        colnames(co) <- rep(x, dim(co)[2])
      else if (class(co) == "list" && class(co[[1]]) == "list")
        names(co) <- rep(x, length(co))
      return(co)
    }
  }
  else if (is.numeric(x)) {
    if (input == "threshold") {
    	# We must match every threshold given by the user to one of our
    	# selected threshold. However we need to be careful to assign 
    	# them to the right one around the exact data point values
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
    	}
    	else {
    		cut_points <- c(rev(cut_points), Inf)
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
    	res <- rbind(
    		threshold = x, # roc$thresholds[thr_idx], # user-supplied vs ours.
    		specificity = roc$specificities[thr_idx],
    		sensitivity = roc$sensitivities[thr_idx]
    	)
    }
    if (input == "specificity") {
    	if (any(x < 0) || any(x > ifelse(roc$percent, 100, 1))) {
    		stop("Input specificity not within the ROC space.")
    	}
    	res <- matrix(nrow=3, ncol=length(x))
    	for (i in seq_along(x)) {
    		sp <- x[i]
    		if (sp %in% roc$sp) {
    			idx <- match(sp, roc$sp)
    			res[, i] <- c(roc$thresholds[idx], roc$sp[idx], roc$se[idx])
    		}
    		else { # need to interpolate
    			idx.next <- match(TRUE, roc$sp > sp)
    			proportion <-  (sp - roc$sp[idx.next - 1]) / (roc$sp[idx.next] - roc$sp[idx.next - 1])
    			int.se <- roc$se[idx.next - 1] - proportion * (roc$se[idx.next - 1] - roc$se[idx.next])
    			res[, i] <- c(NA, sp, int.se)
    		}
    	}
    }
    if (input == "sensitivity") {
    	if (any(x < 0) || any(x > ifelse(roc$percent, 100, 1))) {
    		stop("Input sensitivity not within the ROC space.")
    	}
    	res <- matrix(nrow=3, ncol=length(x))
    	for (i in seq_along(x)) {
    		se <- x[i]
    		if (se %in% roc$se) {
    			idx <- length(roc$se) + 1 - match(TRUE, rev(roc$se) == se)
    			res[, i] <- c(roc$thresholds[idx], roc$sp[idx], roc$se[idx])
    		}
    		else { # need to interpolate
    			idx.next <- match(TRUE, roc$se < se)
    			proportion <- (se - roc$se[idx.next]) / (roc$se[idx.next - 1] - roc$se[idx.next])
    			int.sp <- roc$sp[idx.next] + proportion * (roc$sp[idx.next - 1] - roc$sp[idx.next])
    			res[, i] <- c(NA, int.sp, se)
    		}
    	}
    }
    # Deduce additional tn, tp, fn, fp, npv, ppv
    ncases <- ifelse(methods::is(roc, "smooth.roc"), length(attr(roc, "roc")$cases), length(roc$cases))
    ncontrols <- ifelse(methods::is(roc, "smooth.roc"), length(attr(roc, "roc")$controls), length(roc$controls))
    se <- res[3, ]
    sp <- res[2, ]
    
    substr.percent <- ifelse(roc$percent, 100, 1)
    co <- roc.utils.calc.coords(substr.percent,
    							se, sp, ncases, ncontrols)
    co <- cbind(threshold = res[1, ], co)
    rownames(co) <- x

    if (as.list) {
    	list <- apply(co[, ret, drop=FALSE], 1, as.list)
    	if (drop == TRUE && length(x) == 1) {
    		return(list[[1]])
    	}
    	return(list)
    }
    else {
      res <- t(co)
      colnames(res) <- x
      return(res[ret,, drop=drop])
    }
  }
  else {
    stop("'x' must be a numeric or character vector.")
  }
}
