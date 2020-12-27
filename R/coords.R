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

coords.smooth.roc <- function(smooth.roc,
                              x,
                              input,
                              ret=c("specificity", "sensitivity"),
                              as.list=FALSE,
                              drop=TRUE,
                              best.method=c("youden", "closest.topleft"),
                              best.weights=c(1, 0.5),
                              transpose = FALSE,
                              as.matrix = FALSE,
                              ...) {
  # make sure x was provided
  if (missing(x) || is.null(x) || (length(x) == 0 && !is.numeric(x))) {
    x <- "all"
  }
  else if (length(x) == 0 && is.numeric(x)) {
    stop("Numeric 'x' has length 0")
  }
  
  # match return
  ret <- roc.utils.match.coords.ret.args(ret, threshold = FALSE)

  if (is.character(x)) {
    x <- match.arg(x, c("all", "best")) # no thresholds in smoothed roc: only best or all are possible
    partial.auc <- attr(smooth.roc$auc, "partial.auc")
    
    if (x == "all") {
      # Pre-filter thresholds based on partial.auc
      if (is.null(smooth.roc$auc) || identical(partial.auc, FALSE)) {
        se <- smooth.roc$sensitivities
        sp <- smooth.roc$specificities
      }
      else {
        if (attr(smooth.roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- smooth.roc$sensitivities[smooth.roc$sensitivities <= partial.auc[1] & smooth.roc$sensitivities >= partial.auc[2]]
          sp <- smooth.roc$specificities[smooth.roc$sensitivities <= partial.auc[1] & smooth.roc$sensitivities >= partial.auc[2]]
        }
        else {
          se <- smooth.roc$sensitivities[smooth.roc$specificities <= partial.auc[1] & smooth.roc$specificities >= partial.auc[2]]
          sp <- smooth.roc$specificities[smooth.roc$specificities <= partial.auc[1] & smooth.roc$specificities >= partial.auc[2]]
        }
      }
      if (length(se) == 0) {
        warning("No coordinates found, returning NULL. This is possibly cased by a too small partial AUC interval.")
        return(NULL)
      }
      
      if (any(! ret %in% c("specificity", "sensitivity"))) {
        # Deduce additional tn, tp, fn, fp, npv, ppv
        res <- roc.utils.calc.coords(smooth.roc, NA, se, sp, best.weights)
      }
      else {
        res <- cbind(
          specificity = sp,
          sensitivity = se
        )
      }
    }
    else {
      # cheat: allow the user to pass "topleft"
      best.method <- match.arg(best.method[1], c("youden", "closest.topleft", "topleft"))
      if (best.method == "topleft") {
        best.method <- "closest.topleft"
      }
      optim.crit <- roc.utils.optim.crit(smooth.roc$sensitivities, smooth.roc$specificities,
                                         ifelse(smooth.roc$percent, 100, 1),
                                         best.weights, best.method)
      
      if (is.null(smooth.roc$auc) || identical(partial.auc, FALSE)) {
        se <- smooth.roc$sensitivities[optim.crit==max(optim.crit)]
        sp <- smooth.roc$specificities[optim.crit==max(optim.crit)]
        optim.crit <- optim.crit[optim.crit==max(optim.crit)]
      }
      else {
        if (attr(smooth.roc$auc, "partial.auc.focus") == "sensitivity") {
          optim.crit.partial <- (optim.crit)[smooth.roc$sensitivities <= partial.auc[1] & smooth.roc$sensitivities >= partial.auc[2]]
          se <- smooth.roc$sensitivities[smooth.roc$sensitivities <= partial.auc[1] & smooth.roc$sensitivities >= partial.auc[2]][optim.crit.partial==max(optim.crit.partial)]
          sp <- smooth.roc$specificities[smooth.roc$sensitivities <= partial.auc[1] & smooth.roc$sensitivities >= partial.auc[2]][optim.crit.partial==max(optim.crit.partial)]
          optim.crit <- optim.crit[smooth.roc$sensitivities <= partial.auc[1] & smooth.roc$sensitivities >= partial.auc[2]][optim.crit.partial==max(optim.crit.partial)]
        }
        else {
          optim.crit.partial <- (optim.crit)[smooth.roc$specificities <= partial.auc[1] & smooth.roc$specificities >= partial.auc[2]]
          se <- smooth.roc$sensitivities[smooth.roc$specificities <= partial.auc[1] & smooth.roc$specificities >= partial.auc[2]][optim.crit.partial==max(optim.crit.partial)]
          sp <- smooth.roc$specificities[smooth.roc$specificities <= partial.auc[1] & smooth.roc$specificities >= partial.auc[2]][optim.crit.partial==max(optim.crit.partial)]
          optim.crit <- optim.crit[smooth.roc$specificities <= partial.auc[1] & smooth.roc$specificities >= partial.auc[2]][optim.crit.partial==max(optim.crit.partial)]
        }
      }
      
      if (any(! ret %in% c("specificity", "sensitivity", best.method))) {
      	# Deduce additional tn, tp, fn, fp, npv, ppv
      	res <- roc.utils.calc.coords(smooth.roc, NA, se, sp, best.weights)
      }
      else {
      	res <- cbind(
      		specificity = sp,
      		sensitivity = se,
      		best.method = ifelse(best.method == "youden", 1, -1) * optim.crit
      	)
      	colnames(res)[3] <- best.method
      }
    }

    if (as.list) {
      warning("'as.list' is deprecated and will be removed in a future version.")
    	list <- apply(t(res)[ret, , drop=FALSE], 2, as.list)
    	if (drop == TRUE && length(x) == 1) {
    		return(list[[1]])
    	}
    	return(list)
    }
    else if (transpose) {
      rownames(res) <- NULL
      return(t(res)[ret,, drop=drop])
    }
    else {
      if (missing(drop) ) {
        drop = FALSE
      }
      if (! as.matrix) {
        res <- as.data.frame(res)
      }
      return(res[, ret, drop=drop])
    }
  }
  
  # Adjust drop for downstream call
  if (missing(drop) && ! transpose) {
    drop = FALSE
  }
  
  # match input 
  input <- roc.utils.match.coords.input.args(input, threshold = FALSE)
  
  # use coords.roc
  smooth.roc$thresholds <- rep(NA, length(smooth.roc$specificities))
  return(coords.roc(smooth.roc, x, input, ret, as.list, drop, 
                    transpose = transpose, as.matrix = as.matrix, ...))
}

coords.roc <- function(roc,
                       x,
                       input="threshold",
                       ret=c("threshold", "specificity", "sensitivity"),
                       as.list=FALSE,
                       drop=TRUE,
                       best.method=c("youden", "closest.topleft"),
                       best.weights=c(1, 0.5), 
                       transpose = FALSE,
                       as.matrix = FALSE,
                       ...) {
  # make sure x was provided
  if (missing(x) || is.null(x) || (length(x) == 0 && !is.numeric(x))) {
    x <- "all"
  }
  else if (length(x) == 0 && is.numeric(x)) {
    stop("Numeric 'x' has length 0")
  }
  
  # match input 
  input <- roc.utils.match.coords.input.args(input)
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
        se <- roc$sensitivities
        sp <- roc$specificities
        thres <- roc$thresholds
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- roc$sensitivities[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]]
          sp <- roc$specificities[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]]
          thres <- roc$thresholds[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]]
        }
        else {
          se <- roc$sensitivities[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]]
          sp <- roc$specificities[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]]
          thres <- roc$thresholds[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]]
        }
      }
      if (length(thres) == 0) {
      	warning("No coordinates found, returning NULL. This is possibly cased by a too small partial AUC interval.")
      	return(NULL)
      }
      res <- cbind(
      	threshold = thres,
      	specificity = sp,
      	sensitivity = se
      )
    }
    else if (x == "local maximas") {
      # Pre-filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        se <- roc$sensitivities
        sp <- roc$specificities
        thres <- roc$thresholds
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- roc$sensitivities[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]]
          sp <- roc$specificities[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]]
          thres <- roc$thresholds[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]]
        }
        else {
          se <- roc$sensitivities[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]]
          sp <- roc$specificities[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]]
          thres <- roc$thresholds[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]]
        }
      }
      if (length(thres) == 0) {
      	warning("No coordinates found, returning NULL. This is possibly cased by a too small partial AUC interval.")
      	return(NULL)
      }
      lm.idx <- roc.utils.max.thresholds.idx(thres, sp=sp, se=se)
      res <- cbind(
      	threshold = thres[lm.idx],
      	specificity = sp[lm.idx],
      	sensitivity = se[lm.idx]
      )
    }
    else { # x == "best"
      # cheat: allow the user to pass "topleft"
      best.method <- match.arg(best.method[1], c("youden", "closest.topleft", "topleft"))
      if (best.method == "topleft") {
        best.method <- "closest.topleft"
      }
      optim.crit <- roc.utils.optim.crit(roc$sensitivities, roc$specificities,
      								   ifelse(roc$percent, 100, 1),
      								   best.weights, best.method)

      # Filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
      	se <- roc$sensitivities[optim.crit==max(optim.crit)]
      	sp <- roc$specificities[optim.crit==max(optim.crit)]
        thres <- roc$thresholds[optim.crit==max(optim.crit)]
        optim.crit <- optim.crit[optim.crit==max(optim.crit)]
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          optim.crit <- (optim.crit)[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]]
          se <- roc$sensitivities[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]][optim.crit==max(optim.crit)]
          sp <- roc$specificities[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]][optim.crit==max(optim.crit)]
          thres <- roc$thresholds[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]][optim.crit==max(optim.crit)]
          optim.crit <- optim.crit[roc$sensitivities <= partial.auc[1] & roc$sensitivities >= partial.auc[2]][optim.crit==max(optim.crit)]
        }
        else {
          optim.crit <- (optim.crit)[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]]
          se <- roc$sensitivities[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]][optim.crit==max(optim.crit)]
          sp <- roc$specificities[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]][optim.crit==max(optim.crit)]
          thres <- roc$thresholds[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]][optim.crit==max(optim.crit)]
          optim.crit <- optim.crit[roc$specificities <= partial.auc[1] & roc$specificities >= partial.auc[2]][optim.crit==max(optim.crit)]
        }
      }
      if (length(thres) == 0) {
      	warning("No coordinates found, returning NULL. This is possibly cased by a too small partial AUC interval.")
      	return(NULL)
      }
      res <- cbind(
      	threshold = thres,
      	specificity = sp,
      	sensitivity = se,
      	best.method = ifelse(best.method == "youden", 1, -1) * optim.crit
      )
      colnames(res)[4] <- best.method
    }
  }
  else if (is.numeric(x)) {
    if (input == "threshold") {
    	thr_idx <- roc.utils.thr.idx(roc, x)
    	res <- cbind(
    		threshold = x, # roc$thresholds[thr_idx], # user-supplied vs ours.
    		specificity = roc$specificities[thr_idx],
    		sensitivity = roc$sensitivities[thr_idx]
    	)
    }
    else {
      # Arbitrary coord given in input.
      # We could be tempted to use all_coords directly.
      # However any non monotone coordinate in ret will be inaccurate
      # when interpolated. Therefore it is safer to only interpolate
      # se and sp and re-calculate the remaining coords later.
      res <- cbind(threshold = rep(NA, length(x)),
                        sensitivity = rep(NA, length(x)),
                        specificity = rep(NA, length(x))
                   )
      if (input %in% c("sensitivity", "specificity")) {
        # Shortcut slow roc.utils.calc.coords
        se <- roc$sensitivities
        sp <- roc$specificities
        if (methods::is(roc, "smooth.roc")) {
          thr <- rep(NA, length(roc$sensitivities))
        }
        else {
          thr <- roc$thresholds
        }
        if (input == "sensitivity") {
          input_values <- se
        }
        else {
          input_values <- sp
        }
      }
      else {
        all_coords <- roc.utils.calc.coords(roc, rep(NA, length(roc$sensitivities)), roc$sensitivities, roc$specificities, best.weights)
        input_values <- all_coords[, input]
        se <- all_coords[, "sensitivity"]
        sp <- all_coords[, "specificity"]
        thr <- all_coords[, "threshold"]
      }
      for (i in seq_along(x)) {
        value <- x[i]
        if (value < min(input_values) || value > max(input_values)) {
          stop(sprintf("Input %s (%s) not in range (%s-%s)", input, value, 
                       min(input_values), max(input_values)))
        }
        
        idx <- which(input_values == value)
        if (length(idx) > 1) {
          # More than one to pick from. Need to take best
          # according to sorting
          if (coord.is.decreasing[input]) {
            idx <- idx[length(idx)] # last
          }
          else {
            idx <- idx[1] # first
          }
        }
        if (length(idx) == 1) {
          # Exactly one to pick from
          res[i,] <- c(thr[idx], se[idx], sp[idx])
        }
        else {
          # Need to interpolate
          if (coord.is.decreasing[input]) {
            idx.next <- match(TRUE, input_values < value)
          }
          else {
            idx.next <- match(TRUE, input_values > value)
          }
          proportion <- (value - input_values[idx.next]) / (input_values[idx.next - 1] - input_values[idx.next])
          int.se <- se[idx.next] + proportion * (se[idx.next - 1] - se[idx.next])
          int.sp <- sp[idx.next] + proportion * (sp[idx.next - 1] - sp[idx.next])
          res[i, 2:3] <- c(int.se, int.sp)
        }
      }
    }
  }
  else {
    stop("'x' must be a numeric or character vector.")
  }
  
  if (any(! ret %in% colnames(res))) {
  	# Deduce additional tn, tp, fn, fp, npv, ppv
  	res <- roc.utils.calc.coords(roc, res[, "threshold"], res[, "sensitivity"], res[, "specificity"], best.weights)
  }

  if (as.list) {
    warning("'as.list' is deprecated and will be removed in a future version.")
  	list <- apply(t(res)[ret, , drop=FALSE], 2, as.list)
  	if (drop == TRUE && length(x) == 1) {
  		return(list[[1]])
  	}
  	return(list)
  }
  else if (transpose) {
    rownames(res) <- NULL
    return(t(res)[ret,, drop=drop])
  }
  else {
  	# HACK:
  	# We need an exception for r4lineups that will keep the old drop = TRUE 
  	# behavior, until r4lineups gets updated. This is an ugly hack but allows
  	# us to switch to a better drop = FALSE for everyone else
    if (missing(drop)) {
    	if (sys.nframe() > 2 && 
    	    length(deparse(sys.call(-2))) == 1 && 
    	    deparse(sys.call(-2)) == 'make_rocdata(df_confacc)' && 
    	    length(deparse(sys.call(-1))) == 1 && (
    	      deparse(sys.call(-1)) == 'coords(rocobj, "all", ret = c("tp"))' ||
    	      deparse(sys.call(-1)) == 'coords(rocobj, "all", ret = "fp")'
    	    )) {
    		# We're in r4lineups
    		drop = TRUE
    	}
    	else {
    		drop = FALSE
    	}
    }
    if (! as.matrix) {
      res <- as.data.frame(res)
    }
    return(res[, ret, drop=drop])
  	
  }
}
