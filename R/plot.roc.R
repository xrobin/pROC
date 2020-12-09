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

plot.roc <- function(x, ...) {
  UseMethod("plot.roc")
}

plot.roc.formula <- function(x, data, subset, na.action, ...) {
	data.missing <- missing(data)
	call <- match.call()
	names(call)[2] <- "formula" # forced to be x by definition of plot
	roc.data <- roc.utils.extract.formula(formula=x, data, subset, na.action, ..., 
										  data.missing = data.missing,
										  call = call)
	if (length(roc.data$predictor.name) > 1) {
		stop("Only one predictor supported in 'plot.roc'.")
	}
	response <- roc.data$response
	predictor <- roc.data$predictors[, 1]
	
	roc <- roc(response, predictor, plot=TRUE, ...)
	roc$call <- match.call()
	invisible(roc)
}

plot.roc.default <- function(x, predictor, ...) {
  roc <- roc(x, predictor, plot=TRUE, ...)
  roc$call <- match.call()
  invisible(roc)
}

plot.roc.smooth.roc <- plot.smooth.roc <- function(x, ...) {
  invisible(plot.roc.roc(x, ...)) # force usage of plot.roc.roc: only print.thres not working
}

plot.roc.roc <- function(x,
                         add=FALSE,
                         reuse.auc=TRUE,
                         axes=TRUE,
                         legacy.axes=FALSE,
                         xlim=if(x$percent){c(100, 0)} else{c(1, 0)},
                         ylim=if(x$percent){c(0, 100)} else{c(0, 1)},
                         xlab=ifelse(x$percent, ifelse(legacy.axes, "100 - Specificity (%)", "Specificity (%)"), ifelse(legacy.axes, "1 - Specificity", "Specificity")),
                         ylab=ifelse(x$percent, "Sensitivity (%)", "Sensitivity"),
                         asp=1,
                         mar=c(4, 4, 2, 2)+.1,
                         mgp=c(2.5, 1, 0),
                         # col, lty and lwd for the ROC line only
                         col=par("col"),
                         lty=par("lty"),
                         lwd=2,
                         type="l",
                         # Identity line
                         identity=!add,
                         identity.col="darkgrey",
                         identity.lty=1,
                         identity.lwd=1,
                         # Print the thresholds on the plot
                         print.thres=FALSE,
                         print.thres.pch=20,
                         print.thres.adj=c(-.05,1.25),
                         print.thres.col="black",
                         print.thres.pattern=ifelse(x$percent, "%.1f (%.1f%%, %.1f%%)", "%.3f (%.3f, %.3f)"),
                         print.thres.cex=par("cex"),
                         print.thres.pattern.cex=print.thres.cex,
                         print.thres.best.method=NULL,
                         print.thres.best.weights=c(1, 0.5),
                         # Print the AUC on the plot
                         print.auc=FALSE,
                         print.auc.pattern=NULL,
                         print.auc.x=ifelse(x$percent, 50, .5), 
                         print.auc.y=ifelse(x$percent, 50, .5),
                         print.auc.adj=c(0,1),
                         print.auc.col=col,
                         print.auc.cex=par("cex"),
                         # Grid
                         grid=FALSE,
                         grid.v={
                           if(is.logical(grid) && grid[1]==TRUE){seq(0, 1, 0.1) * ifelse(x$percent, 100, 1)}
                           else if(is.numeric(grid)) {seq(0, ifelse(x$percent, 100, 1), grid[1])}
                           else {NULL}
                         },
                         grid.h={
                           if (length(grid) == 1) {grid.v}
                           else if (is.logical(grid) && grid[2]==TRUE){seq(0, 1, 0.1) * ifelse(x$percent, 100, 1)}
                           else if(is.numeric(grid)) {seq(0, ifelse(x$percent, 100, 1), grid[2])}
                           else {NULL}
                         },
                         # for grid.lty, grid.lwd and grid.col, a length 2 value specifies both values for vertical (1) and horizontal (2) grid
                         grid.lty=3,
                         grid.lwd=1,
                         grid.col="#DDDDDD",
                         # Polygon for the auc
                         auc.polygon=FALSE,
                         auc.polygon.col="gainsboro", # Other arguments can be passed to polygon() using "..." (for these two we cannot)
                         auc.polygon.lty=par("lty"),
                         auc.polygon.density=NULL,
                         auc.polygon.angle=45,
                         auc.polygon.border=NULL,
                         # Should we show the maximum possible area as another polygon?
                         max.auc.polygon=FALSE,
                         max.auc.polygon.col="#EEEEEE", # Other arguments can be passed to polygon() using "..." (for these two we cannot)
                         max.auc.polygon.lty=par("lty"),
                         max.auc.polygon.density=NULL,
                         max.auc.polygon.angle=45,
                         max.auc.polygon.border=NULL,
                         # Confidence interval
                         ci=!is.null(x$ci),
                         ci.type=c("bars", "shape", "no"),
                         ci.col=ifelse(ci.type=="bars", par("fg"), "gainsboro"),
                         ...
                         ) {
  percent <- x$percent
  
  if (max.auc.polygon | auc.polygon | print.auc) {# we need the auc here
    if (is.null(x$auc) | !reuse.auc)
      x$auc <- auc(x, ...)
    partial.auc <- attr(x$auc, "partial.auc")
    partial.auc.focus <- attr(x$auc, "partial.auc.focus")
  }

  # compute a reasonable default for print.auc.pattern if required
  if (print.auc & is.null(print.auc.pattern)) {
    print.auc.pattern <- ifelse(identical(partial.auc, FALSE), "AUC: ", "Partial AUC: ")
    print.auc.pattern <- paste(print.auc.pattern, ifelse(percent, "%.1f%%", "%.3f"), sep="")
    if (ci && methods::is(x$ci, "ci.auc"))
      print.auc.pattern <- paste(print.auc.pattern, " (", ifelse(percent, "%.1f%%", "%.3f"), "\u2013", ifelse(percent, "%.1f%%", "%.3f"), ")",sep="")
  }
    
  # get and sort the sensitivities and specificities
  se <- sort(x$sensitivities, decreasing=TRUE)
  sp <- sort(x$specificities, decreasing=FALSE) 
  if (!add) {
    opar <- par(mar=mar, mgp=mgp)
    on.exit(par(opar))
    # type="n" to plot background lines and polygon shapes first. We will add the line later. axes=FALSE, we'll add them later according to legacy.axis
    suppressWarnings(plot(x$specificities, x$sensitivities, xlab=xlab, ylab=ylab, type="n", axes=FALSE, xlim=xlim, ylim=ylim, lwd=lwd, asp=asp, ...))

    # As we had axes=FALSE we need to add them again unless axes=FALSE
    if (axes) {
      box()
      # axis behave differently when at and labels are passed (no decimals on 1 and 0), 
      # so handle each case separately and consistently across axes
      if (legacy.axes) {
      	lab.at <- axTicks(side=1)
      	lab.labels <- format(ifelse(x$percent, 100, 1) - lab.at)
      	suppressWarnings(axis(side=1, at=lab.at, labels=lab.labels, ...))
      	lab.at <- axTicks(side=2)
      	suppressWarnings(axis(side=2, at=lab.at, labels=format(lab.at), ...))
      } 
      else {
      	suppressWarnings(axis(side=1, ...))
      	suppressWarnings(axis(side=2, ...))
      }
    }
  }

  # Plot the grid
  # make sure grid.lty, grid.lwd and grid.col are at least of length 2
  grid.lty <- rep(grid.lty, length.out=2)
  grid.lwd <- rep(grid.lwd, length.out=2)
  grid.col <- rep(grid.col, length.out=2)
  if (!is.null(grid.v)) {
    suppressWarnings(abline(v=grid.v, lty=grid.lty[1], col=grid.col[1], lwd=grid.lwd[1], ...))
  }
  if (!is.null(grid.h)) {
    suppressWarnings(abline(h=grid.h, lty=grid.lty[2], col=grid.col[2], lwd=grid.lwd[2], ...))
  }

  # Plot the polygon displaying the maximal area
  if (max.auc.polygon) {
    if (identical(partial.auc, FALSE)) {
      map.y <- c(0, 1, 1, 0) * ifelse(percent, 100, 1)
      map.x <- c(1, 1, 0, 0) * ifelse(percent, 100, 1)
    }
    else {
      if (partial.auc.focus=="sensitivity") {
        map.y <- c(partial.auc[2], partial.auc[2], partial.auc[1], partial.auc[1]) 
        map.x <- c(0, 1, 1, 0) * ifelse(percent, 100, 1) 
      }
      else {
        map.y <- c(0, 1, 1, 0) * ifelse(percent, 100, 1) 
        map.x <- c(partial.auc[2], partial.auc[2], partial.auc[1], partial.auc[1])
      }
    }
    suppressWarnings(polygon(map.x, map.y, col=max.auc.polygon.col, lty=max.auc.polygon.lty, border=max.auc.polygon.border, density=max.auc.polygon.density, angle=max.auc.polygon.angle, ...))
  }
  # Plot the ci shape
  if (ci && ! methods::is(x$ci, "ci.auc")) {
    ci.type <- match.arg(ci.type)
    if (ci.type=="shape")
      plot(x$ci, type="shape", col=ci.col, no.roc=TRUE, ...)
  }
  # Plot the polygon displaying the actual area
  if (auc.polygon) {
    if (identical(partial.auc, FALSE)) {
      suppressWarnings(polygon(c(sp, 0), c(se, 0), col=auc.polygon.col, lty=auc.polygon.lty, border=auc.polygon.border, density=auc.polygon.density, angle=auc.polygon.angle, ...))
    }
    else {
      if (partial.auc.focus == "sensitivity") {
        x.all <- rev(se)
        y.all <- rev(sp)
      }
      else {
        x.all <- sp
        y.all <- se
      }
      # find the SEs and SPs in the interval
      x.int <- x.all[x.all <= partial.auc[1] & x.all >= partial.auc[2]]
      y.int <- y.all[x.all <= partial.auc[1] & x.all >= partial.auc[2]]
      # if the upper limit is not exactly present in SPs, interpolate
      if (!(partial.auc[1] %in% x.int)) {
        x.int <- c(x.int, partial.auc[1])
        # find the limit indices
        idx.out <- match(FALSE, x.all < partial.auc[1])
        idx.in <- idx.out - 1
        # interpolate y
        proportion.start <- (partial.auc[1] - x.all[idx.out]) / (x.all[idx.in] - x.all[idx.out])
        y.start <- y.all[idx.out] - proportion.start * (y.all[idx.out] - y.all[idx.in])
        y.int <- c(y.int, y.start)
      }
      # if the lower limit is not exactly present in SPs, interpolate
      if (!(partial.auc[2] %in% x.int)) {
        x.int <- c(partial.auc[2], x.int)
        # find the limit indices
        idx.out <- length(x.all) - match(TRUE, rev(x.all) < partial.auc[2]) + 1
        idx.in <- idx.out + 1
        # interpolate y
        proportion.end <- (x.all[idx.in] - partial.auc[2]) / (x.all[idx.in] - x.all[idx.out])
        y.end <- y.all[idx.in] + proportion.end * (y.all[idx.out] - y.all[idx.in])
        y.int <- c(y.end, y.int)
      }
      # anchor to baseline
      x.int <- c(partial.auc[2], x.int, partial.auc[1])
      y.int <- c(0, y.int, 0)
      if (partial.auc.focus == "sensitivity") {
        # for SE, invert x and y again
        suppressWarnings(polygon(y.int, x.int, col=auc.polygon.col, lty=auc.polygon.lty, border=auc.polygon.border, density=auc.polygon.density, angle=auc.polygon.angle, ...))
      }
      else {
        suppressWarnings(polygon(x.int, y.int, col=auc.polygon.col, lty=auc.polygon.lty, border=auc.polygon.border, density=auc.polygon.density, angle=auc.polygon.angle, ...))
      }
    }
  }
  # Identity line
  if (identity) suppressWarnings(abline(ifelse(percent, 100, 1), -1, col=identity.col, lwd=identity.lwd, lty=identity.lty, ...))
  # Actually plot the ROC curve
  suppressWarnings(lines(sp, se, type=type, lwd=lwd, col=col, lty=lty, ...))
  # Plot the ci bars
  if (ci && !methods::is(x$ci, "ci.auc")) {
    if (ci.type=="bars")
      plot(x$ci, type="bars", col=ci.col, ...)
  }
  # Print the thresholds on the curve if print.thres is TRUE
  if (isTRUE(print.thres))
    print.thres <- "best"
  if (is.character(print.thres))
    print.thres <- match.arg(print.thres, c("no", "all", "local maximas", "best"))
  if (methods::is(x, "smooth.roc")) {
    if (is.numeric(print.thres))
      stop("Numeric 'print.thres' unsupported on a smoothed ROC plot.")
    else if (print.thres == "all" || print.thres == "local maximas")
      stop("'all' and 'local maximas' 'print.thres' unsupported on a smoothed ROC plot.") 
    else if (print.thres == "best") {
      co <- coords(x, print.thres, best.method=print.thres.best.method, best.weights=print.thres.best.weights, transpose=FALSE)
        suppressWarnings(points(co$specificity, co$sensitivity, pch=print.thres.pch, cex=print.thres.cex, col=print.thres.col, ...))
        suppressWarnings(text(co$specificity, co$sensitivity, sprintf(print.thres.pattern, NA, co$specificity, co$sensitivity), adj=print.thres.adj, cex=print.thres.pattern.cex, col=print.thres.col, ...))
     
    } # else print.thres == no > do nothing
  }
  else if (is.numeric(print.thres) || is.character(print.thres)) {
    if (is.character(print.thres) && print.thres == "no") {} # do nothing
    else {
      co <- coords(x, print.thres, best.method=print.thres.best.method, best.weights=print.thres.best.weights, transpose=FALSE)
      suppressWarnings(points(co$specificity, co$sensitivity, pch=print.thres.pch, cex=print.thres.cex, col=print.thres.col, ...))
      suppressWarnings(text(co$specificity, co$sensitivity, sprintf(print.thres.pattern, co$threshold, co$specificity, co$sensitivity), adj=print.thres.adj, cex=print.thres.pattern.cex, col=print.thres.col, ...))
    }
  }

  # Print the AUC on the plot
  if (print.auc) {
    if (ci && methods::is(x$ci, "ci.auc")) {
      labels <- sprintf(print.auc.pattern, x$auc, x$ci[1], x$ci[3])
      suppressWarnings(text(print.auc.x, print.auc.y, labels, adj=print.auc.adj, cex=print.auc.cex, col=print.auc.col, ...))
    }
    else
      labels <- sprintf(print.auc.pattern, x$auc)
    suppressWarnings(text(print.auc.x, print.auc.y, labels, adj=print.auc.adj, cex=print.auc.cex, col=print.auc.col, ...))
  }
  
  invisible(x)
}
