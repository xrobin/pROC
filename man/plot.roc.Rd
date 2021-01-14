\encoding{UTF-8}
\name{plot.roc}
\alias{plot}
\alias{plot.roc}
\alias{plot.roc.roc}
\alias{plot.smooth.roc}
\alias{plot.roc.smooth.roc}
\alias{plot.roc.default}
\alias{plot.roc.formula}

\title{
  Plot a ROC curve
}
\description{
  This function plots a ROC curve. It can accept many arguments to
  tweak the appearance of the plot. Two syntaxes are possible: one
  object of class \dQuote{\link{roc}}, or either two vectors (response, predictor) or a
  formula (response~predictor) as in the \code{\link{roc}} function. 
}
\usage{
\S3method{plot}{roc}(x, ...)
\S3method{plot}{smooth.roc}(x, ...)
\S3method{plot.roc}{roc}(x, add=FALSE, reuse.auc=TRUE,
axes=TRUE, legacy.axes=FALSE,
# Generic arguments for par:
xlim=if(x$percent){c(100, 0)} else{c(1, 0)},
ylim=if(x$percent){c(0, 100)} else{c(0, 1)},
xlab=ifelse(x$percent, ifelse(legacy.axes, "100 - Specificity (\%)", "Specificity (\%)"),
            ifelse(legacy.axes, "1 - Specificity", "Specificity")),
ylab=ifelse(x$percent, "Sensitivity (\%)", "Sensitivity"),
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
print.thres.pattern=ifelse(x$percent, "\%.1f (\%.1f\%\%, \%.1f\%\%)", "\%.3f (\%.3f, \%.3f)"),
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
grid.v={if(is.logical(grid) && grid[1]==TRUE)
          {seq(0, 1, 0.1) * ifelse(x$percent, 100, 1)}
        else if(is.numeric(grid)) 
          {seq(0, ifelse(x$percent, 100, 1), grid[1])} else {NULL}},
grid.h={if (length(grid) == 1) {grid.v} 
        else if (is.logical(grid) && grid[2]==TRUE)
          {seq(0, 1, 0.1) * ifelse(x$percent, 100, 1)} 
        else if(is.numeric(grid))
          {seq(0, ifelse(x$percent, 100, 1), grid[2])} else {NULL}},
grid.lty=3,
grid.lwd=1,
grid.col="#DDDDDD",
# Polygon for the AUC
auc.polygon=FALSE,
auc.polygon.col="gainsboro",
auc.polygon.lty=par("lty"),
auc.polygon.density=NULL,
auc.polygon.angle=45,
auc.polygon.border=NULL,
# Polygon for the maximal AUC possible                           
max.auc.polygon=FALSE,
max.auc.polygon.col="#EEEEEE", 
max.auc.polygon.lty=par("lty"),
max.auc.polygon.density=NULL,
max.auc.polygon.angle=45,
max.auc.polygon.border=NULL,
# Confidence interval
ci=!is.null(x$ci),
ci.type=c("bars", "shape", "no"),
ci.col=ifelse(ci.type=="bars", par("fg"), "gainsboro"),
...)
\S3method{plot.roc}{formula}(x, data, subset, na.action, ...)
\S3method{plot.roc}{default}(x, predictor, ...)
\S3method{plot.roc}{smooth.roc}(x, ...)

}
		   
\arguments{
  \item{x}{a roc object from the \link{roc} function (for plot.roc.roc),
    a formula (for plot.roc.formula) or a response vector (for
    plot.roc.default).
  }
  \item{predictor, data}{arguments for the \link{roc} function.}
  \item{subset,na.action}{arguments for \code{\link{model.frame}}}
  \item{add}{if TRUE, the ROC curve will be added to an existing
    plot. If FALSE (default), a new plot will be created.
  }
  \item{reuse.auc}{if \code{TRUE} (default) and the \dQuote{roc} object
    contains an \dQuote{auc} field, re-use these specifications for the
    plot (specifically \code{print.auc}, \code{auc.polygon} and
    \code{max.auc.polygon} arguments). See details.
  }
  \item{axes}{a logical indicating if the plot axes must be drawn.}
  \item{legacy.axes}{a logical indicating if the specificity axis (x
    axis) must be plotted as as decreasing \dQuote{specificity}
    (\code{FALSE}, the default) or increasing \dQuote{1 - specificity}
    (\code{TRUE}) as in most legacy software. This affects only the
    axis, not the plot coordinates.}
  \item{xlim, ylim, xlab, ylab, asp, mar, mgp}{Generic arguments for the
	plot. See \link{plot} and \link{plot.window} for more details. Only
	used if \code{add=FALSE}.
  }
  \item{col,lty, lwd}{color, line type and line width for the ROC
    curve. See \link{par} for more details.
  }
  \item{type}{type of plotting as in \code{\link{plot}}.}
  \item{identity}{logical: whether or not the identity line (no discrimination
    line) must be displayed. Default: only on new plots.
  }
  \item{identity.col, identity.lty, identity.lwd}{color, line type and
    line width for the identity line. Used only if identity=TRUE. See
    \link{par} for more details.
  }
  \item{print.thres}{Should a selected set of thresholds be displayed on
	the ROC curve? \code{FALSE}, \code{NULL} or \dQuote{no}: no threshold is
	displayed. \code{TRUE} or \dQuote{best}: the threshold with the
	highest sum sensitivity + specificity is plotted (this might be more
	than one threshold). \dQuote{all}: all the points of the ROC
	curve. \dQuote{local maximas}: all the local maximas. Numeric
	vector: direct definition of the thresholds to display.
	Note that on a smoothed ROC curve, only \dQuote{best} is supported.
  }
  \item{print.thres.pch, print.thres.adj, print.thres.col,
    print.thres.cex}{the plotting character (pch), text string
    adjustment (adj), color (col) and character expansion factor (cex)
    parameters for the printing of the thresholds. See \link{points} and
    \link{par} for more details.
  }
  \item{print.thres.pattern}{the text pattern for the thresholds, as a
    \link{sprintf} format. Three numerics are passed to sprintf:
    threshold, specificity, sensitivity.
  }
  \item{print.thres.pattern.cex}{the character expansion factor (cex) for the
    threshold text pattern. See \link{par} for more details.
  }
  \item{print.thres.best.method, print.thres.best.weights}{if
    \code{print.thres="best"} or \code{print.thres=TRUE}, what method must be used to determine which
    threshold is the best. See argument \code{best.method} and \code{best.weights} to
    \code{\link{coords}} for more details.
  }
  \item{print.auc}{boolean. Should the numeric value of AUC be printed
    on the plot?
  }
  \item{print.auc.pattern}{the text pattern for the AUC, as a
    \link{sprintf} format. If NULL, a reasonable value is computed that
    takes partial AUC, CI and percent into account. If the CI
    of the AUC was computed, three numerics are passed to
    sprintf: AUC, lower CI bound, higher CI bound. Otherwise, only AUC is
    passed.
  }
  \item{print.auc.x, print.auc.y}{x and y position for the printing of
    the AUC.
  }
  \item{print.auc.adj, print.auc.cex, print.auc.col}{the text
    adjustment, character expansion factor and color for the printing of
    the AUC. See \link{par} for more details.
  }
  \item{grid}{boolean or numeric vector of length 1 or 2. Should a
    background grid be added to the plot? Numeric: show a grid with the
    specified interval between each line; Logical: show the grid or
    not. Length 1: same values are taken for horizontal and vertical
    lines. Length 2: grid value for vertical (grid[1]) and horizontal
    (grid[2]). Note that these values are used to compute grid.v and
    grid.h. Therefore if you specify a grid.h and grid.v, it will be
    ignored. 
  }
  \item{grid.v, grid.h}{numeric. The x and y values at which a vertical
    or horizontal line (respectively) must be drawn. NULL if no line
    must be added.
  }
  \item{grid.lty, grid.lwd, grid.col}{the line type (lty), line width (lwd) and
    color (col) of the lines of the grid. See \link{par} for
    more details. Note that you can pass vectors of length 2, in which
    case it specifies the vertical (1) and horizontal (2) lines.
  }
  \item{auc.polygon}{boolean. Whether or not to display the area as a
    polygon.
  }
  \item{auc.polygon.col, auc.polygon.lty, auc.polygon.density,
    auc.polygon.angle, auc.polygon.border}{color (col), line type
    (lty), density, angle and border for the AUC polygon. See
    \code{\link{polygon}} and \code{\link{par}} for more details.
  }
  \item{max.auc.polygon}{boolean. Whether or not to display the maximal
    possible area as a polygon.
  }
  \item{max.auc.polygon.col, max.auc.polygon.lty, max.auc.polygon.density,
    max.auc.polygon.angle, max.auc.polygon.border}{color (col), line type
    (lty), density, angle and border for the maximum AUC polygon. See
    \code{\link{polygon}} and \code{\link{par}} for more details.
  }
  \item{ci}{boolean. Should we plot the confidence intervals?}
  \item{ci.type, ci.col}{\code{type} and \code{col} arguments for
	\code{\link{plot.ci}}. The special value \dQuote{no} disables the plotting of confidence intervals.
  }
  \item{\dots}{further arguments passed to or from other methods,
    especially arguments for \code{\link{roc}} and \code{plot.roc.roc} when calling
    \code{plot.roc.default} or \code{plot.roc.formula}. Note that the
    \code{plot} argument for \code{\link{roc}} is not allowed.
    Arguments for \code{\link{auc}} and graphical functions
    \code{\link{plot}}, \code{\link{abline}}, \code{\link{polygon}},
    \code{\link{points}}, \code{\link{text}} and \code{\link{plot.ci}}
    if applicable.
  }
}

\details{
  This function is typically called from \code{\link{roc}} when plot=TRUE (not by
  default). \code{plot.roc.formula} and \code{plot.roc.default} are convenience methods
  that build the ROC curve (with the \code{\link{roc}} function) before
  calling \code{plot.roc.roc}. You can pass them arguments for both
  \code{\link{roc}} and \code{plot.roc.roc}. Simply use \code{plot.roc}
  that will dispatch to the correct method.

  The plotting is done in the following order:
  \enumerate{
    \item A new plot is created if \code{add=FALSE}.
    \item The grid is added if \code{grid.v} and \code{grid.h} are not NULL.
    \item The maximal AUC polygon is added if \code{max.auc.polygon=TRUE}.
	\item The CI shape is added if \code{ci=TRUE}, \code{ci.type="shape"} and \code{x$ci} isn't a \dQuote{ci.auc}.
    \item The AUC polygon is added if \code{auc.polygon=TRUE}.
    \item The identity line if \code{identity=TRUE}.
    \item The actual ROC line is added.
	\item The CI bars are added if \code{ci=TRUE}, \code{ci.type="bars"} and \code{x$ci} isn't a \dQuote{ci.auc}.
    \item The selected thresholds are printed if \code{print.thres} is \code{TRUE} or numeric.
    \item The AUC is printed if \code{print.auc=TRUE}.
  }

  Graphical functions are called with \link{suppressWarnings}.
}

\section{AUC specification}{
  For \code{print.auc}, \code{auc.polygon} and \code{max.auc.polygon}
  arguments, an AUC specification is
  required. By default, the total AUC is plotted, but you may want a
  partial AUCs. The specification is defined by:
  \enumerate{
    \item the \dQuote{auc} field in the \dQuote{\link{roc}} object if
    \code{reuse.auc} is set to \code{TRUE} (default). It is naturally
    inherited from any call to \code{\link{roc}} and fits most cases.
    \item passing the specification to \code{\link{auc}} with \dots
    (arguments \code{partial.auc}, \code{partial.auc.correct} and
    \code{partial.auc.focus}). In this case, you must ensure either that
    the \code{\link{roc}} object do not contain an \code{auc} field (if
    you called \code{\link{roc}} with \code{auc=FALSE}), or set
    \code{reuse.auc=FALSE}.
  }

  If \code{reuse.auc=FALSE} the \code{\link{auc}} function will always
  be called with \code{\dots} to determine the specification, even if
  the \dQuote{\link{roc}} object do contain an \code{auc} field.

  As well if the \dQuote{\link{roc}} object do not contain an \code{auc}
  field, the \code{\link{auc}} function will always be called with
  \code{\dots} to determine the specification.

  Warning: if the roc object passed to plot.roc contains an \code{auc}
  field and \code{reuse.auc=TRUE}, \link{auc} is not called and
  arguments such as \code{partial.auc} are silently ignored.
  
}

\value{
This function returns a list of class \dQuote{roc} invisibly. See \link{roc} for more details.
}

\references{
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
}

\seealso{
\code{\link{roc}}, \code{\link{auc}}, \code{\link{ci}}
}
\examples{

# Create a few ROC curves:
data(aSAH)
roc.s100b <- roc(aSAH$outcome, aSAH$s100b)
roc.wfns <- roc(aSAH$outcome, aSAH$wfns)
roc.ndka <- roc(aSAH$outcome, aSAH$wfns)

# Simple example:
plot(roc.s100b)

# Add a smoothed ROC:
plot(smooth(roc.s100b), add=TRUE, col="blue")
legend("bottomright", legend=c("Empirical", "Smoothed"),
       col=c(par("fg"), "blue"), lwd=2)

# With more options:
plot(roc.s100b, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue", print.thres=TRUE)

# To plot a different partial AUC, we need to ignore the existing value
# with reuse.auc=FALSE:
plot(roc.s100b, print.auc=TRUE, auc.polygon=TRUE, partial.auc=c(1, 0.8),
     partial.auc.focus="se", grid=c(0.1, 0.2), grid.col=c("green", "red"),
     max.auc.polygon=TRUE, auc.polygon.col="lightblue", 
     print.thres=TRUE, print.thres.adj = c(1, -1),
     reuse.auc=FALSE)

# Add a second ROC curve to the previous plot:
plot(roc.wfns, add=TRUE)

# Plot some thresholds, add them to the same plot
plot(roc.ndka, print.thres="best", print.thres.best.method="youden")
plot(roc.ndka, print.thres="best", print.thres.best.method="closest.topleft",
                                   add = TRUE)
plot(roc.ndka, print.thres="best", print.thres.best.method="youden",
                                   print.thres.best.weights=c(50, 0.2),
                                   print.thres.adj = c(1.1, 1.25),
                                   add = TRUE)
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}

