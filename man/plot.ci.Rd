\encoding{UTF-8}
\name{plot.ci}
\alias{plot.ci}
\alias{plot.ci.thresholds}
\alias{plot.ci.se}
\alias{plot.ci.sp}

\title{
  Plot confidence intervals
}
\description{
  This function adds confidence intervals to a ROC curve plot, either as
  bars or as a confidence shape.
}
\usage{
\S3method{plot}{ci.thresholds}(x, length=.01*ifelse(attr(x,
  "roc")$percent, 100, 1), col=par("fg"), ...)
\S3method{plot}{ci.sp}(x, type=c("bars", "shape"), length=.01*ifelse(attr(x,
"roc")$percent, 100, 1), col=ifelse(type=="bars", par("fg"),
"gainsboro"), no.roc=FALSE, ...)
\S3method{plot}{ci.se}(x, type=c("bars", "shape"), length=.01*ifelse(attr(x,
"roc")$percent, 100, 1), col=ifelse(type=="bars", par("fg"),
"gainsboro"), no.roc=FALSE, ...)
}
		   
\arguments{
  \item{x}{a confidence interval object from the functions
    \code{\link{ci.thresholds}}, \code{\link{ci.se}} or \code{\link{ci.sp}}.
  }
  \item{type}{type of plot, \dQuote{bars} or \dQuote{shape}. Can be
    shortened to \dQuote{b} or \dQuote{s}. \dQuote{shape} is only available for
    \code{ci.se} and \code{ci.sp}, not for \code{ci.thresholds}.
  }
  \item{length}{the length (as plot coordinates) of the bar ticks. Only
    if \code{type="bars"}.
  }
  \item{no.roc}{
    if \code{FALSE}, the ROC line is re-added over the
    shape. Otherwise if \code{TRUE}, only the shape is plotted. Ignored
    if \code{type="bars"}
  }
  \item{col}{color of the bars or shape.}
  \item{\dots}{further arguments for \code{\link{segments}} (if
    \code{type="bars"}) or \code{\link{polygon}} (if
    \code{type="shape"}).
  }
}

\details{
  This function adds confidence intervals to a ROC curve plot, either as
  bars or as a confidence shape, depending on the state of the
  \code{type} argument. The shape is plotted over the ROC curve, so that
  the curve is re-plotted unless \code{no.roc=TRUE}.
  
  Graphical functions are called with \link{suppressWarnings}.
}

\section{Warnings}{
  With \code{type="shape"}, the warning \dQuote{Low definition shape} is
  issued when the shape is defined by less than 15 confidence
  intervals. In such a case, the shape is not well defined and the ROC
  curve could pass outside the shape. To get a better shape, increase
  the number of intervals, for example with:
  \preformatted{plot(ci.sp(rocobj, sensitivities=seq(0, 1, .01)), type="shape")}
}

\value{
This function returns the confidence interval object invisibly.
}

\references{
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
}

\seealso{
\code{\link{plot.roc}}, \code{\link{ci.thresholds}}, \code{\link{ci.sp}}, \code{\link{ci.se}}
}
\examples{
data(aSAH)
\dontrun{
# Start a ROC plot
rocobj <- plot.roc(aSAH$outcome, aSAH$s100b)
plot(rocobj)
# Thresholds
ci.thresolds.obj <- ci.thresholds(rocobj)
plot(ci.thresolds.obj)
# Specificities
plot(rocobj) # restart a new plot
ci.sp.obj <- ci.sp(rocobj, boot.n=500)
plot(ci.sp.obj)
# Sensitivities
plot(rocobj) # restart a new plot
ci.se.obj <- ci(rocobj, of="se", boot.n=500)
plot(ci.se.obj)

# Plotting a shape. We need more
ci.sp.obj <- ci.sp(rocobj, sensitivities=seq(0, 1, .01), boot.n=100)
plot(rocobj) # restart a new plot
plot(ci.sp.obj, type="shape", col="blue")

# Direct syntax (response, predictor):
plot.roc(aSAH$outcome, aSAH$s100b,
         ci=TRUE, of="thresholds")
}
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}

