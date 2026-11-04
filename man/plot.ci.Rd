\encoding{UTF-8}
\name{plot.ci}
\alias{plot.ci}
\alias{plot.ci.thresholds}
\alias{plot.ci.se}
\alias{plot.ci.sp}
\alias{plot.ci.coords}

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
\S3method{plot}{ci.coords}(x, type=c("bars", "shape"), length=NULL,
col=ifelse(type=="bars", par("fg"), "gainsboro"), ...)
}
		   
\arguments{
  \item{x}{a confidence interval object from the functions
    \code{\link{ci.thresholds}}, \code{\link{ci.se}}, \code{\link{ci.sp}} or
    \code{\link{ci.coords}}.
  }
  \item{type}{type of plot, \dQuote{bars} or \dQuote{shape}. Can be
    shortened to \dQuote{b} or \dQuote{s}. \dQuote{shape} is only available for
    \code{ci.se}, \code{ci.sp} and \code{ci.coords}, not for \code{ci.thresholds}.
    For \code{ci.se} and \code{ci.sp}, \dQuote{shape} is superseded by
    \code{\link{polygon_ci}} (see \link{all_in_one_calls}).
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

  With a \code{\link{ci.coords}} object, the coordinates given as
  \code{x} are drawn on the horizontal axis and the confidence interval of
  the first \code{ret} coordinate on the vertical axis (a warning is
  issued if there are several): on a ROC plot this matches the axes only
  for \code{input="specificity"} and \code{ret="sensitivity"}. \code{x}
  must be numeric:
  \dQuote{best} or level labels of an ordered predictor give an error.

  These functions must be called on an existing plot
  (\code{\link{dev.cur}() == 1} is an error).

  Graphical functions are called with \link{suppressWarnings}.
}

\section{Warnings}{
  With \code{type="shape"}, the warning \dQuote{Low definition shape} is
  issued when the shape is defined by less than 15 confidence
  intervals. In such a case, the shape is not well defined and the ROC
  curve could pass outside the shape. To get a better shape, increase
  the number of intervals, for example with:
  \preformatted{plot(ci.sp(roc1, sensitivities=seq(0, 1, .01)), type="shape")}
}

\value{
This function returns the confidence interval object invisibly.
}

\references{
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{12}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
}

\seealso{
\code{\link{plot.roc}}, \code{\link{ci.thresholds}}, \code{\link{ci.sp}}, \code{\link{ci.se}},
\code{\link[=text.auc]{text}} for \code{\link{ci.auc}}, \code{\link{polygon_ci}}
}
\examples{
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)

# CI of thresholds
\dontrun{
ci_th <- ci.thresholds(roc1)}\dontshow{
ci_th <- ci.thresholds(roc1, boot.n = 10)}
plot(roc1)
plot(ci_th)

# CI of specificities, as bars
\dontrun{
ci_sp <- ci.sp(roc1)}\dontshow{
ci_sp <- ci.sp(roc1, boot.n = 10)}
plot(roc1) # restart a new plot
plot(ci_sp)

# As a shape. We need more points
\dontrun{
ci_sp <- ci.sp(roc1, sensitivities = seq(0, 1, .01))}\dontshow{
ci_sp <- ci.sp(roc1, sensitivities = seq(0, 1, .01), boot.n = 10)}
plot(roc1)
plot(ci_sp, type = "shape", col = "blue")

# CI of a precision-recall curve
co <- coords(roc1, x = "all", input = "recall", ret = c("recall", "precision"))
\dontrun{
ci_pr <- ci.coords(roc1, x = seq(0, 1, .05), input = "recall", ret = "precision")}\dontshow{
ci_pr <- ci.coords(roc1, x = seq(0, 1, .05), input = "recall", ret = "precision",
                   boot.n = 10)}
plot(co, type = "l", ylim = c(0, 1))
plot(ci_pr, type = "shape")
lines(co)
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}

