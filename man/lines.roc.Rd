\encoding{UTF-8}
\name{lines.roc}
\alias{lines.roc}
\alias{lines.roc.roc}
\alias{lines.smooth.roc}
\alias{lines.roc.smooth.roc}
\alias{lines.roc.default}
\alias{lines.roc.formula}

\title{
  Add a ROC line to a ROC plot
}
\description{
  This convenience function adds a ROC line to a ROC curve.
}
\usage{
\S3method{lines}{roc}(x, ...)
\S3method{lines}{smooth.roc}(x, ...)
\S3method{lines.roc}{roc}(x, lwd=2, ...)
\S3method{lines.roc}{formula}(x, data, ...)
\S3method{lines.roc}{default}(x, predictor, ...)
\S3method{lines.roc}{smooth.roc}(x, ...)

}
		   
\arguments{
  \item{x}{a roc object from the \link{roc} function (for plot.roc.roc),
    a formula (for plot.roc.formula) or a response vector (for
    plot.roc.default).
  }
  \item{predictor, data}{arguments for the \link{roc} function.}
  \item{lwd}{line width (see \code{\link{par}}).}
  \item{\dots}{graphical parameters for \code{\link{lines}}, and
    especially \code{type} (see \code{\link{plot.default}}) and
    arguments for \code{\link{par}} such as \code{col} (color),
    \code{lty} (line type) or line characteristics \code{lend},
    \code{ljoin} and \code{lmitre}.
  }
}

\value{
This function returns a list of class \dQuote{roc} invisibly. See \code{\link{roc}} for more details.
}

\seealso{
\code{\link{roc}}, \code{\link{plot.roc}}
}
\examples{
data(aSAH)

rocobj <- plot.roc(aSAH$outcome, aSAH$s100b, type="n")
lines(rocobj, type="b", pch=21, col="blue", bg="grey")

# Without using 'lines':
rocobj <- plot.roc(aSAH$outcome, aSAH$s100b, type="b", pch=21, col="blue", bg="grey")

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}

