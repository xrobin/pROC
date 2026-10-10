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
\S3method{lines.roc}{formula}(x, data, subset, na.action, ...)
\S3method{lines.roc}{default}(x, predictor, ...)
\S3method{lines.roc}{smooth.roc}(x, ...)

}
		   
\arguments{
  \item{x}{a roc object from the \link{roc} function (for lines.roc.roc),
    a smooth.roc object from the \link[=smooth.roc]{smooth} function (for
    lines.roc.smooth.roc), a formula (for lines.roc.formula) or a response
    vector (for lines.roc.default).
  }
  \item{predictor, data}{arguments for the \link{roc} function.}
  \item{subset,na.action}{arguments for \code{\link{model.frame}}}
  \item{lwd}{line width (see \code{\link{par}}).}
  \item{\dots}{graphical parameters for \code{\link{lines}}, and
    especially \code{type} (see \code{\link{plot.default}}) and
    arguments for \code{\link{par}} such as \code{col} (color),
    \code{lty} (line type) or line characteristics \code{lend},
    \code{ljoin} and \code{lmitre}.
  }
}

\value{
This function returns the ROC curve invisibly: a list of class \dQuote{roc},
or \dQuote{smooth.roc} for a smoothed curve. See \code{\link{roc}} and
\code{\link[=smooth.roc]{smooth}} for more details.
}

\seealso{
\code{\link{roc}}, \code{\link{plot.roc}}
}
\examples{

# Create a few ROC curves:
data(aSAH)
roc.s100b <- roc(aSAH$outcome, aSAH$s100b)
roc.wfns <- roc(aSAH$outcome, aSAH$wfns)

# We need a plot to be ready
plot(roc.s100b, type = "n") # but don't actually plot the curve

# Add the line
lines(roc.s100b, type="b", pch=21, col="blue", bg="grey")

# Add the line of an other ROC curve
lines(roc.wfns, type="o", pch=19, col="red")


# Without using 'lines':
rocobj <- roc(aSAH$outcome, aSAH$s100b)
plot(rocobj, type="b", pch=21, col="blue", bg="grey")

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}

