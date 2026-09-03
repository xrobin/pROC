\encoding{UTF-8}
\name{geom_roc_threshold}
\alias{geom_roc_threshold}
\alias{geom_roc_threshold.roc}
\alias{geom_roc_threshold.smooth.roc}

\title{
  Print thresholds on a ggroc plot
}
\description{
  Mark selected thresholds on a ggroc plot, matching
  \code{\link{plot.roc}(print.thres = ...)}.
}
\usage{
\S3method{geom_roc_threshold}{roc}(data, thresholds = "best", pattern = NULL,
      hjust = -0.05, vjust = 1.25, size = 2, ...)
\S3method{geom_roc_threshold}{smooth.roc}(data, thresholds = "best", pattern = NULL,
      hjust = -0.05, vjust = 1.25, size = 2, ...)
}

\arguments{
  \item{data}{a \code{\link{roc}} or \code{\link{smooth.roc}} object.
  }
  \item{thresholds}{passed as \code{x} to \code{\link{coords}}:
    \code{"best"}, \code{"all"}, \code{"local maximas"}, numeric
    cut-offs, or ordered/character level labels on an ordered ROC.
  }
  \item{pattern}{sprintf pattern as in \code{print.thres.pattern}.
    The default uses \code{\%s} for ordered thresholds and
    \code{\%.3f} (or \code{\%.1f} with \code{percent}) otherwise.
  }
  \item{hjust, vjust}{text justification, matching
    \code{print.thres.adj = c(-0.05, 1.25)} by default.
  }
  \item{size}{point size.
  }
  \item{\dots}{additional aesthetics for the point and text layers.
  }
}

\details{
  Follows the parent \code{\link{ggroc}} plot's axes; there is no
  \code{legacy.axes} argument.
}

\seealso{
  \code{\link{ggroc}}, \code{\link{coords}}, \code{\link{plot.roc}}
}
\examples{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
if (require(ggplot2)) {
  ggroc(rocobj) + geom_roc_threshold(rocobj, thresholds = "best")
}
}
