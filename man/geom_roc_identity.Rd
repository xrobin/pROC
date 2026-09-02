\encoding{UTF-8}
\name{geom_roc_identity}
\alias{geom_roc_identity}
\alias{geom_roc_identity.roc}
\alias{geom_roc_identity.smooth.roc}

\title{
  Add the identity line to a ggroc plot
}
\description{
  Add the no-discrimination (identity) line to a ggroc plot, matching
  \code{\link{plot.roc}(identity = TRUE)}.
}
\usage{
\S3method{geom_roc_identity}{roc}(data, colour = "darkgrey", ...)
\S3method{geom_roc_identity}{smooth.roc}(data, colour = "darkgrey", ...)
}

\arguments{
  \item{data}{a \code{\link{roc}} or \code{\link{smooth.roc}} object,
    same as the one passed to \code{\link{ggroc}}.
  }
  \item{colour}{line colour. Default: \code{"darkgrey"}.
  }
  \item{\dots}{additional aesthetics for \code{\link[ggplot2:geom_segment]{geom_segment}}.
  }
}

\details{
  The segment is drawn in data coordinates so it is not clipped by
  \code{\link[ggplot2:coord_fixed]{coord_fixed}()} on a reversed specificity axis.
  It is inserted behind the ROC curve. It follows the parent
  \code{\link{ggroc}} plot's axes; there is no \code{legacy.axes} argument.
}

\seealso{
  \code{\link{ggroc}}, \code{\link{plot.roc}}
}
\examples{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
if (require(ggplot2)) {
  ggroc(rocobj) + geom_roc_identity(rocobj)
}
}
