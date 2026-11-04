\encoding{UTF-8}
\name{geom_polygon_ci}
\alias{geom_polygon_ci}
\alias{geom_polygon_ci.ci.se}
\alias{geom_polygon_ci.ci.sp}

\title{
  Add a CI band to a ggroc plot
}
\description{
  Plot a \code{\link{ci.se}} or \code{\link{ci.sp}} object as a filled
  band, matching \code{\link{plot.ci}(type = "shape")}.
}
\usage{
\S3method{geom_polygon_ci}{ci.se}(data, fill = "gainsboro", colour = NA, ...)
\S3method{geom_polygon_ci}{ci.sp}(data, fill = "gainsboro", colour = NA, ...)
}

\arguments{
  \item{data}{a \code{\link{ci.se}} or \code{\link{ci.sp}} object.
  }
  \item{fill, colour}{fill and border for the band.
  }
  \item{\dots}{additional aesthetics for \code{\link[ggplot2:geom_polygon]{geom_polygon}}.
  }
}

\details{
  Inserted behind the ROC curve. The default fill is \code{"gainsboro"},
  as in \code{\link{plot.ci}}. It follows the parent \code{\link{ggroc}}
  plot's axes; there is no \code{legacy.axes} argument. Like
  \code{\link{plot.ci}}, a filled band is not available for
  \code{\link{ci.thresholds}}.
}

\seealso{
  \code{\link{ggroc}}, \code{\link{geom_errorbar_ci}},
  \code{\link{ci.se}}, \code{\link{ci.sp}}, \code{\link{plot.ci}}
}
\examples{
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)
\dontrun{
ci_se <- ci.se(roc1, specificities = seq(0, 100, 5))}\dontshow{
ci_se <- ci.se(roc1, specificities = seq(0, 100, 5), boot.n = 10)}
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggroc(roc1) + geom_polygon_ci(ci_se, fill = "#1c61b6", alpha = 0.5)
}
}
