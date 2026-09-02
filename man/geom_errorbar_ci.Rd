\encoding{UTF-8}
\name{geom_errorbar_ci}
\alias{geom_errorbar_ci}
\alias{geom_errorbar_ci.ci.se}
\alias{geom_errorbar_ci.ci.sp}
\alias{geom_errorbar_ci.ci.thresholds}

\title{
  Add CI error bars to a ggroc plot
}
\description{
  Plot a \code{\link{ci.se}}, \code{\link{ci.sp}}, or
  \code{\link{ci.thresholds}} object as error bars, matching
  \code{\link{plot.ci}(type = "bars")} and
  \code{\link{plot.ci.thresholds}}.
}
\usage{
\S3method{geom_errorbar_ci}{ci.se}(data, ...)
\S3method{geom_errorbar_ci}{ci.sp}(data, ...)
\S3method{geom_errorbar_ci}{ci.thresholds}(data, ...)
}

\arguments{
  \item{data}{a \code{\link{ci.se}}, \code{\link{ci.sp}}, or
    \code{\link{ci.thresholds}} object.
  }
  \item{\dots}{additional aesthetics for \code{\link[ggplot2:geom_errorbar]{geom_errorbar}}.
  }
}

\details{
  Drawn in front of the ROC curve. \code{ci.se} bars are vertical,
  \code{ci.sp} bars are horizontal, and \code{ci.thresholds} draws a
  sensitivity/specificity cross. It follows the parent
  \code{\link{ggroc}} plot's axes; there is no \code{legacy.axes}
  argument.
}

\seealso{
  \code{\link{ggroc}}, \code{\link{geom_polygon_ci}},
  \code{\link{geom_roc_threshold}}, \code{\link{ci.se}},
  \code{\link{ci.sp}}, \code{\link{ci.thresholds}}, \code{\link{plot.ci}}
}
\examples{
\dontrun{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)
ci_sp <- ci.sp(rocobj, sensitivities = seq(0, 100, 5), boot.n = 100)
ci_th <- ci.thresholds(rocobj, thresholds = "best", boot.n = 100)
if (require(ggplot2)) {
  ggroc(rocobj) + geom_errorbar_ci(ci_sp)
  ggroc(rocobj) +
    geom_errorbar_ci(ci_th) +
    geom_roc_threshold(rocobj, thresholds = "best")
}
}
}
