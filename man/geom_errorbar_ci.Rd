\encoding{UTF-8}
\name{geom_errorbar_ci}
\alias{geom_errorbar_ci}
\alias{geom_errorbar_ci.ci.se}
\alias{geom_errorbar_ci.ci.sp}
\alias{geom_errorbar_ci.ci.thresholds}
\alias{geom_errorbar_ci.list}

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
geom_errorbar_ci(data, ...)
\S3method{geom_errorbar_ci}{ci.se}(data, mapping = NULL, ...)
\S3method{geom_errorbar_ci}{ci.sp}(data, mapping = NULL, ...)
\S3method{geom_errorbar_ci}{ci.thresholds}(data, mapping = NULL, ...)
\S3method{geom_errorbar_ci}{list}(data, ...)
}

\arguments{
  \item{data}{a \code{\link{ci.se}}, \code{\link{ci.sp}}, or
    \code{\link{ci.thresholds}} object.
    Or a list of these, named as in \code{\link{ggroc}}, to draw one layer
    per curve (see \sQuote{Lists and facets} in \code{\link{ggroc}}).
  }
  \item{mapping}{additional aesthetic mappings, created with
    \code{\link[ggplot2:aes]{aes}()}: see \sQuote{Lists and facets} in
    \code{\link{ggroc}}.
  }
  \item{\dots}{additional aesthetics for \code{\link[ggplot2:geom_errorbar]{geom_errorbar}}.
    The default \code{colour} is that of the curve: see
    \sQuote{Colours of the layers} in \code{\link{ggroc}}.
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
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)
\dontrun{
ci_sp <- ci.sp(roc1, sensitivities = seq(0, 100, 5))
ci_th <- ci.thresholds(roc1, thresholds = "best")}\dontshow{
ci_sp <- ci.sp(roc1, sensitivities = seq(0, 100, 5), boot.n = 10)
ci_th <- ci.thresholds(roc1, thresholds = "best", boot.n = 10)}
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggroc(roc1) + geom_errorbar_ci(ci_sp)
  ggroc(roc1) +
    geom_errorbar_ci(ci_th) +
    geom_roc_threshold(roc1, thresholds = "best")
}
}
