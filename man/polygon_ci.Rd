\encoding{UTF-8}
\name{polygon_ci}
\alias{polygon_ci}
\alias{polygon_ci.ci.se}
\alias{polygon_ci.ci.sp}

\title{
  Add a CI band to a ROC plot
}
\description{
  Plot a \code{\link{ci.se}} or \code{\link{ci.sp}} object as a filled
  band, matching \code{\link{plot.ci}(type = "shape")}. The
  \code{panel.first}-safe counterpart to \code{\link{polygon_auc}}: it
  always draws with \code{no.roc = TRUE}, so it never redraws the curve
  and never loses a custom \code{col}/\code{lwd}.
}
\usage{
\S3method{polygon_ci}{ci.se}(x, col = "gainsboro", ...)
\S3method{polygon_ci}{ci.sp}(x, col = "gainsboro", ...)
}

\arguments{
  \item{x}{a \code{\link{ci.se}} or \code{\link{ci.sp}} object.
  }
  \item{col}{fill for the band. See \code{\link{polygon}} and
    \code{\link{par}} for more details.
  }
  \item{\dots}{additional arguments passed to \code{\link{plot.ci}}
    (and from there to \code{\link{polygon}}).
  }
}

\details{
  Like \code{\link{polygon_auc}}, this must be called on an existing
  plot (\code{\link{dev.cur}() == 1} is an error), and is meant to be
  used from a \code{\link{plot.roc}(panel.first = ...)} hook so the band
  lands underneath the curve. Calling it after \code{plot()} instead
  draws it on top.

  \code{plot.ci(type = "shape")} on its own (\code{no.roc = FALSE}, the
  default there) draws the band and then redraws the curve on top of it
  with default styling, to avoid hiding it; that redraw does not
  preserve a customized \code{col}/\code{lwd} on the original curve.
  \code{polygon_ci} sidesteps this entirely by drawing only the band,
  under \code{panel.first}, so the curve is drawn once, with whatever
  styling you gave it.

  There is no \code{polygon_ci} for \code{\link{ci.thresholds}}: a
  filled band is not defined for it, same as \code{\link{plot.ci}} and
  \code{\link{geom_polygon_ci}}. Its error bars (\code{type = "bars"})
  are drawn in front of the curve by design (narrow marks, not a
  full-coverage shape) and need no \code{panel.first} equivalent; call
  \code{\link{plot.ci}} after \code{\link{plot}} as usual.
}

\value{
  Returns \code{x} invisibly.
}

\seealso{
  \code{\link{plot.roc}}, \code{\link{plot.ci}}, \code{\link{ci.se}},
  \code{\link{ci.sp}}, \code{\link{polygon_auc}},
  \code{\link{geom_polygon_ci}}
}
\examples{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
ci_se <- ci.se(rocobj, specificities = seq(0, 1, .05), boot.n = 100)

# Custom curve styling is preserved: polygon_ci() never redraws the curve
plot(rocobj, col = "red", lwd = 3,
     panel.first = polygon_ci(ci_se, col = "lightblue"))
}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}
