\encoding{UTF-8}
\name{points.roc}
\alias{points.roc}
\alias{points.smooth.roc}
\alias{points.auc}

\title{
  Mark thresholds on a ROC plot
}
\description{
  Mark selected thresholds on an existing base graphics ROC plot,
  matching \code{\link{plot.roc}(print.thres = ...)}. These are
  \code{\link{points}} methods, so simply call \code{points()} on a
  \code{\link{roc}}, \code{\link{smooth.roc}} or \code{\link{auc}}
  object. This is the base graphics twin of
  \code{\link{geom_roc_threshold}}.
}
\usage{
\S3method{points}{roc}(x, thresholds = "best", pattern = NULL, pch = 20,
       adj = c(-.05, 1.25), col = "black", cex = par("cex"),
       pattern.cex = cex, best.method = c("youden", "closest.topleft"),
       best.weights = c(1, 0.5), ...)
\S3method{points}{smooth.roc}(x, thresholds = "best", pattern = NULL, pch = 20,
       adj = c(-.05, 1.25), col = "black", cex = par("cex"),
       pattern.cex = cex, best.method = c("youden", "closest.topleft"),
       best.weights = c(1, 0.5), ...)
\S3method{points}{auc}(x, ...)
}

\arguments{
  \item{x}{a \code{\link{roc}} or \code{\link{smooth.roc}} object, or an
    \code{\link{auc}} object, whose curve is used with its (partial) AUC
    specification: with a partial AUC, the thresholds (including
    \dQuote{best}) are searched within it.
  }
  \item{thresholds}{passed as \code{x} to \code{\link{coords}}:
    \code{"best"}, \code{"all"}, \code{"local maximas"}, numeric
    cut-offs, or ordered/character level labels on an ordered ROC.
    Only \code{"best"} is available on smoothed curves, which have no
    thresholds.
  }
  \item{pattern}{\link{sprintf} pattern for the labels, with the
    threshold, specificity and sensitivity, as
    \code{print.thres.pattern} in \code{\link{plot.roc}}, with the same
    default.
  }
  \item{pch, cex, col}{the plotting character, its expansion factor
    and the color of the points and labels. See \code{\link{par}}.
  }
  \item{adj, pattern.cex}{the justification and expansion factor of
    the labels.
  }
  \item{best.method, best.weights}{if \code{thresholds="best"}, the
    method and weights used to determine which threshold is the best.
    Passed to \code{\link{coords}}; see its arguments
    \code{best.method} and \code{best.weights} for details.
  }
  \item{\dots}{additional arguments passed to \code{\link{points}} and
    \code{\link{text}}. For \code{points.auc}, the arguments of
    \code{points.roc}.
  }
}

\details{
  These functions must be called on an existing plot
  (\code{\link{dev.cur}() == 1} is an error).
}

\value{
  Returns \code{x} invisibly.
}

\seealso{
  \code{\link{plot.roc}}, \code{\link{coords}},
  \code{\link{geom_roc_threshold}}, \code{\link{text.auc}}
}
\examples{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
plot(rocobj)
points(rocobj, "best")

# The best threshold within a partial AUC
pauc <- auc(rocobj, partial.auc = c(1, 0.9))
points(pauc, "best", col = "red")

# Several thresholds
plot(rocobj)
points(rocobj, c(0.2, 0.5), pattern = "\%.1f (\%.2f, \%.2f)")
}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}
