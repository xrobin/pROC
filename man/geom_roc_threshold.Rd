\encoding{UTF-8}
\name{geom_roc_threshold}
\alias{geom_roc_threshold}
\alias{geom_roc_threshold.roc}
\alias{geom_roc_threshold.smooth.roc}
\alias{geom_roc_threshold.auc}
\alias{geom_roc_threshold.list}

\title{
  Print thresholds on a ggroc plot
}
\description{
  Mark selected thresholds on a ggroc plot, matching
  \code{\link{plot.roc}(print.thres = ...)}.
}
\usage{
geom_roc_threshold(data, ...)
\S3method{geom_roc_threshold}{roc}(data, thresholds = "best", pattern = NULL,
      hjust = -0.05, vjust = 1.25, size = 2,
      best.method = c("youden", "closest.topleft"), best.weights = c(1, 0.5),
      mapping = NULL, ...)
\S3method{geom_roc_threshold}{smooth.roc}(data, thresholds = "best", pattern = NULL,
      hjust = -0.05, vjust = 1.25, size = 2,
      best.method = c("youden", "closest.topleft"), best.weights = c(1, 0.5),
      mapping = NULL, ...)
\S3method{geom_roc_threshold}{auc}(data, ...)
\S3method{geom_roc_threshold}{list}(data, ...)
}

\arguments{
  \item{data}{a \code{\link{roc}} or \code{\link{smooth.roc}} object,
    or an \code{\link{auc}} object, whose curve is used with its
    (partial) AUC specification: with a partial AUC, the thresholds
    (including \dQuote{best}) are searched within it.
    Or a list of these, named as in \code{\link{ggroc}}, to draw one layer
    per curve (see \sQuote{Lists and facets} in \code{\link{ggroc}}).
  }
  \item{thresholds}{passed as \code{x} to \code{\link{coords}}:
    \code{"best"}, \code{"all"}, \code{"local maximas"}, numeric
    cut-offs, or ordered/character level labels on an ordered ROC.
  }
  \item{pattern}{sprintf pattern, as \code{print.thres.pattern} in
    \code{\link{plot.roc}}, with the same default.
  }
  \item{hjust, vjust}{text justification, matching
    \code{print.thres.adj = c(-0.05, 1.25)} by default.
  }
  \item{size}{point size.
  }
  \item{best.method, best.weights}{if \code{thresholds="best"}, the
    method and weights used to determine which threshold is the best,
    as \code{print.thres.best.method} and \code{print.thres.best.weights}
    in \code{\link{plot.roc}}. Passed to \code{\link{coords}}; see its
    arguments \code{best.method} and \code{best.weights} for details.
  }
  \item{mapping}{additional aesthetic mappings, created with
    \code{\link[ggplot2:aes]{aes}()}: see \sQuote{Lists and facets} in
    \code{\link{ggroc}}.
  }
  \item{\dots}{additional aesthetics for the point and text layers.
    For \code{geom_roc_threshold.auc}, the arguments of
    \code{geom_roc_threshold.roc}.
  }
}

\details{
  Follows the parent \code{\link{ggroc}} plot's axes.
}

\seealso{
  \code{\link{ggroc}}, \code{\link{coords}}, \code{\link{plot.roc}},
  \code{\link{points.roc}}
}
\examples{
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggroc(roc1) + geom_roc_threshold(roc1, thresholds = "best")
}
}
