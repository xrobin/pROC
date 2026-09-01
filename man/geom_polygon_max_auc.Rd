\encoding{UTF-8}
\name{geom_polygon_max_auc}
\alias{geom_polygon_max_auc}
\alias{geom_polygon_max_auc.auc}
\alias{geom_polygon_max_auc.roc}
\alias{geom_polygon_max_auc.smooth.roc}

\title{
  Add the maximal AUC rectangle to a ggroc plot
}
\description{
  Add the 100\% AUC (or partial AUC) rectangle, matching
  \code{\link{plot.roc}(max.auc.polygon = TRUE)}.
}
\usage{
\S3method{geom_polygon_max_auc}{auc}(data, fill = "#EEEEEE", colour = NA, ...)
\S3method{geom_polygon_max_auc}{roc}(data, ...)
\S3method{geom_polygon_max_auc}{smooth.roc}(data, ...)
}

\arguments{
  \item{data}{a \code{\link{roc}} or \code{\link{auc}} object. Partial
    AUC bounds are taken from the AUC.
  }
  \item{fill, colour}{fill and border for the rectangle.
  }
  \item{\dots}{additional aesthetics for \code{\link[ggplot2:geom_polygon]{geom_polygon}}.
  }
}

\details{
  Inserted behind the ROC curve. For a partial AUC, the rectangle covers
  the pAUC band (specificity or sensitivity focus). It follows the parent
  \code{\link{ggroc}} plot's axes; there is no \code{legacy.axes} argument.
  Add this layer before \code{\link{geom_polygon_auc}} so the maximal
  rectangle sits behind the actual AUC fill.
}

\seealso{
  \code{\link{ggroc}}, \code{\link{geom_polygon_auc}}, \code{\link{plot.roc}}
}
\examples{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
auc_sp <- auc(rocobj, partial.auc = c(0.9, 1))
if (require(ggplot2)) {
  ggroc(rocobj) +
    geom_polygon_max_auc(auc_sp) +
    geom_polygon_auc(auc_sp)
}
}
