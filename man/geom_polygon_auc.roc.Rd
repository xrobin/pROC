\encoding{UTF-8}
\name{geom_polygon_auc}
\alias{geom_polygon_auc.auc}
\alias{geom_polygon_auc.roc}
\alias{geom_polygon_auc.smooth.roc}
\alias{geom_polygon_auc}
\alias{ggplot_add.ggroc_layer}

\title{
  Add an AUC polygon to a ggroc plot
}
\description{
  EXPERIMENTAL - Add an AUC polygon to a ggroc plot.
}
\usage{
\S3method{geom_polygon_auc}{auc}(data, ...)
\S3method{geom_polygon_auc}{roc}(data, ...)
\S3method{geom_polygon_auc}{smooth.roc}(data, ...)
}

\arguments{
  \item{data}{a roc object from the \link{roc} function, or its
    \code{auc}, same as the one used to build the ggroc initially.
  }
  \item{...}{additional aesthetics for \code{\link[ggplot2:geom_polygon]{geom_polygon}}
    to set: \code{alpha}, \code{colour}, \code{linetype} and \code{linewidth}.
  }
}

\details{
  This is highly experimental and may change in the future.

  The polygon is inserted immediately before the ROC curve, matching
  \code{\link{plot.roc}(auc.polygon = TRUE)}. It follows the parent
  \code{\link{ggroc}} plot's axes; there is no \code{legacy.axes} argument.
}


\seealso{
\code{\link{ggroc}}, \code{\link{plot.roc}}
}
\examples{

# Create a ROC curve:
data(aSAH)
roc.s100b <- roc(aSAH$outcome, aSAH$s100b)
roc.s100b.percent <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)

if (require(ggplot2)) {
ggroc(roc.s100b) + geom_polygon_auc(roc.s100b$auc)

ggroc(roc.s100b.percent, legacy.axes=TRUE) + geom_polygon_auc(roc.s100b.percent)

# Partial AUCs
auc.s100b.partial.sp <- auc(roc.s100b, partial.auc = c(0.9, 1))
auc.s100b.partial.se <- auc(roc.s100b, partial.auc = c(0.8, 0.9), partial.auc.focus="se")

ggroc(roc.s100b) + geom_polygon_auc(auc.s100b.partial.sp)
ggroc(roc.s100b) + geom_polygon_auc(auc.s100b.partial.se)
}
}
