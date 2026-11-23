\encoding{UTF-8}
\name{geom_polygon_max_auc}
\alias{geom_polygon_max_auc}
\alias{geom_polygon_max_auc.auc}
\alias{geom_polygon_max_auc.roc}
\alias{geom_polygon_max_auc.smooth.roc}
\alias{geom_polygon_max_auc.list}

\title{
  Add the maximal AUC rectangle to a ggroc plot
}
\description{
  Add the 100\% AUC (or partial AUC) rectangle, matching
  \code{\link{plot.roc}(max.auc.polygon = TRUE)}.
}
\usage{
geom_polygon_max_auc(data, ...)
\S3method{geom_polygon_max_auc}{auc}(data, fill = NULL, colour = NA, mapping = NULL, ...)
\S3method{geom_polygon_max_auc}{roc}(data, ...)
\S3method{geom_polygon_max_auc}{smooth.roc}(data, ...)
\S3method{geom_polygon_max_auc}{list}(data, ...)
}

\arguments{
  \item{data}{a \code{\link{roc}} or \code{\link{auc}} object. Partial
    AUC bounds are taken from the AUC.
    Or a list of these, named as in \code{\link{ggroc}}, to draw one layer
    per curve (see \sQuote{Lists and facets} in \code{\link{ggroc}}).
  }
  \item{fill, colour}{fill and border for the rectangle. If \code{fill} is
    \code{NULL}, the theme's ink and paper colours
    (\code{\link[ggplot2:element_geom]{element_geom}}) are mixed 80\% of
    the way from ink to paper.
  }
  \item{mapping}{additional aesthetic mappings, created with
    \code{\link[ggplot2:aes]{aes}()}: see \sQuote{Lists and facets} in
    \code{\link{ggroc}}.
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
roc1 <- roc(aSAH$outcome, aSAH$s100b)
auc_sp <- auc(roc1, partial.auc = c(0.9, 1))
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggroc(roc1) +
    geom_polygon_max_auc(auc_sp) +
    geom_polygon_auc(auc_sp)
}
}
