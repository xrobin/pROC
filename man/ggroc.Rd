\encoding{UTF-8}
\name{ggroc.roc}
\alias{ggroc.roc}
\alias{ggroc.smooth.roc}
\alias{ggroc.list}
\alias{ggroc}

\title{
  Plot a ROC curve with ggplot2
}
\description{
  This function plots a ROC curve with ggplot2.
}
\usage{
\S3method{ggroc}{roc}(data, legacy.axes = FALSE, ...)
\S3method{ggroc}{smooth.roc}(data, legacy.axes = FALSE, ...)
\S3method{ggroc}{list}(data, aes = c("colour", "alpha", "linetype", "linewidth", "group"),
      legacy.axes = FALSE, ...)
}
		   
\arguments{
  \item{data}{
    a roc object from the \link{roc} or \link{smooth.roc} function,
    or a list of \code{roc} and/or \code{smooth.roc} objects.
  }
  \item{aes}{the name(s) of the aesthetics for \code{\link[ggplot2:geom_path]{geom_line}} 
    to map to the different ROC curves supplied.
  }
  \item{legacy.axes}{a logical indicating if the specificity axis (x
    axis) must be plotted as as decreasing \dQuote{specificity}
    (\code{FALSE}, the default) or increasing \dQuote{1 - specificity}
    (\code{TRUE}) as in most legacy software.
  }
  \item{...}{additional aesthetics for \code{\link[ggplot2:geom_path]{geom_line}}
    to set: \code{alpha}, \code{colour}, \code{linetype} and \code{linewidth}.
    \code{size} is deprecated and converted to \code{linewidth} with a warning.
    The \code{group} aesthetic is always active since version 1.18.5 and is
    kept for backwards compatibility.
  }
}

\details{
  This function initializes a ggplot object from a ROC curve (or multiple if a list is passed).
  
  It returns the ggplot with a line layer and axis labels.
  You can print it directly or add your own layers and theme elements.

  It also sets \code{\link[ggplot2:coord_fixed]{coord_fixed}()}, so one
  unit of specificity equals one unit of sensitivity and distances to
  the identity line stay meaningful.
  Add another coordinate system to replace it (e.g.
  \code{\link[ggplot2:coord_cartesian]{coord_cartesian}()}).
}

\section{Lists and facets}{
  With a list, the curves are named by the list (unnamed elements by
  their index) and the plot has a \code{name} column, so that
  \code{facet_grid(. ~ name)} draws each curve in its own panel.

  The \code{geom_*} layers also accept a list, named the same way. Their
  data get the same \code{name} column, so each curve's polygon, point
  or label is drawn in the panel of its curve. The other arguments of
  the layer apply to all the curves, and its \code{mapping} can map
  aesthetics to \code{name} through the scales of the plot, for instance
  \code{aes(colour = name)}. Build lists of \code{\link{auc}},
  \code{\link{ci.auc}} or \code{\link{ci.se}} objects from a list of
  curves with \code{\link{lapply}}, as in the examples.
}

\section{Colours of the layers}{
  The polygons and error bars of a curve (\code{\link{geom_polygon_auc}},
  \code{\link{geom_polygon_ci}} and \code{\link{geom_errorbar_ci}})
  take the colour of their curve: the colour that \code{ggroc} maps to
  the curves of a list, otherwise the fixed colour of the ROC line,
  otherwise the theme's ink
  (\code{\link[ggplot2:element_geom]{element_geom}}). Polygons are
  filled with it at 25\% opacity. The threshold points and AUC text
  (\code{\link{geom_roc_threshold}} and \code{\link{geom_text_auc}})
  are in the theme's ink, which stands out over the polygons. A
  \code{colour} (\code{fill} for the polygons), set or mapped, replaces
  these defaults.
}


\seealso{
\code{\link{roc}}, \code{\link{plot.roc}}, \code{\link{geom_polygon_auc}},
\code{\link{geom_roc_identity}}, \code{\link{geom_polygon_max_auc}},
\code{\link{geom_text_auc}}, \code{\link{geom_roc_threshold}},
\code{\link{geom_polygon_ci}}, \code{\link{geom_errorbar_ci}},
\pkg{\link[ggplot2:ggplot2-package]{ggplot2}}
}
\examples{
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)

if (requireNamespace("ggplot2", quietly = TRUE)) {
library(ggplot2)

# A single curve
g <- ggroc(roc1)
g

# With fixed aesthetics
ggroc(roc1, alpha = 0.5, colour = "red", linetype = 2, linewidth = 2)

# With 1 - specificity on the x axis
ggroc(roc1, legacy.axes = TRUE)

# The result is a regular ggplot, with your own theme, title, etc.
g + theme_minimal() + ggtitle("My ROC curve")

# All the layers that add to a ggroc plot: the partial AUC between 90 and
# 100% specificity, the identity line, the best threshold within it, and
# confidence intervals (see their help pages)
auc1 <- auc(roc1, partial.auc = c(1, 0.9))
ci_auc <- ci(auc1, boot.n = 100)
ci_se <- ci.se(roc1, specificities = seq(0, 1, 0.05), boot.n = 100)
ci_th <- ci.thresholds(auc1, thresholds = "best", boot.n = 100)
ggroc(roc1) +
  geom_polygon_max_auc(auc1) +
  geom_polygon_auc(auc1) +
  geom_roc_identity(roc1) +
  geom_polygon_ci(ci_se, fill = "#1c61b6", alpha = 0.3) +
  geom_errorbar_ci(ci_th) +
  geom_roc_threshold(auc1, thresholds = "best") +
  geom_text_auc(ci_auc)

# Multiple curves, from a list or from a formula
rocs <- roc(outcome ~ s100b + ndka + wfns, data = aSAH)
g.list <- ggroc(rocs)
g.list
ggroc(list(empirical = roc1, binormal = smooth(roc1)))

# Map other aesthetics to the curves
ggroc(rocs, aes = "linetype", colour = "red")
ggroc(rocs, aes = c("linetype", "colour"))

# Change the scales as with any ggplot
g.list + scale_colour_brewer(palette = "RdGy")

# Facets
g.list + facet_grid(. ~ name) + theme(legend.position = "none")

# Layers in facets: lists of curves, AUCs or CIs, built with lapply.
# The polygons take the colour of their curve.
aucs <- lapply(rocs, auc, partial.auc = c(1, 0.9))
ci_aucs <- lapply(aucs, ci, boot.n = 100)
ci_ses <- lapply(rocs, ci.se, specificities = seq(0, 1, 0.05), boot.n = 100)
g.list + facet_grid(. ~ name) +
  geom_polygon_max_auc(aucs) +
  geom_polygon_auc(aucs) +
  geom_polygon_ci(ci_ses) +
  geom_roc_threshold(aucs, thresholds = "best") +
  geom_text_auc(ci_aucs)
}
}
