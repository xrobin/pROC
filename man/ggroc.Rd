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

# All the layers that add to a ggroc plot: the AUC, the identity line,
# a threshold, and confidence intervals (see their help pages)
auc1 <- auc(roc1)
ci_auc <- ci(roc1)
ci_se <- ci.se(roc1, specificities = seq(0, 1, 0.05), boot.n = 100)
ci_th <- ci.thresholds(roc1, thresholds = "best", boot.n = 100)
ggroc(roc1) +
  geom_polygon_max_auc(auc1) +
  geom_polygon_auc(auc1) +
  geom_roc_identity(roc1) +
  geom_polygon_ci(ci_se, fill = "#1c61b6", alpha = 0.3) +
  geom_errorbar_ci(ci_th) +
  geom_roc_threshold(roc1, thresholds = "best") +
  geom_text_auc(auc1, ci = ci_auc)

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
}
}
