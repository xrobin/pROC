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
\S3method{ggroc}{list}(data, aes = c("colour", "alpha", "linetype", "size", "group"),
      legacy.axes = FALSE, ...)
}
		   
\arguments{
  \item{data}{a roc object from the \link{roc} function, or a list of roc objects.
  }
  \item{aes}{the name(s) of the aesthetics for \code{\link[ggplot2:geom_path]{geom_line}} 
    to map to the different ROC curves supplied. Use \dQuote{group} if
    you want the curves to appear with the same aestetic, for instance if you are 
    faceting instead.
  }
  \item{legacy.axes}{a logical indicating if the specificity axis (x
    axis) must be plotted as as decreasing \dQuote{specificity}
    (\code{FALSE}, the default) or increasing \dQuote{1 - specificity}
    (\code{TRUE}) as in most legacy software.
  }
  \item{...}{additional aesthetics for \code{\link[ggplot2:geom_path]{geom_line}}
    to set: \code{alpha}, \code{colour}, \code{linetype} and \code{size}.
  }
}

\details{
  This function initializes a ggplot object from a ROC curve (or multiple if a list is passed). It returns the ggplot with a line layer on it. You can print it directly or add your own layers and theme elements.
}


\seealso{
\code{\link{roc}}, \code{\link{plot.roc}}, \pkg{\link[ggplot2:ggplot2-package]{ggplot2}}
}
\examples{
# Create a basic roc object
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
rocobj2 <- roc(aSAH$outcome, aSAH$wfns)

if (require(ggplot2)) {
g <- ggroc(rocobj)
g
# with additional aesthetics:
ggroc(rocobj, alpha = 0.5, colour = "red", linetype = 2, size = 2)

# You can then your own theme, etc.
g + theme_minimal() + ggtitle("My ROC curve") + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")

# And change axis labels to FPR/FPR
gl <- ggroc(rocobj, legacy.axes = TRUE)
gl
gl + xlab("FPR") + ylab("TPR") + 
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")

# Multiple curves:
g2 <- ggroc(list(s100b=rocobj, wfns=rocobj2, ndka=roc(aSAH$outcome, aSAH$ndka)))
g2

# This is equivalent to using roc.formula:
roc.list <- roc(outcome ~ s100b + ndka + wfns, data = aSAH)
g.list <- ggroc(roc.list)
g.list

# You can change the aesthetics as you normally would with ggplot2:
g.list + scale_colour_brewer(palette="RdGy")
g.list + scale_colour_manual(values = c("red", "blue", "black"))

# with additional aesthetics:
g3 <- ggroc(roc.list, linetype=2)
g3
g4 <- ggroc(roc.list, aes="linetype", color="red")
g4
# changing multiple aesthetics:
g5 <- ggroc(roc.list, aes=c("linetype", "color"))
g5

# OR faceting
g.list + facet_grid(.~name) + theme(legend.position="none")
# To have all the curves of the same color, use aes="group":
g.group <- ggroc(roc.list, aes="group")
g.group
g.group + facet_grid(.~name)
}
}
