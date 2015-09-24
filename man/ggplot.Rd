\encoding{UTF-8}
\name{ggplot.roc}
\alias{ggplot.roc}
\alias{ggplot.list}
\alias{ggplot}

\title{
  Plot a ROC curve with ggplot2
}
\description{
  This function plots a ROC curve with ggplot2.
}
\usage{
\S3method{ggplot}{roc}(data, ...)
\S3method{ggplot}{list}(data, aes = c("colour", "alpha", "linetype", "size"), ...)
}
		   
\arguments{
  \item{data}{a roc object from the \link{roc} function, or a list of roc objects.
  }
  \item{aes}{the name of the aesthetics for \code{\link{geom_line}} to map to the different ROC curves supplied}
  \item{...}{additional aesthetics for \code{\link{geom_line}} to set: \code{alpha}, \code{colour}, \code{linetype} and \code{size}.}
}

\details{
  This function initializes a ggplot object from a ROC curve (or multiple if a list is passed). It returns the ggplot with a line layer on it. You can print it directly or add your own layers and theme elements.
}


\seealso{
\code{\link{roc}}, \code{\link{plot.roc}}, \code{\link{ggplot}}
}
\examples{
# Create a basic roc object
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
rocobj2 <- roc(aSAH$outcome, aSAH$wfns)

library(ggplot2)
g <- ggplot(rocobj)
g
# with additional aesthetics:
ggplot(rocobj, alpha = 0.5, colour = "red", linetype = 2, size = 2)

# You can then your own theme, etc.
g + theme_minimal() + ggtitle("My ROC curve")

# Multiple curves:
g2 <- ggplot(list(s100b=rocobj, wfns=rocobj2))
g2
# with additional aesthetics:
g2 <- ggplot(list(s100b=rocobj, wfns=rocobj2), linetype=2)
ggplot(list(s100b=rocobj, wfns=rocobj2), aes="linetype", color="red")

}


