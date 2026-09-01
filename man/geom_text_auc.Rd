\encoding{UTF-8}
\name{geom_text_auc}
\alias{geom_text_auc}
\alias{geom_text_auc.auc}
\alias{geom_text_auc.roc}
\alias{geom_text_auc.smooth.roc}

\title{
  Print the AUC on a ggroc plot
}
\description{
  Add an AUC (or partial AUC) label, matching \code{\link{plot.roc}(print.auc = TRUE)}.
}
\usage{
\S3method{geom_text_auc}{auc}(data, ci = NULL, x = NULL, y = NULL, pattern = NULL,
      hjust = 0, vjust = 1, ...)
\S3method{geom_text_auc}{roc}(data, ...)
\S3method{geom_text_auc}{smooth.roc}(data, ...)
}

\arguments{
  \item{data}{a \code{\link{roc}} or \code{\link{auc}} object.
  }
  \item{ci}{optional \code{\link{ci.auc}} object. Taken from
    \code{data$ci} when \code{data} is a ROC curve.
  }
  \item{x, y}{position in plot coordinates. Defaults to the centre of
    the ROC space (0.5 or 50).
  }
  \item{pattern}{sprintf pattern as in \code{\link{plot.roc}}. \code{NULL}
    selects the same default as \code{print.auc.pattern}.
  }
  \item{hjust, vjust}{text justification (ggplot2), corresponding to
    \code{print.auc.adj = c(0, 1)} by default.
  }
  \item{\dots}{additional arguments to \code{\link[ggplot2:annotate]{annotate}("text")}.
  }
}

\seealso{
  \code{\link{ggroc}}, \code{\link{plot.roc}}
}
\examples{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b, ci = TRUE)
if (require(ggplot2)) {
  ggroc(rocobj) + geom_text_auc(rocobj)
}
}
