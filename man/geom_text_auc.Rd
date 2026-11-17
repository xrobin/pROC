\encoding{UTF-8}
\name{geom_text_auc}
\alias{geom_text_auc}
\alias{geom_text_auc.auc}
\alias{geom_text_auc.roc}
\alias{geom_text_auc.smooth.roc}
\alias{geom_text_auc.ci.auc}
\alias{geom_text_auc.list}

\title{
  Print the AUC on a ggroc plot
}
\description{
  Add an AUC (or partial AUC) label, matching \code{\link{plot.roc}(print.auc = TRUE)}.
}
\usage{
\S3method{geom_text_auc}{auc}(data, ci = NULL, x = NULL, y = NULL, pattern = NULL,
      hjust = 0, vjust = 1, mapping = NULL, ...)
\S3method{geom_text_auc}{roc}(data, ci = NULL, ...)
\S3method{geom_text_auc}{smooth.roc}(data, ci = NULL, ...)
\S3method{geom_text_auc}{ci.auc}(data, ...)
\S3method{geom_text_auc}{list}(data, ci = NULL, pattern = NULL, name = FALSE, ...)
}

\arguments{
  \item{data}{a \code{\link{roc}}, \code{\link{auc}} or
    \code{\link{ci.auc}} object (the AUC and its CI).
    Or a list of these, named as in \code{\link{ggroc}}, to label each
    curve (see \sQuote{Lists and facets} in \code{\link{ggroc}}).
  }
  \item{ci}{optional \code{\link{ci.auc}} object. Not with a list:
    put the \code{\link{ci.auc}} objects in the list instead.
  }
  \item{x, y}{position in plot coordinates. Defaults to the centre of
    the ROC space (0.5 or 50).
  }
  \item{pattern}{sprintf pattern as in \code{\link{plot.roc}}. \code{NULL}
    selects the same default as \code{print.auc.pattern}.
  }
  \item{name}{with a list, whether to start each label with the name
    of its curve.
  }
  \item{hjust, vjust}{text justification (ggplot2), corresponding to
    \code{print.auc.adj = c(0, 1)} by default.
  }
  \item{mapping}{additional aesthetic mappings, created with
    \code{\link[ggplot2:aes]{aes}()}: see \sQuote{Lists and facets} in
    \code{\link{ggroc}}.
  }
  \item{\dots}{additional arguments to \code{\link[ggplot2:geom_text]{geom_text}}.
  }
}

\details{
  With a list, the labels that share a panel are stacked from
  (\code{x}, \code{y}) down, 1.5 lines apart, which needs a numeric
  \code{vjust}: in facets by \code{name}, each panel shows the label
  of its curve at (\code{x}, \code{y}).

  The default \code{colour} is the theme's ink: see
  \sQuote{Colours of the layers} in \code{\link{ggroc}}.
}

\seealso{
  \code{\link{ggroc}}, \code{\link{plot.roc}}
}
\examples{
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)
ci1 <- ci(roc1)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggroc(roc1) + geom_text_auc(roc1, ci = ci1)
  ggroc(roc1) + geom_text_auc(ci1)

  # Multiple curves, with or without facets
  rocs <- roc(outcome ~ s100b + ndka + wfns, data = aSAH)
  cis <- lapply(rocs, ci)
  ggroc(rocs) + geom_text_auc(cis)
  ggroc(rocs) + geom_text_auc(cis) + ggplot2::facet_grid(. ~ name)
  # Curves that differ by linetype only: name the labels
  ggroc(rocs, aes = "linetype") + geom_text_auc(cis, name = TRUE)
}
}
