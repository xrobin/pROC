\encoding{UTF-8}
\name{text.auc}
\alias{text.auc}
\alias{text.ci.auc}
\alias{text.roc}
\alias{text.smooth.roc}

\title{
  Print the AUC on a ROC plot
}
\description{
  Add an AUC (or partial AUC), optionally with a confidence interval,
  label to an existing base graphics ROC plot, matching
  \code{\link{plot.roc}(print.auc = TRUE)}. These are \code{\link{text}}
  methods, so simply call \code{text()} on an \code{\link{auc}},
  \code{\link{ci.auc}}, \code{\link{roc}} or \code{\link{smooth.roc}}
  object.
}
\usage{
\S3method{text}{auc}(x, ci = NULL, pattern = NULL, xy = NULL, adj = c(0, 1),
     col = par("col"), cex = par("cex"), ...)
\S3method{text}{ci.auc}(x, ...)
\S3method{text}{roc}(x, ...)
\S3method{text}{smooth.roc}(x, ...)
}

\arguments{
  \item{x}{an \code{\link{auc}}, \code{\link{ci.auc}}, \code{\link{roc}}
    or \code{\link{smooth.roc}} object. The \code{ci.auc} method draws
    \code{attr(x, "auc")} with \code{x} as its confidence interval. The
    \code{roc}/\code{smooth.roc} methods use \code{x$auc} and \code{x$ci}
    (compute an AUC first, with \code{\link{auc}} or \code{roc(...,
    auc = TRUE)}, if missing).
  }
  \item{ci}{optional \code{\link{ci.auc}} object to include in the
    label. Ignored (no interval shown) if it is not of class
    \dQuote{ci.auc}, for example a \code{\link{ci.se}}.
  }
  \item{pattern}{\link{sprintf} pattern for the label. \code{NULL}
    selects the same default as \code{plot.roc(print.auc.pattern =
    NULL)}: \code{"AUC: "} or \code{"Partial AUC: "}, the value, and the
    confidence interval if \code{ci} is given.
  }
  \item{xy}{label position, as \code{c(x, y)}. Defaults to the centre of
    the ROC space (\code{c(50, 50)} if \code{percent}, \code{c(.5, .5)}
    otherwise), matching \code{plot.roc(print.auc.x, print.auc.y)}'s
    default.
  }
  \item{adj, col, cex}{text justification (adj), color (col) and
    character expansion factor (cex). See \code{\link{par}} and
    \code{\link{text}} for more details.
  }
  \item{\dots}{additional arguments passed to \code{\link{text}}.
  }
}

\details{
  Unlike \code{\link{polygon_auc}}, this needs no \code{panel.first}
  hook: the label is drawn on top, so calling it any time after
  \code{\link{plot}}/\code{\link{plot.roc}} is enough.

  These functions must be called on an existing plot
  (\code{\link{dev.cur}() == 1} is an error).

  The color is always taken from \code{col} (default \code{par("col")})
  and never inherited from a curve's own \code{col}; repeat it
  explicitly if you want the label to match a custom curve color.
}

\value{
  Returns the label (a character string) invisibly.
}

\seealso{
  \code{\link{plot.roc}}, \code{\link{auc}}, \code{\link{ci.auc}},
  \code{\link{polygon_auc}}, \code{\link{geom_text_auc}}
}
\examples{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b, ci = TRUE, boot.n = 100)
plot(rocobj)
text(rocobj)                 # "AUC: 0.XXX (0.XXX-0.XXX)"

aucobj <- auc(rocobj)
plot(rocobj)
text(aucobj, xy = c(40, 40), pattern = "Area under the curve: \%.2f")
}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}
