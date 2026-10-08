\encoding{UTF-8}
\name{polygon_auc}
\alias{polygon_auc}
\alias{polygon_auc.auc}
\alias{polygon_auc.roc}
\alias{polygon_auc.smooth.roc}
\alias{polygon_max_auc}
\alias{polygon_max_auc.auc}
\alias{polygon_max_auc.roc}
\alias{polygon_max_auc.smooth.roc}

\title{
  Add an AUC polygon to a ROC plot
}
\description{
  Add the AUC (or partial AUC) area, or the 100\% (or partial) maximal
  area rectangle, to an existing base graphics ROC plot. These
  functions let you compose the plot that \code{\link{plot.roc}}
  otherwise draws in one call with \code{auc.polygon = TRUE} and
  \code{max.auc.polygon = TRUE}, from separately-built \code{\link{roc}}
  and \code{\link{auc}} objects.
}
\usage{
\S3method{polygon_auc}{auc}(x, col = "gainsboro", lty = par("lty"), density = NULL,
            angle = 45, border = NULL, ...)
\S3method{polygon_auc}{roc}(x, ...)
\S3method{polygon_auc}{smooth.roc}(x, ...)

\S3method{polygon_max_auc}{auc}(x, col = "#EEEEEE", ...)
\S3method{polygon_max_auc}{roc}(x, ...)
\S3method{polygon_max_auc}{smooth.roc}(x, ...)
}

\arguments{
  \item{x}{a \code{\link{roc}}, \code{\link{smooth.roc}} or \code{\link{auc}}
    object. The \code{roc}/\code{smooth.roc} methods use \code{x$auc}
    (compute one first, with \code{\link{auc}} or \code{roc(..., auc =
    TRUE)}, if missing). The curve used to draw \code{polygon_auc} is
    taken from \code{attr(x, "roc")} on the \code{auc} object, which may
    itself be a \code{smooth.roc}. Partial AUC bounds, for both
    functions, are taken from the \code{auc} object's own attributes.
  }
  \item{col, lty, density, angle, border}{color (col), line type (lty),
    density, angle and border for the polygon. See \code{\link{polygon}}
    and \code{\link{par}} for more details. \code{polygon_max_auc} takes
    only \code{col} directly; pass the others through \code{\dots}.
  }
  \item{\dots}{additional arguments passed to \code{\link{polygon}}.
  }
}

\details{
  These functions must be called on an existing plot
  (\code{\link{dev.cur}() == 1} is an error): either after
  \code{\link{plot}}/\code{\link{plot.roc}}, in which case the polygon
  is drawn \emph{on top} of whatever is already there (it will cover
  the curve, the identity line and any threshold points if they overlap
  the polygon's area), or from a \code{\link{plot.roc}(panel.first =
  ...)} hook, in which case it is drawn \emph{underneath} them, exactly
  where \code{plot.roc(auc.polygon = TRUE, max.auc.polygon = TRUE)}
  would draw it. The latter is almost always what you want; see
  \code{\link{plot.roc}} and the examples below.

  When overlaying several AUC regions on one curve (for example a
  specificity-focused and a sensitivity-focused partial AUC), add
  \code{polygon_max_auc} before \code{polygon_auc} for each region so
  the maximal rectangles stay behind the actual AUC fills.

  \code{multiclass.auc} and \code{mv.multiclass.auc} objects (from
  \code{\link{auc}} applied to a \code{\link{multiclass.roc}}) are not
  supported and raise an error.
}

\value{
  Returns \code{x} invisibly.
}

\seealso{
  \code{\link{plot.roc}}, \code{\link{auc}}, \code{\link[=text.auc]{text}},
  \code{\link{geom_polygon_auc}}, \code{\link{geom_polygon_max_auc}}
}
\examples{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
aucobj <- auc(rocobj)

# Preferred: panel.first draws the polygon underneath the curve
plot(rocobj, panel.first = {
  polygon_max_auc(aucobj)
  polygon_auc(aucobj)
})

# Two partial AUC regions on the same curve
auc_sp <- auc(rocobj, partial.auc = c(1, 0.9), partial.auc.correct = TRUE)
auc_se <- auc(rocobj, partial.auc = c(1, 0.9), partial.auc.correct = TRUE,
              partial.auc.focus = "se")
plot(rocobj, panel.first = {
  polygon_max_auc(auc_se, col = "#00860022")
  polygon_max_auc(auc_sp, col = "#1c61b622")
  polygon_auc(auc_se, col = "#008600")
  polygon_auc(auc_sp, col = "#1c61b6")
})
}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}
