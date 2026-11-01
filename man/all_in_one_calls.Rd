\encoding{UTF-8}
\name{all_in_one_calls}
\alias{all_in_one_calls}

\title{All-in-one calls and their modular replacements}

\description{This help page describes the all-in-one calls of pROC, which
are superseded, and the modular calls that replace them.}

\section{Background information}{
Since its first version, pROC can perform many operations in a single
call. \code{\link{roc}} can smooth the curve, compute a partial AUC and its
confidence interval, and plot the result, with one \code{\dots} shared by
\code{\link[=smooth.roc]{smooth}}, \code{\link{auc}}, \code{\link{ci}} and
\code{\link{plot.roc}}:

\preformatted{
roc(aSAH$outcome, aSAH$s100b, percent = TRUE,
    # smoothing
    smooth = TRUE, smooth.method = "density", kernel = "epanechnikov",
    # partial AUC
    auc = TRUE, partial.auc = c(100, 90), partial.auc.focus = "sensitivity",
    partial.auc.correct = TRUE,
    # confidence interval
    ci = TRUE, boot.n = 100, conf.level = 0.9, boot.stratified = FALSE,
    # plot
    plot = TRUE, grid = TRUE, auc.polygon = TRUE, max.auc.polygon = TRUE,
    print.auc = TRUE, print.thres = "best", col = "darkblue", main = "S100B")
}

Likewise, \code{\link{auc}}, \code{\link{ci}}, \code{\link{roc.test}},
\code{\link{plot.roc}} and other functions can build the curve from raw
data, and settings such as the partial AUC or the CI are stored in the
curve and reused later by other functions.

With these calls, an argument can reach the wrong function, or be
partially matched to an argument of another function. Settings stored in
the curve also take precedence over explicit arguments:

\preformatted{
r1 <- roc(aSAH$outcome, aSAH$s100b)
ci.auc(r1, partial.auc = c(1, 0.8))
# 95\% CI: 0.6301-0.8326 (DeLong)
}

This is the CI of the full AUC (0.731), not of the partial AUC: \code{r1}
already contains the full AUC, which \code{ci.auc} reuses
(\code{reuse.auc = TRUE} by default), and \code{partial.auc} is ignored.
}

\section{The modular calls}{
Each function does one thing:
\itemize{
  \item \code{\link{roc}} builds the curve, and computes its full AUC.
  \item \code{\link[=smooth.roc]{smooth}} smooths it.
  \item Only \code{\link{auc}} takes the partial AUC arguments. A partial
    AUC is an \code{auc} object, which is passed on to \code{\link{ci}},
    \code{\link{var}}, \code{\link{cov}}, \code{\link{roc.test}},
    \code{\link{power.roc.test}}, \code{\link{coords}},
    \code{\link{ci.coords}} or \code{\link{ci.thresholds}}.
  \item \code{\link{plot.roc}} draws the curve. Add-on functions draw
    the rest: \code{\link{polygon_auc}}, \code{\link{polygon_max_auc}},
    \code{\link[=text.auc]{text}}, \code{\link[=points.roc]{points}},
    \code{\link{polygon_ci}} and \code{\link[=plot.ci]{plot}} on a CI. Given an
    \code{auc} object, they use its partial AUC.
}

The call above becomes:

\preformatted{
rocobj <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)
smoothed <- smooth(rocobj, method = "density", kernel = "epanechnikov")
pauc <- auc(smoothed, partial.auc = c(100, 90),
            partial.auc.focus = "sensitivity", partial.auc.correct = TRUE)
ciobj <- ci(pauc, boot.n = 100, conf.level = 0.9, boot.stratified = FALSE)
plot(smoothed, grid = TRUE, col = "darkblue", main = "S100B",
     panel.first = {
       polygon_max_auc(pauc)
       polygon_auc(pauc)
     })
points(pauc, "best")
text(ciobj, col = "darkblue")
}

and the partial AUC CI of the second example is:

\preformatted{
ci.auc(auc(r1, partial.auc = c(1, 0.8)))
}
}

\section{Replacements}{
\describe{
  \item{\code{roc(smooth = TRUE)}, with \code{smooth.method},
    \code{smooth.n} and \code{density}}{\code{smooth(rocobj, method =,
    n =, density =)}.}
  \item{\code{roc(ci = TRUE)}, with \code{ci.method}}{\code{ci(rocobj,
    method =)}, or \code{ci(aucobj)} for a partial AUC.}
  \item{\code{roc(plot = TRUE)}}{\code{plot(rocobj)}.}
  \item{The partial AUC arguments (\code{partial.auc},
    \code{partial.auc.focus} and \code{partial.auc.correct}) of
    \code{roc}, \code{multiclass.roc}, \code{ci}, \code{ci.auc},
    \code{var}, \code{cov}, \code{roc.test}, \code{power.roc.test} and
    \code{plot.roc}}{\code{aucobj <- auc(rocobj, partial.auc =)}, then
    pass \code{aucobj} to these functions.}
  \item{\code{reuse.auc}}{pass an \code{auc} object.}
  \item{\code{reuse.auc} and \code{reuse.ci} of
    \code{smooth}}{call \code{auc} and \code{ci} on the smoothed curve.}
  \item{The default and formula methods of \code{auc}, \code{ci},
    \code{ci.*}, \code{roc.test}, \code{plot.roc} and
    \code{lines.roc}}{build the curve with \code{\link{roc}} first.}
  \item{\code{print.auc} of \code{plot.roc}}{\code{text(aucobj)} or
    \code{text(ciobj)}.}
  \item{\code{auc.polygon} and \code{max.auc.polygon} of
    \code{plot.roc}}{\code{polygon_auc(aucobj)} and
    \code{polygon_max_auc(aucobj)}, in \code{panel.first} to draw them
    under the curve.}
  \item{\code{print.thres} of \code{plot.roc}}{\code{points(rocobj,
    thresholds)}, or \code{points(aucobj, thresholds)} to search the
    thresholds within a partial AUC.}
  \item{\code{ci}, \code{ci.type} and \code{ci.col} of
    \code{plot.roc}, and \code{plot(ci, type = "shape")}}{\code{plot(ciobj)},
    \code{polygon_ci(ciobj)} and \code{text(ciobj)}.}
}
The styling arguments of the display arguments (such as
\code{print.auc.cex} or \code{auc.polygon.col}) are the arguments of the
add-on functions (\code{cex} of \code{text}, \code{col} of
\code{polygon_auc}).
}

\section{Reuse of the AUC}{
With a \dQuote{roc} or \dQuote{smooth.roc} curve instead of an \code{auc}
object, \code{\link{ci.auc}}, \code{\link{var}}, \code{\link{cov}},
\code{\link{roc.test}}, \code{\link{power.roc.test}} and
\code{\link{plot.roc}} use the \dQuote{auc} field of the curve if
\code{reuse.auc=TRUE} (default). If the curve has no \dQuote{auc} field, or
with \code{reuse.auc=FALSE}, they call \code{\link{auc}} with their
\code{\dots}, where the partial AUC arguments (\code{partial.auc},
\code{partial.auc.focus} and \code{partial.auc.correct}) can be given.
When the \dQuote{auc} field is used, these arguments are ignored, without
a warning.
}

\section{Changes in 1.20}{
\enumerate{
  \item \code{\link{power.roc.test}}, \code{\link{ci.coords}} and
    \code{\link{ci.thresholds}} accept \code{auc} objects.
  \item \code{\link[=points.roc]{points}} methods mark thresholds, as
    \code{print.thres} does. \code{\link{geom_roc_threshold}} accepts
    \code{auc} objects.
  \item \code{\link[=text.auc]{text}} and \code{\link{geom_text_auc}} on a
    curve no longer use a CI stored in it by default.
  \item The all-in-one calls listed in \sQuote{Replacements} are
    superseded: they keep working without a warning.
}
}

\section{Changes in future versions}{
\enumerate{
  \item The all-in-one calls listed in \sQuote{Replacements}, except the
    display arguments, will be deprecated: they will keep working, with a warning that shows
    the replacement call.
  \item The display arguments of \code{plot.roc}
    (\code{print.auc}, \code{auc.polygon}, \code{max.auc.polygon} and
    \code{print.thres}, with their styling arguments) and \code{plot(ci,
    type = "shape")} will stay superseded, without a warning. They will
    keep working for the full curve.
  \item Later, the deprecated arguments and methods will be removed.
}
}

\section{See also}{
\href{https://github.com/xrobin/pROC/issues/149}{The GitHub issue tracking the changes described in this manual page}.
}
