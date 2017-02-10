\encoding{UTF-8}
\name{ci}
\alias{ci}
\alias{ci.default}
\alias{ci.formula}
\alias{ci.roc}
\alias{ci.smooth.roc}

\title{
  Compute the confidence interval of a ROC curve
}
\description{
  This function computes the confidence interval (CI) of a ROC curve. The
  \code{of} argument controls the type of CI that will be computed.
  By default, the 95\% CI are computed with 2000 stratified bootstrap
  replicates.
}
\usage{
ci(...)
\S3method{ci}{roc}(roc, of = c("auc", "thresholds", "sp", "se", "coords"), ...)
\S3method{ci}{smooth.roc}(smooth.roc, of = c("auc", "sp", "se", "coords"), ...)
\S3method{ci}{formula}(formula, data, ...)
\S3method{ci}{default}(response, predictor, ...)
}
		   
\arguments{
  \item{roc, smooth.roc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, or a \dQuote{smooth.roc} object from the
	\code{\link[=smooth.roc]{smooth}} function.
  }
  \item{response, predictor}{arguments for the \code{\link{roc}} function.}
  \item{formula, data}{a formula (and possibly a data object) of type
	response~predictor for the \code{\link{roc}} function.
  }
  \item{of}{The type of confidence interval. One of \dQuote{auc},
	\dQuote{thresholds}, \dQuote{sp}, \dQuote{se} or \dQuote{coords}. Note that
	confidence interval on \dQuote{thresholds} are not available for
	smoothed ROC curves.
  }
  \item{\dots}{further arguments passed to or from other methods,
    especially \code{\link{auc}}, \code{\link{roc}}, and the specific
	\code{ci} functions \code{\link{ci.auc}}, \code{\link{ci.se}},
	\code{\link{ci.sp}} and \code{\link{ci.thresholds}}.
  }
}

\details{
  \code{ci.formula} and \code{ci.default} are convenience methods
  that build the ROC curve (with the \code{\link{roc}} function) before
  calling \code{ci.roc}. You can pass them arguments for both
  \code{\link{roc}} and \code{ci.roc}. Simply use \code{ci}
  that will dispatch to the correct method.
  
  This function is typically called from \code{\link{roc}} when \code{ci=TRUE} (not by
  default). Depending on the \code{of} argument, the specific
  \code{ci} functions \code{\link{ci.auc}}, \code{\link{ci.thresholds}},
  \code{\link{ci.sp}}, \code{\link{ci.se}} or \code{\link{ci.coords}} are called.
  
  When the ROC curve has an \code{\link{auc}} of 1 (or 100\%), the confidence interval will always be null
  (there is no interval). This is true for both \dQuote{delong} and \dQuote{bootstrap} methods that can
  not properly assess the variance in this case. This result is misleading, as the variance is of course not null.
  A \code{\link{warning}} will be displayed to inform of this condition, and of the misleading output.
}

\value{
The return value of the specific \code{ci} functions
\code{\link{ci.auc}}, \code{\link{ci.thresholds}}, \code{\link{ci.sp}}, \code{\link{ci.se}} or \code{\link{ci.coords}}, depending on the
\code{of} argument.
}

\references{
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \href{http://dx.doi.org/10.1186/1471-2105-12-77}{10.1186/1471-2105-12-77}.
}

\seealso{
  \code{\link{roc}}, \code{\link{auc}}, \code{\link{ci.auc}},
  \code{\link{ci.thresholds}}, \code{\link{ci.sp}}, \code{\link{ci.se}}, \code{\link{ci.coords}}
}

\examples{
data(aSAH)

# Syntax (response, predictor):
ci(aSAH$outcome, aSAH$s100b)

# With a roc object:
rocobj <- roc(aSAH$outcome, aSAH$s100b)

# Of an AUC 
ci(rocobj)
ci(rocobj, of="auc")
# this is strictly equivalent to:
ci.auc(rocobj)

# Of thresholds, sp, se...
\dontrun{
ci(rocobj, of="thresholds")
ci(rocobj, of="thresholds", thresholds=0.51)
ci(rocobj, of="thresholds", thresholds="all")
ci(rocobj, of="sp", sensitivities=c(.95, .9, .85))
ci(rocobj, of="se")
}

# Alternatively, you can get the CI directly from roc():
rocobj <- roc(aSAH$outcome, aSAH$s100b, ci=TRUE, of="auc")
rocobj$ci

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
