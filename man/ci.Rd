\encoding{UTF-8}
\name{ci}
\alias{ci}
\alias{ci.default}
\alias{ci.formula}
\alias{ci.roc}
\alias{ci.smooth.roc}
\alias{ci.multiclass.roc}
\alias{ci.multiclass.auc}

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
\S3method{ci}{multiclass.roc}(multiclass.roc, of = "auc", ...)
\S3method{ci}{multiclass.auc}(multiclass.auc, of = "auc", ...)
\S3method{ci}{formula}(formula, data, ...)
\S3method{ci}{default}(response, predictor, ...)
}
		   
\arguments{
  \item{roc, smooth.roc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, or a \dQuote{smooth.roc} object from the
	\code{\link[=smooth.roc]{smooth}} function.
  }
  \item{multiclass.roc, multiclass.auc}{not implemented.}
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
  
  CI of multiclass ROC curves and AUC is not implemented yet. Attempting to call these
  methods returns an error.
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
  DOI: \doi{10.1186/1471-2105-12-77}.
}

\seealso{
  \code{\link{roc}}, \code{\link{auc}}, \code{\link{ci.auc}},
  \code{\link{ci.thresholds}}, \code{\link{ci.sp}}, \code{\link{ci.se}}, \code{\link{ci.coords}}
}

\examples{

# Create a ROC curve:
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)


## AUC ## 
ci(roc1)
# this is equivalent to:
ci(roc1, of = "auc")
# or:
ci.auc(roc1)


## Coordinates ##
\dontrun{
# Thresholds
ci(roc1, of = "thresholds")
ci(roc1, of = "thresholds", thresholds = "all")
ci(roc1, of = "thresholds", thresholds = 0.51)
# equivalent to:
ci.thresholds(roc1, thresholds = 0.51)

# SE/SP
ci(roc1, of = "sp", sensitivities = c(.95, .9, .85))
ci.sp(roc1)
ci(roc1, of = "se")
ci.se(roc1)

# Arbitrary coordinates
ci(roc1, of = "coords", "best")
ci.coords(roc1, 0.51, "threshold")}

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
