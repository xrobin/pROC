\encoding{UTF-8}
\name{auc}
\alias{auc}
\alias{auc.default}
\alias{auc.formula}
\alias{auc.roc}
\alias{auc.smooth.roc}
\alias{auc.multiclass.roc}

\title{
  Compute the area under the ROC curve
}
\description{
  This function computes the numeric value of area under the ROC curve
  (AUC) with the trapezoidal rule. Two syntaxes are possible: one object of class \dQuote{\link{roc}}, or either
  two vectors (response, predictor) or a formula (response~predictor) as
  in the \code{\link{roc}} function.
  By default, the total AUC is computed, but a portion of the ROC curve
  can be specified with \code{partial.auc}.
}
\usage{
auc(...)
\S3method{auc}{roc}(roc, partial.auc=FALSE, partial.auc.focus=c("specificity",
"sensitivity"), partial.auc.correct=FALSE, 
allow.invalid.partial.auc.correct = FALSE, ...)
\S3method{auc}{smooth.roc}(smooth.roc, ...)
\S3method{auc}{multiclass.roc}(multiclass.roc, ...)
\S3method{auc}{formula}(formula, data, ...)
\S3method{auc}{default}(response, predictor, ...)
}

\arguments{
  \item{roc, smooth.roc, multiclass.roc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, a \dQuote{smooth.roc} object from the
	\code{\link[=smooth.roc]{smooth}} function, or a \dQuote{multiclass.roc}
	or \dQuote{mv.multiclass.roc} from the \code{\link{multiclass.roc}} function.
  }
  \item{response, predictor}{arguments for the \code{\link{roc}} function.}
  \item{formula, data}{a formula (and possibly a data object) of type response~predictor for the \code{\link{roc}} function.}
  \item{partial.auc}{either \code{FALSE} (default: consider total area) or a
    numeric vector of length 2: boundaries of the AUC to consider in
    [0,1] (or [0,100] if percent is \code{TRUE}).
  }
  \item{partial.auc.focus}{if \code{partial.auc} is not \code{FALSE} and a partial
    AUC is computed, specifies if \code{partial.auc} specifies the bounds in
    terms of specificity (default) or sensitivity. Can be shortened to spec/sens
    or even sp/se. Ignored if \code{partial.auc=FALSE}.
  }
  \item{partial.auc.correct}{logical indicating if the correction of
    AUC must be applied in order to have a maximal AUC of 1.0 and a
    non-discriminant AUC of 0.5 whatever the \code{partial.auc}
    defined. Ignored if \code{partial.auc=FALSE}. Default: \code{FALSE}.
  }
  \item{allow.invalid.partial.auc.correct}{logical indicating if
  the correction must return \code{NA} (with a \code{\link{warning}})
  when attempting to correct a pAUC below the diagonal.
  Set to \code{TRUE} to return a (probably invalid) corrected AUC.
  This is useful especially to avoid introducing a bias against low
  pAUCs in bootstrap operations.}
  \item{\dots}{further arguments passed to or from other methods,
    especially arguments for \code{\link{roc}} when calling
	\code{auc.default}, \code{auc.formula}, \code{auc.smooth.roc}.
	Note that the \code{auc}
	argument of \code{\link{roc}} is not allowed. Unused in
	\code{auc.roc}. 
  }
}

\details{
  This function is typically called from \code{\link{roc}} when \code{auc=TRUE}
  (default). It is also used by \code{\link{ci}}. When it is called with
  two vectors (response, predictor) or a formula (response~predictor)
  arguments, the \code{\link{roc}} function is called and only the AUC is
  returned.

  By default the total area under the curve is computed, but a partial AUC (pAUC)
  can be
  specified with the \code{partial.auc} argument. It specifies the bounds of
  specificity or sensitivity (depending on \code{partial.auc.focus}) between
  which the AUC will be computed. As it specifies specificities or
  sensitivities, you must adapt it in relation to the 'percent'
  specification (see details in \code{\link{roc}}).
  
  \code{partial.auc.focus} is ignored if
  \code{partial.auc=FALSE} (default). If a partial AUC is computed,
  \code{partial.auc.focus} specifies if the bounds specified in
  \code{partial.auc} must be interpreted as sensitivity or
  specificity. Any other value will produce an error. It is recommended to
  \code{\link[=plot.roc]{plot}} the ROC curve with \code{auc.polygon=TRUE} in order to 
  make sure the specification is correct.
  
  If a pAUC is defined, it can be standardized (corrected). This correction is 
  controled by the \code{partial.auc.correct} argument. If \code{partial.auc.correct=TRUE},
  the correction by McClish will be applied:
  
  \deqn{\frac{1+\frac{auc-min}{max-min}}{2}}{(1+(auc-min)/(max-min))/2}
  
  where auc is the uncorrected pAUC computed in the region defined by \code{partial.auc},
  min is the value of the non-discriminant AUC (with an AUC of 0.5 or 50%)
  in the region and max is the maximum possible AUC in the region. With this correction, the AUC
  will be 0.5 if non discriminant and 1.0 if maximal, whatever the region
  defined. This correction is fully compatible with \code{percent}.
  
  Note that this correction is undefined for curves below the diagonal (auc < min). Attempting
  to correct such an AUC will return \code{NA} with a warning.
}
\section{Smoothed ROC curves}{
  There is no difference in the computation of the area under a smoothed
  ROC curve, except for curves smoothed with \code{method="binomial"}. In this case
  and only if a full AUC is requested, the classical binormal AUC formula is applied:

  \deqn{auc=\phi\frac{a}{\sqrt{1 + b^2}}.}{pnorm(a/sqrt(1+b^2).}

  If the ROC curve is smoothed with any other \code{method} or if a partial AUC
  is requested, the empirical AUC described in the previous section is applied.
}
\section{Multi-class AUCs}{
  With an object of class \dQuote{multiclass.roc}, a multi-class AUC is
  computed as an average AUC as defined by Hand and Till (equation 7).

  \deqn{auc=\frac{2}{c(c-1)}\sum{aucs}}{2/(count * (count -
    1))*sum(aucs)}

  with aucs all the pairwise roc curves.

}

\value{
  The numeric AUC value, of class \code{c("auc", "numeric")} (or
  \code{c("multiclass.auc", "numeric")} or \code{c("mv.multiclass.auc", "numeric")}
  if a \dQuote{multiclass.roc} was supplied), in
  fraction of the area or in percent if \code{percent=TRUE}, with the
  following attributes:
  \item{partial.auc}{if the AUC is full (FALSE) or partial (and in this
    case the bounds), as defined in argument.}
  \item{partial.auc.focus}{only for a partial AUC, if the bound specifies
    the sensitivity or specificity, as defined in argument.}
  \item{partial.auc.correct}{only for a partial AUC, was it corrected?
    As defined in argument.}
  \item{percent}{whether the AUC is given in percent or fraction.}
  \item{roc}{the original ROC curve, as a \dQuote{\link{roc}},
    \dQuote{\link{smooth.roc}} or \dQuote{\link{multiclass.roc}} object.}
}
\references{
  Tom Fawcett (2006) ``An introduction to ROC analysis''. \emph{Pattern
    Recognition Letters} \bold{27}, 861--874. DOI:
  \doi{10.1016/j.patrec.2005.10.010}.

  David J. Hand and Robert J. Till (2001). A Simple Generalisation of
  the Area Under the ROC Curve for Multiple Class Classification
  Problems. \emph{Machine Learning} \bold{45}(2), p. 171--186. DOI:
  \doi{10.1023/A:1010920819831}.

  Donna Katzman McClish (1989) ``Analyzing a Portion of the ROC
  Curve''. \emph{Medical Decision Making} \bold{9}(3), 190--195. DOI:
  \doi{10.1177/0272989X8900900307}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
}

\seealso{
  \code{\link{roc}}, \code{\link{ci.auc}}
}
\examples{

# Create a ROC curve:
data(aSAH)
roc.s100b <- roc(aSAH$outcome, aSAH$s100b)

# Get the full AUC
auc(roc.s100b)

# Get the partial AUC:
auc(roc.s100b, partial.auc=c(1, .8), partial.auc.focus="se", partial.auc.correct=TRUE)

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
