\encoding{UTF-8}
\name{ci.coords}
\alias{ci.coords}
\alias{ci.coords.default}
\alias{ci.coords.formula}
\alias{ci.coords.roc}
\alias{ci.coords.smooth.roc}

\title{
  Compute the confidence interval of arbitrary coordinates
}

\description{
  This function computes the confidence interval (CI) of the coordinates
  of a ROC curves with the \code{\link{coords}} function.
  By default, the 95\% CI are computed with 2000 stratified bootstrap replicates.
}

\usage{
# ci.coords(...)
\S3method{ci.coords}{roc}(roc, x,
input=c("threshold", "specificity", "sensitivity"),
ret=c("threshold", "specificity", "sensitivity"),
best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5),
conf.level=0.95, boot.n=2000,
boot.stratified=TRUE,
progress=getOption("pROCProgress")$name, ...) 
\S3method{ci.coords}{formula}(formula, data, ...)
\S3method{ci.coords}{smooth.roc}(smooth.roc, x,
input=c("specificity", "sensitivity"), ret=c("specificity", "sensitivity"),
best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5),
conf.level=0.95, boot.n=2000,
boot.stratified=TRUE,
progress=getOption("pROCProgress")$name, ...)
\S3method{ci.coords}{default}(response, predictor, ...)
}
		   
\arguments{
  \item{roc, smooth.roc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, or a \dQuote{smooth.roc} object from the
	\code{\link{smooth.roc}} function.
  }
  \item{response, predictor}{arguments for the \code{\link{roc}} function.}
  \item{formula, data}{a formula (and possibly a data object) of type
	response~predictor for the \code{\link{roc}} function.
  }
  \item{x, input, ret, best.method, best.weights}{Arguments passed to \code{\link{coords}}.
    See there for more details. The only difference is on the \code{x} argument which cannot be
    \dQuote{all} or \dQuote{local maximas}.
  }
  \item{conf.level}{the width of the confidence interval as [0,1], never
    in percent. Default: 0.95, resulting in a 95\% CI.
  }
  \item{boot.n}{the number of bootstrap replicates. Default: 2000.}
  \item{boot.stratified}{should the bootstrap be stratified (default, same number
	of cases/controls in each replicate than in the original sample) or
	not.
  }
  \item{progress}{the name of progress bar to display. Typically
    \dQuote{none}, \dQuote{win}, \dQuote{tk} or \dQuote{text} (see the
    \code{name} argument to \code{\link[plyr]{create_progress_bar}} for
    more information), but a list as returned by \code{\link[plyr]{create_progress_bar}}
    is also accepted. See also the \dQuote{Progress bars} section of
    \link[=pROC-package]{this package's documentation}.
  }
  \item{\dots}{further arguments passed to or from other methods,
    especially arguments for \code{\link{roc}} and \code{ci.coords.roc}
    when calling \code{ci.coords.default} or \code{ci.coords.formula}.
    Arguments for \code{\link{txtProgressBar}} (only 
    \code{char} and \code{style}) if applicable.
  }
}

\details{
  \code{ci.coords.formula} and \code{ci.coords.default} are convenience methods
  that build the ROC curve (with the \code{\link{roc}} function) before
  calling \code{ci.coords.roc}. You can pass them arguments for both
  \code{\link{roc}} and \code{ci.coords.roc}. Simply use \code{ci.coords}
  that will dispatch to the correct method.

  This function creates \code{boot.n} bootstrap replicate of the ROC
  curve, and evaluates the coordinates specified by the \code{x}, \code{input},
  \code{ret}, \code{best.method} and \code{best.weights} arguments. Then it computes the
  confidence interval as the percentiles given by \code{conf.level}.

  For more details about the bootstrap, see the Bootstrap section in
  \link[=pROC-package]{this package's documentation}.
}

\section{Warnings}{
  If \code{boot.stratified=FALSE} and the sample has a large imbalance between
  cases and controls, it could happen that one or more of the replicates
  contains no case or control observation, producing a \code{NA} area.
  The warning \dQuote{NA value(s) produced during bootstrap were ignored.}
  will be issued and the observation will be ignored. If you have a large
  imbalance in your sample, it could be safer to keep
  \code{boot.stratified=TRUE}.
}

\value{
  A matrix of class \dQuote{ci.coords}, with the confidence
  intervals of the CI. The matrix has 3 columns (lower bound, median and upper bound) and
  as many rows as confidence intervals computed. If more than one element was provided in \code{x},
  
  Additionally, the list has the following attributes: 
  \item{conf.level}{the width of the CI, in fraction.}
  \item{boot.n}{the number of bootstrap replicates.}
  \item{boot.stratified}{whether or not the bootstrapping was stratified.}
  \item{roc}{the object of class \dQuote{\link{roc}} that was used to
	compute the CI.
  }
}

\references{
  James Carpenter and John Bithell (2000) ``Bootstrap condence intervals:
  when, which, what? A practical guide for medical statisticians''.
  \emph{Statistics in Medicine} \bold{19}, 1141--1164.
  DOI: \href{http://dx.doi.org/10.1002/(SICI)1097-0258(20000515)19:9<1141::AID-SIM479>3.0.CO;2-F}{10.1002/(SICI)1097-0258(20000515)19:9<1141::AID-SIM479>3.0.CO;2-F}.

  Tom Fawcett (2006) ``An introduction to ROC analysis''. \emph{Pattern
    Recognition Letters} \bold{27}, 861--874. DOI:
  \href{http://dx.doi.org/10.1016/j.patrec.2005.10.010}{10.1016/j.patrec.2005.10.010}.
  
  Hadley Wickham (2011) ``The Split-Apply-Combine Strategy for Data Analysis''. \emph{Journal of Statistical Software}, \bold{40}, 1--29.
  URL: \href{http://www.jstatsoft.org/v40/i01}{www.jstatsoft.org/v40/i01}.
}

\seealso{
  \code{\link{roc}},
  \code{\link{coords}},
  \code{\link{ci}}
  
  CRAN package \pkg{plyr}, employed in this function.
}

\examples{
data(aSAH)

\dontrun{
# Syntax (response, predictor):
ci.coords(aSAH$outcome, aSAH$s100b, x="best", input = "threshold", ret=c("specificity", "ppv", "tp"))

# With a roc object:
rocobj <- roc(aSAH$outcome, aSAH$s100b)
ci.coords(rocobj, x=0.9, input = "sensitivity", ret=c("specificity", "ppv", "tp"))

# With a smoothed roc:
rocobj <- roc(aSAH$outcome, aSAH$s100b)
ci.coords(smooth(rocobj), x=0.9, input = "sensitivity", ret=c("specificity", "ppv", "tp"))

# Return everything we can:
rets <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", 
          "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
ci.coords(rocobj, x="best", input = "threshold", ret=rets)
}
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
