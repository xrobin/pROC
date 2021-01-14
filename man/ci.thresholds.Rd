\encoding{UTF-8}
\name{ci.thresholds}
\alias{ci.thresholds}
\alias{ci.thresholds.default}
\alias{ci.thresholds.formula}
\alias{ci.thresholds.roc}
\alias{ci.thresholds.smooth.roc}

\title{
  Compute the confidence interval of thresholds
}

\description{
  This function computes the confidence interval (CI) of the sensitivity
  and specificity of the thresholds given in argument.
  By default, the 95\% CI are computed with 2000 stratified bootstrap replicates.
}

\usage{
# ci.thresholds(...)
\S3method{ci.thresholds}{roc}(roc, conf.level=0.95, boot.n=2000,
boot.stratified=TRUE, thresholds = "local maximas",
progress=getOption("pROCProgress")$name, parallel=FALSE, ...) 
\S3method{ci.thresholds}{formula}(formula, data, ...)
\S3method{ci.thresholds}{smooth.roc}(smooth.roc, ...)
\S3method{ci.thresholds}{default}(response, predictor, ...)
}
		   
\arguments{
  \item{roc}{a \dQuote{roc} object from the \code{\link{roc}} function.}
  \item{smooth.roc}{not available for \link[=smooth.roc]{smoothed} ROC
    curves, available only to catch the error and provide a clear error
    message.
  }
  \item{response, predictor}{arguments for the \code{\link{roc}} function.}
  \item{formula, data}{a formula (and possibly a data object) of type
	response~predictor for the \code{\link{roc}} function.
  }
  \item{conf.level}{the width of the confidence interval as [0,1], never
  in percent. Default: 0.95, resulting in a 95\% CI.
  }
  \item{boot.n}{the number of bootstrap replicates. Default: 2000.}
  \item{boot.stratified}{should the bootstrap be stratified (default, same number
	of cases/controls in each replicate than in the original sample) or
	not.
  }
  \item{thresholds}{on which thresholds to evaluate the CI. Either the
	numeric values of the thresholds, a logical vector (as index of
	\code{roc$thresholds}) or a character \dQuote{all}, \dQuote{local
	  maximas} or \dQuote{best} that will be used to determine the threshold(s) 
	  on the supplied curve with \code{\link{coords}} (not on the resampled curves).
  }
  \item{progress}{the name of progress bar to display. Typically
    \dQuote{none}, \dQuote{win}, \dQuote{tk} or \dQuote{text} (see the
    \code{name} argument to \code{\link[plyr]{create_progress_bar}} for
    more information), but a list as returned by \code{\link[plyr]{create_progress_bar}}
    is also accepted. See also the \dQuote{Progress bars} section of
    \link[=pROC-package]{this package's documentation}.
  }
  \item{parallel}{if TRUE, the bootstrap is processed in parallel, using 
    parallel backend provided by plyr (foreach).
  }
  \item{\dots}{further arguments passed to or from other methods,
    especially arguments for \code{\link{roc}} and \code{ci.thresholds.roc}
    when calling \code{ci.thresholds.default} or \code{ci.thresholds.formula}.
    Arguments for \code{\link{txtProgressBar}} (only 
    \code{char} and \code{style}) if applicable. 
    Arguments \code{best.method} and \code{best.weights} to \code{\link{coords}}.
  }
}

\details{
  \code{ci.thresholds.formula} and \code{ci.thresholds.default} are convenience methods
  that build the ROC curve (with the \code{\link{roc}} function) before
  calling \code{ci.thresholds.roc}. You can pass them arguments for both
  \code{\link{roc}} and \code{ci.thresholds.roc}. Simply use \code{ci.thresholds}
  that will dispatch to the correct method.

  This function creates \code{boot.n} bootstrap replicate of the ROC
  curve, and evaluates the sensitivity and specificity at thresholds
  given by the \code{thresholds} argument. Then it computes the
  confidence interval as the percentiles given by \code{conf.level}.
  
  A threshold given as a \code{logical} vector or \code{character} is converted to the corresponding numeric vector once
  \emph{using the supplied ROC curve}, and not at each bootstrap iteration. See \code{\link{ci.coords}} for the latter behaviour.

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
  A list of length 2 and class \dQuote{ci.thresholds}, \dQuote{ci} and \dQuote{list} (in this order), with the confidence
  intervals of the CI and the following items:
  \item{specificity}{a matrix of CI for the specificity. Row (names) are the
	thresholds, the first column the lower bound, the 2nd column the
	median and the 3rd column the upper bound.
  }
  \item{sensitivity}{same than specificity.}
  Additionally, the list has the following attributes: 
  \item{conf.level}{the width of the CI, in fraction.}
  \item{boot.n}{the number of bootstrap replicates.}
  \item{boot.stratified}{whether or not the bootstrapping was stratified.}
  \item{thresholds}{the thresholds, as given in argument.}
  \item{roc}{the object of class \dQuote{\link{roc}} that was used to
	compute the CI.
  }
}

\references{
  James Carpenter and John Bithell (2000) ``Bootstrap condence intervals:
  when, which, what? A practical guide for medical statisticians''.
  \emph{Statistics in Medicine} \bold{19}, 1141--1164.
  DOI: \doi{10.1002/(SICI)1097-0258(20000515)19:9<1141::AID-SIM479>3.0.CO;2-F}.

  Tom Fawcett (2006) ``An introduction to ROC analysis''. \emph{Pattern
    Recognition Letters} \bold{27}, 861--874. DOI:
  \doi{10.1016/j.patrec.2005.10.010}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
  
  Hadley Wickham (2011) ``The Split-Apply-Combine Strategy for Data Analysis''. \emph{Journal of Statistical Software}, \bold{40}, 1--29.
  URL: \href{https://www.jstatsoft.org/v40/i01}{www.jstatsoft.org/v40/i01}.
}

\seealso{
  \code{\link{roc}},
  \code{\link{ci}}
}

\examples{
data(aSAH)

# Create a ROC curve:
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)


## Basic example ##
# Compute CI of all local maxima thresholds
\dontrun{
ci.thresholds(roc1)}\dontshow{ci.thresholds(roc1, boot.n = 10)}


## More options ##
# Customized bootstrap and thresholds:
\dontrun{
ci.thresholds(roc1,
			  thresholds=c(0.5, 1, 2),
              boot.n=10000, conf.level=0.9, stratified=FALSE)}\dontshow{
ci.thresholds(roc1,
			  thresholds=c(0.5, 1, 2),
              boot.n=10, conf.level=0.9, stratified=FALSE)}


## Plotting the CI ##
\dontrun{
ci1 <- ci.thresholds(roc1)}\dontshow{
ci1 <- ci.thresholds(roc1, boot.n = 10)}
plot(roc1)
plot(ci1)
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
