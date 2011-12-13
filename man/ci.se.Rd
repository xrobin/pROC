\encoding{UTF-8}
\name{ci.se}
\alias{ci.se}
\alias{ci.se.default}
\alias{ci.se.formula}
\alias{ci.se.roc}
\alias{ci.se.smooth.roc}

\title{
  Compute the confidence interval of sensitivities at given specificities
}

\description{
  This function computes the confidence interval (CI) of the sensitivity
  at the given specificity points.
  By default, the 95\% CI are computed with 2000 stratified bootstrap replicates.
}

\usage{
# ci.se(...)
\S3method{ci.se}{roc}(roc, specificities = seq(0, 1, .1) * ifelse(roc$percent,
100, 1), conf.level=0.95, boot.n=2000, boot.stratified=TRUE,
progress=getOption("pROCProgress")$name, ...) 
\S3method{ci.se}{smooth.roc}(smooth.roc, specificities = seq(0, 1, .1) *
ifelse(smooth.roc$percent, 100, 1), conf.level=0.95, boot.n=2000,
boot.stratified=TRUE, progress=getOption("pROCProgress")$name, ...)
\S3method{ci.se}{formula}(formula, data, ...)
\S3method{ci.se}{default}(response, predictor, ...)
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
  \item{specificities}{on which specificities to evaluate the CI.
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
    especially arguments for \code{\link{roc}} and \code{ci.se.roc}
    when calling \code{ci.se.default} or \code{ci.se.formula}.
    Arguments for \code{\link{txtProgressBar}} (only 
    \code{char} and \code{style}) if applicable.
  }
}

\details{
  \code{ci.se.formula} and \code{ci.se.default} are convenience methods
  that build the ROC curve (with the \code{\link{roc}} function) before
  calling \code{ci.se.roc}. You can pass them arguments for both
  \code{\link{roc}} and \code{ci.se.roc}. Simply use \code{ci.se}
  that will dispatch to the correct method.
  
  The \code{ci.se.roc} function creates \code{boot.n} bootstrap replicate of the ROC
  curve, and evaluates the sensitivity at specificities
  given by the \code{specificities} argument. Then it computes the
  confidence interval as the percentiles given by \code{conf.level}.

  For more details about the bootstrap, see the Bootstrap section in
  \link[=pROC-package]{this package's documentation}.

  For \link[=smooth.roc]{smoothed ROC curves}, smoothing is performed again at each
  bootstrap replicate with the parameters originally provided.
  If a density smoothing was performed with user-provided
  \code{density.cases} or \code{density.controls} the bootstrap cannot
  be performed and an error is issued.
}

\section{Warnings}{
  If \code{boot.stratified=FALSE} and the sample has a large imbalance between
  cases and controls, it could happen that one or more of the replicates
  contains no case or control observation, or that there are not enough
  points for smoothing, producing a \code{NA} area.
  The warning \dQuote{NA value(s) produced during bootstrap were ignored.}
  will be issued and the observation will be ignored. If you have a large
  imbalance in your sample, it could be safer to keep
  \code{boot.stratified=TRUE}.
}

\section{Errors}{
  If \code{density.cases} and \code{density.controls} were provided
  for smoothing, the error \dQuote{Cannot compute the statistic on ROC
  curves smoothed with density.controls and density.cases.} is issued.
}

\value{
  A matrix of CI for the given sensitivities. Row (names) are the
  specificities, the first column the lower bound, the 2nd column the
  median and the 3rd column the upper bound.

  Additionally, the list has the following attributes: 
  \item{conf.level}{the width of the CI, in fraction.}
  \item{boot.n}{the number of bootstrap replicates.}
  \item{boot.stratified}{whether or not the bootstrapping was stratified.}
  \item{specificities}{the specificities as given in argument.}
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
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \href{http://dx.doi.org/10.1186/1471-2105-12-77}{10.1186/1471-2105-12-77}.
}

\seealso{
  \code{\link{roc}},
  \code{\link{ci}},
  \code{\link{ci.sp}},
  \code{\link{plot.ci}}
}

\examples{
data(aSAH)

\dontrun{
# Syntax (response, predictor):
ci.se(aSAH$outcome, aSAH$s100b)

# With a roc object and less bootstrap:
rocobj <- roc(aSAH$outcome, aSAH$s100b)
ci.se(rocobj, boot.n=100)

# Customized bootstrap and specific specificities:
ci.se(rocobj, c(.95, .9, .85), boot.n=500, conf.level=0.9, stratified=FALSE)
}

# Alternatively, you can get the CI directly from roc():
rocobj <- roc(aSAH$outcome,
              aSAH$s100b, ci=TRUE, of="se", boot.n=100)
rocobj$ci

# Plotting the CI
plot(rocobj)
plot(rocobj$ci)

\dontrun{
# On a smoothed ROC, the CI is re-computed automatically
smooth(rocobj)
# Or you can compute a new one:
ci.se(smooth(rocobj, method="density", reuse.ci=FALSE), boot.n=100)
}
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
