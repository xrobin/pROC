\encoding{UTF-8}
\name{ci.auc}
\alias{ci.auc}
\alias{ci.auc.auc}
\alias{ci.auc.default}
\alias{ci.auc.formula}
\alias{ci.auc.roc}
\alias{ci.auc.smooth.roc}
\alias{ci.auc.multiclass.roc}
\alias{ci.auc.multiclass.auc}

\title{
  Compute the confidence interval of the AUC
}

\description{
  This function computes the confidence interval (CI) of an area under
  the curve (AUC). By default, the 95\% CI is computed with 2000
  stratified bootstrap replicates.
}

\usage{
# ci.auc(...)
\S3method{ci.auc}{roc}(roc, conf.level=0.95, method=c("delong",
"bootstrap"), boot.n = 2000, boot.stratified = TRUE, reuse.auc=TRUE,
progress = getOption("pROCProgress")$name, parallel=FALSE, ...)
\S3method{ci.auc}{smooth.roc}(smooth.roc, conf.level=0.95, boot.n=2000,
boot.stratified=TRUE, reuse.auc=TRUE,
progress=getOption("pROCProgress")$name, parallel=FALSE, ...)
\S3method{ci.auc}{auc}(auc, ...)
\S3method{ci.auc}{multiclass.roc}(multiclass.roc, ...)
\S3method{ci.auc}{multiclass.auc}(multiclass.auc, ...)
\S3method{ci.auc}{auc}(auc, ...)
\S3method{ci.auc}{formula}(formula, data, ...)
\S3method{ci.auc}{default}(response, predictor, ...)
}
		   
\arguments{
  \item{roc, smooth.roc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, or a \dQuote{smooth.roc} object from the
	\code{\link[=smooth.roc]{smooth}} function.
  }
  \item{auc}{an \dQuote{auc} object from the \code{\link{auc}} function.}
  \item{multiclass.roc, multiclass.auc}{not implemented.}
  \item{response, predictor}{arguments for the \code{\link{roc}} function.}
  \item{formula, data}{a formula (and possibly a data object) of type
	response~predictor for the \code{\link{roc}} function.
  }
  \item{conf.level}{the width of the confidence interval as [0,1], never
  in percent. Default: 0.95, resulting in a 95\% CI.
  }
  \item{method}{the method to use, either \dQuote{delong} or
    \dQuote{bootstrap}. The first letter is sufficient. If omitted, the
    appropriate method is selected as explained in details.
  }
  \item{boot.n}{the number of bootstrap replicates. Default: 2000.}
  \item{boot.stratified}{should the bootstrap be stratified (default, same number
	of cases/controls in each replicate than in the original sample) or
	not.
  }
  \item{reuse.auc}{if \code{TRUE} (default) and the \dQuote{roc} object
    contains an \dQuote{auc} field, re-use these specifications for the
    test. If false, use optional \code{\dots} arguments to
	\code{\link{auc}}. See details.
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
    especially arguments for \code{\link{roc}} and \code{roc.test.roc}
    when calling \code{roc.test.default} or \code{roc.test.formula}.
    Arguments for \code{\link{auc}}
    and \code{\link{txtProgressBar}} (only \code{char} and \code{style})
    if applicable.
  }
}

\details{
  This function computes the CI of an AUC. Two methods are available:
  \dQuote{delong} and \dQuote{bootstrap} with the parameters defined in \dQuote{roc$auc} to
  compute a CI. When it is called with two vectors (response, predictor)
  or a formula (response~predictor) arguments, the \code{\link{roc}}
  function is called to build the ROC curve first.

  The default is to use
  \dQuote{delong} method except for comparison of partial AUC and smoothed
  curves, where \code{bootstrap} is used. Using \dQuote{delong} for
  partial AUC and smoothed ROCs is not supported.

  With \code{method="bootstrap"}, the function calls \code{\link{auc}}
  \code{boot.n} times. For more details about the bootstrap, see the Bootstrap section in
  \link[=pROC-package]{this package's documentation}.

  For \link[=smooth.roc]{smoothed ROC curves}, smoothing is performed again at each
  bootstrap replicate with the parameters originally provided.
  If a density smoothing was performed with user-provided
  \code{density.cases} or \code{density.controls} the bootstrap cannot
  be performed and an error is issued.

  With \code{method="delong"}, the variance of the AUC is computed as
  defined by DeLong \emph{et al.} (1988) using the algorithm by Sun and Xu (2014) 
  and the CI is deduced with \code{\link{qnorm}}.
  
  CI of multiclass ROC curves and AUC is not implemented yet. Attempting to call these
  methods returns an error.
}

\section{AUC specification}{
  The comparison of the CI needs a specification of the AUC. This allows
  to compute the CI for full or partial AUCs. The specification is defined by:
  \enumerate{
    \item the \dQuote{auc} field in the \dQuote{\link{roc}} object if
    \code{reuse.auc} is set to \code{TRUE} (default). It is naturally
    inherited from any call to \code{\link{roc}} and fits most cases.
    \item passing the specification to \code{\link{auc}} with \dots
    (arguments \code{partial.auc}, \code{partial.auc.correct} and
    \code{partial.auc.focus}). In this case, you must ensure either that
    the \code{\link{roc}} object do not contain an \code{auc} field (if
    you called \code{\link{roc}} with \code{auc=FALSE}), or set
    \code{reuse.auc=FALSE}.
  }

  If \code{reuse.auc=FALSE} the \code{\link{auc}} function will always
  be called with \code{\dots} to determine the specification, even if
  the \dQuote{\link{roc}} object do contain an \code{auc} field.

  As well if the \dQuote{\link{roc}} object do not contain an \code{auc}
  field, the \code{\link{auc}} function will always be called with
  \code{\dots} to determine the specification.

  Warning: if the roc object passed to ci contains an \code{auc}
  field and \code{reuse.auc=TRUE}, \link{auc} is not called and
  arguments such as \code{partial.auc} are silently ignored.
  
}

\section{Warnings}{
  If \code{method="delong"} and the AUC specification specifies a
  partial AUC, the warning \dQuote{Using DeLong's test for partial AUC is
    not supported. Using bootstrap test instead.} is issued. The
  \code{method} argument is ignored and \dQuote{bootstrap} is used
  instead.

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
A numeric vector of length 3 and class \dQuote{ci.auc}, \dQuote{ci} and \dQuote{numeric} (in this order), with the lower
bound, the median and the upper bound of the CI, and the following attributes:
\item{conf.level}{the width of the CI, in fraction.}
\item{method}{the method employed.}
\item{boot.n}{the number of bootstrap replicates.}
\item{boot.stratified}{whether or not the bootstrapping was stratified.}
\item{auc}{an object of class \dQuote{\link{auc}} stored for reference about the
  compued AUC details (partial, percent, ...)}

The \code{aucs} item is not included in this list since version 1.2 for
consistency reasons.
}
\references{
  James Carpenter and John Bithell (2000) ``Bootstrap condence intervals:
  when, which, what? A practical guide for medical statisticians''.
  \emph{Statistics in Medicine} \bold{19}, 1141--1164.
  DOI: \doi{10.1002/(SICI)1097-0258(20000515)19:9<1141::AID-SIM479>3.0.CO;2-F}.

  Elisabeth R. DeLong, David M. DeLong and Daniel L. Clarke-Pearson
  (1988) ``Comparing the areas under two or more correlated receiver
  operating characteristic curves: a nonparametric
  approach''. \emph{Biometrics} \bold{44}, 837--845.
  
  Xu Sun and Weichao Xu (2014) ``Fast Implementation of DeLongs Algorithm for Comparing
  the Areas Under Correlated Receiver Operating Characteristic Curves''. \emph{IEEE Signal
  Processing Letters}, \bold{21}, 1389--1393. 
  DOI: \doi{10.1109/LSP.2014.2337313}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
  
  Hadley Wickham (2011) ``The Split-Apply-Combine Strategy for Data Analysis''. \emph{Journal of Statistical Software}, \bold{40}, 1--29.
  URL: \href{https://www.jstatsoft.org/v40/i01}{www.jstatsoft.org/v40/i01}.
}

\seealso{
  \code{\link{roc}}, \code{\link{auc}}, \code{\link{ci}}
}

\examples{
# Create a ROC curve:
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)


## Basic example ##
ci.auc(roc1)

# You can also write:
ci(roc1)
ci(auc(roc1))


## More options ##
# Partial AUC and customized bootstrap:
\dontrun{
ci.auc(roc1,
       conf.level=0.9,
	   partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10000, stratified=FALSE)}\dontshow{
ci.auc(roc1,
       conf.level=0.9,
	   partial.auc=c(1, .8),  partial.auc.focus="se", partial.auc.correct=TRUE,
       boot.n=10, stratified=FALSE)}

# Note that the following will NOT give a CI of the partial AUC:
\dontrun{
ci.auc(roc1,
       partial.auc=c(1, .8), partial.auc.focus="se", partial.auc.correct=FALSE)}
# This is because rocobj$auc is not a partial AUC and reuse.auc = TRUE by default.
# You can overcome this problem by passing an AUC instead:
auc1 <- auc(roc1, partial.auc=c(1, .8), partial.auc.focus="se", 
            partial.auc.correct=FALSE)
\dontrun{
ci.auc(auc1)}\dontshow{ci.auc(auc1, boot.n = 10)}


## On smoothed ROC curves with bootstrap ##
\dontrun{
ci.auc(smooth(roc1, method="density"))}\dontshow{
ci.auc(smooth(roc1, method="density"), boot.n = 10)}
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
