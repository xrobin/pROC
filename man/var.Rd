\encoding{UTF-8}
\name{var.roc}
\alias{var}
\alias{var.roc}
\alias{var.auc}
\alias{var.smooth.roc}
\alias{var.default}
\title{
 Variance of a ROC curve
}
\description{
 These functions compute the variance of the AUC of a ROC curve.
}
\usage{
var(...)
\S3method{var}{default}(...)
\S3method{var}{auc}(auc, ...)
\S3method{var}{roc}(roc, method=c("delong", "bootstrap"), 
boot.n = 2000, boot.stratified = TRUE, reuse.auc=TRUE, 
progress = getOption("pROCProgress")$name, ...)
\S3method{var}{smooth.roc}(smooth.roc, ...)
}

\arguments{
  \item{roc, smooth.roc, auc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, a \dQuote{smooth.roc} object from the
	\code{\link{smooth.roc}} function or an \dQuote{auc} object from
	the \code{\link{auc}} function.
  }
  \item{method}{the method to use, either \dQuote{delong} or
    \dQuote{bootstrap}. The first letter is
    sufficient. If omitted, the appropriate method is selected as
    explained in details.
  }
  \item{reuse.auc}{if \code{TRUE} (default) and the \dQuote{roc} objects
    contain an \dQuote{auc} field, re-use these specifications for the
    test. See details.
  }
  \item{boot.n}{for \code{method="bootstrap"} only: the number of
    bootstrap replicates or permutations. Default: \var{2000}.
  }
  \item{boot.stratified}{for \code{method="bootstrap"} only:
    should the bootstrap be stratified (same number 
    of cases/controls in each replicate than in the original sample) or
    not. Default: \var{TRUE}.
  }
  \item{progress}{the name of progress bar to display. Typically
    \dQuote{none}, \dQuote{win}, \dQuote{tk} or \dQuote{text} (see the
    \code{name} argument to \code{\link[plyr]{create_progress_bar}} for
    more information), but a list as returned by \code{\link[plyr]{create_progress_bar}}
    is also accepted. See also the \dQuote{Progress bars} section of
    \link[=pROC-package]{this package's documentation}.
  }
  \item{\dots}{further arguments passed to or from other methods,
    especially arguments for \code{var.roc} when calling \code{var},
    \code{var.auc} and \code{var.smooth.roc}. Arguments for
    \code{\link{auc}} (if \code{reuse.auc=FALSE}) and
    \code{\link{txtProgressBar}} (only \code{char} and \code{style}) if
    applicable.
  }
}

\details{
  The \code{var} function computes the variance of the AUC of a ROC
  curve. It is typically called with the \link{roc} object of
  interest. Two methods are available: \dQuote{delong} and
  \dQuote{bootstrap} (see \dQuote{Computational
  details} section below).

  The default is to use \dQuote{delong} method except for with
  partial AUC and smoothed curves where \dQuote{bootstrap} is employed.
  Using \dQuote{delong} for partial AUC and smoothed ROCs is not
  supported (a warning is produced and \dQuote{bootstrap} is employed
  instead).

  For \link[=smooth.roc]{smoothed ROC curves}, smoothing is performed again at each
  bootstrap replicate with the parameters originally provided.
  If a density smoothing was performed with user-provided
  \code{density.cases} or \code{density.controls} the bootstrap cannot
  be performed and an error is issued.

  \code{var.default} forces the usage of the
  \code{\link[stats]{var}} function in the \pkg{stats} package, so
  that other code relying on \code{var} should continue to function
  normally.
}

\section{AUC specification}{
  \code{var} needs a specification of the AUC to compute
  the variance of the AUC of the ROC curve.
  The specification is defined by:
  \enumerate{
    \item the \dQuote{auc} field in the \dQuote{\link{roc}} objects if
    \code{reuse.auc} is set to \code{TRUE} (default)
    \item passing the specification to \code{\link{auc}} with \dots
    (arguments \code{partial.auc}, \code{partial.auc.correct} and
    \code{partial.auc.focus}). In this case, you must ensure either that
    the \code{\link{roc}} object do not contain an \code{auc} field (if
    you called \code{\link{roc}} with \code{auc=FALSE}), or set
    \code{reuse.auc=FALSE}.
  }

  If \code{reuse.auc=FALSE} the \code{\link{auc}} function will always
  be called with \code{\dots} to determine the specification, even if
  the \dQuote{\link{roc}} objects do contain an \code{auc} field.

  As well if the \dQuote{\link{roc}} objects do not contain an \code{auc}
  field, the \code{\link{auc}} function will always be called with
  \code{\dots} to determine the specification.

  Warning: if the roc object passed to roc.test contains an \code{auc}
  field and \code{reuse.auc=TRUE}, \link{auc} is not called and
  arguments such as \code{partial.auc} are silently ignored.
}


\section{Computation details}{
  With \code{method="bootstrap"}, the processing is done as follow:
  \enumerate{
    \item \code{boot.n} bootstrap replicates are drawn from the
    data. If \code{boot.stratified} is \var{TRUE}, each replicate contains
    exactly the same number of controls and cases than the original
    sample, otherwise if \var{FALSE} the numbers can vary.
    \item for each bootstrap replicate, the AUC of the ROC curve
    is computed and stored.
    \item the variance of the resampled AUCs are computed and returned.
  }

  With \code{method="delong"}, the processing is done as described in
  Hanley and Hajian-Tilaki (1997).
}

\value{
  The numeric value of the variance.
}

\section{Warnings}{  
  If \code{method="delong"} and the AUC specification specifies a
  partial AUC, the warning \dQuote{Using DeLong for partial AUC is
    not supported. Using bootstrap test instead.} is issued. The
  \code{method} argument is ignored and \dQuote{bootstrap} is used instead.

  If \code{method="delong"} and the ROC
  curve is smoothed, the warning \dQuote{Using DeLong for
  smoothed ROCs is not supported. Using bootstrap test instead.} is
  issued. The \code{method} argument is ignored and \dQuote{bootstrap}
  is used instead.

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
  for smoothing, the error \dQuote{Cannot compute the covariance on ROC
    curves smoothed with density.controls and density.cases.} is
  issued.
}

\references{
  Elisabeth R. DeLong, David M. DeLong and Daniel L. Clarke-Pearson
  (1988) ``Comparing the areas under two or more correlated receiver
  operating characteristic curves: a nonparametric
  approach''. \emph{Biometrics} \bold{44}, 837--845.

  James A. Hanley and Karim O. Hajian-Tilaki (1997) ``Sampling
  variability of nonparametric estimates of the areas under receiver
  operating characteristic curves: An update''. \emph{Academic
  Radiology} \bold{4}, 49--58. DOI:
  \href{http://dx.doi.org/10.1016/S1076-6332(97)80161-4}{10.1016/S1076-6332(97)80161-4}.
}

\seealso{
 \code{\link{roc}}, \code{\link{cov.roc}}
}

\examples{
data(aSAH)

##  Basic example
roc1 <- roc(aSAH$outcome, aSAH$s100b)
roc2 <- roc(aSAH$outcome, aSAH$wfns)
var(roc1)
var(roc2)

# We could also write it in one line:
var(roc(aSAH$outcome, aSAH$s100b))

\dontrun{
# The latter used Delong. To use bootstrap:
var(roc1, method="bootstrap")
# Decrease boot.n for a faster execution
var(roc1,method="bootstrap", boot.n=1000)
}

\dontrun{
# Variance of smoothed ROCs:
# Smoothing is re-done at each iteration, and execution is slow
var(smooth(roc1))
}

# or from an AUC (no smoothing)
var(auc(roc1))

## Test data from Hanley and Hajian-Tilaki, 1997
disease.present <- c("Yes", "No", "Yes", "No", "No", "Yes", "Yes", "No", "No", "Yes", "No", "No", "Yes", "No", "No")
field.strength.1 <- c(1, 2, 5, 1, 1, 1, 2, 1, 2, 2, 1, 1, 5, 1, 1)
field.strength.2 <- c(1, 1, 5, 1, 1, 1, 4, 1, 2, 2, 1, 1, 5, 1, 1)
roc3 <- roc(disease.present, field.strength.1)
roc4 <- roc(disease.present, field.strength.2)
# Assess the variance:
var(roc3)
var(roc4)

\dontrun{
# With bootstrap:
var(roc3, method="bootstrap")
var(roc4, method="bootstrap")
}

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
