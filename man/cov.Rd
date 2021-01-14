\encoding{UTF-8}
\name{cov.roc}
\alias{cov}
\alias{cov.default}
\alias{cov.auc}
\alias{cov.smooth.roc}
\alias{cov.roc}

\title{
  Covariance of two paired ROC curves
}
\description{
  This function computes the covariance between the AUC of two correlated (or paired) ROC
  curves. 
}
\usage{
cov(...)
\S3method{cov}{default}(...)
\S3method{cov}{auc}(roc1, roc2, ...)
\S3method{cov}{smooth.roc}(roc1, roc2, ...)
\S3method{cov}{roc}(roc1, roc2, method=c("delong", "bootstrap", "obuchowski"),
  reuse.auc=TRUE, boot.n=2000, boot.stratified=TRUE, boot.return=FALSE,
  progress=getOption("pROCProgress")$name, parallel=FALSE, ...)
}

\arguments{
  \item{roc1, roc2}{the two ROC curves on which to compute the covariance. Either
	\dQuote{\link{roc}}, \dQuote{\link{auc}} or
	\dQuote{\link{smooth.roc}} objects (types can be mixed as long as
	the original ROC curve are paired).
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
  \item{boot.return}{
    if \var{TRUE} and \code{method="bootstrap"}, also return the
    bootstrapped values. See the \dQuote{Value} section for
    more details.
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
    especially arguments for \code{cov.roc} when calling \code{cov},
    \code{cov.auc} or \code{cov.smooth.roc}. Arguments for
    \code{\link{auc}} (if \code{reuse.auc=FALSE}) and
    \code{\link{txtProgressBar}} (only \code{char} and \code{style}) if
    applicable.
  }
}

\details{
  This function computes the covariance between the AUC of two
  correlated (or paired, according to the detection of \code{\link{are.paired}}) ROC
  curves. It is typically called with the two \link{roc} objects of
  interest. Two methods are available: \dQuote{delong} and
  \dQuote{bootstrap} (see \dQuote{Computational
  details} section below).

  The default is to use \dQuote{delong} method except with
  partial AUC and smoothed curves where \dQuote{bootstrap} is employed.
  Using \dQuote{delong} for partial AUC and smoothed ROCs is not
  supported.

  For \link[=smooth.roc]{smoothed ROC curves}, smoothing is performed again at each
  bootstrap replicate with the parameters originally provided.
  If a density smoothing was performed with user-provided
  \code{density.cases} or \code{density.controls} the bootstrap cannot
  be performed and an error is issued.

  \code{cov.default} forces the usage of the
  \code{\link[stats:cor]{cov}} function in the \pkg{stats} package, so
  that other code relying on \code{cov} should continue to function
  normally.
}

\section{AUC specification}{
  To compute the covariance of the AUC of the ROC curves, \code{cov} needs a specification of the
  AUC. The specification is defined by:
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
    \item for each bootstrap replicate, the AUC of the two ROC curves
    are computed and stored.
    \item the variance (as per \code{\link{var.roc}}) of the resampled
    AUCs and their covariance are assessed in a single bootstrap pass.
    \item The following formula is used to compute the final covariance:
    \eqn{Var[AUC1] + Var[AUC2] - 2cov[AUC1,AUC2]}
  }

  With \code{method="delong"}, the processing is done as described in
  Hanley and Hajian-Tilaki (1997) using the algorithm by Sun and Xu (2014). 

  With \code{method="obuchowski"}, the processing is done as described
  in Obuchowski and McClish (1997), Table 1 and Equation 5, p. 1531. The
  computation of \eqn{g} for partial area under the ROC curve is
  modified as: 

  \deqn{expr1 * (2 * pi * expr2) ^ {(-1)} * (-expr4) - A * B * expr1 * (2 * pi * expr2^3) ^ {(-1/2)} * expr3}.
}

\section{Binormality assumption}{
  The \dQuote{obuchowski} method makes the assumption that the data is binormal.
  If the data shows a deviation from this assumption, it might help to
  normalize the data first (that is, before calling \code{\link{roc}}),
  for example with quantile normalization: 
  \preformatted{
    norm.x <- qnorm(rank(x)/(length(x)+1))
    cov(roc(response, norm.x, ...), ...)
  }

  \dQuote{delong} and \dQuote{bootstrap} methods make no such assumption.
}

\value{
  The numeric value of the covariance.

  If \code{boot.return=TRUE} and \code{method="bootstrap"}, an attribute
  \code{resampled.values} is set with the resampled (bootstrapped)
  values. It contains a matrix with the columns representing the two ROC
  curves, and the rows the \code{boot.n} bootstrap replicates.
}

\section{Errors}{
  If \code{density.cases} and \code{density.controls} were provided
  for smoothing, the error \dQuote{Cannot compute the covariance on ROC
    curves smoothed with density.controls and density.cases.} is
  issued.
}

\section{Warnings}{
  If \dQuote{auc} specifications are different in both roc objects, the
  warning \dQuote{Different AUC specifications in the ROC
    curves. Enforcing the inconsistency, but unexpected results may be
    produced.} is issued. Unexpected results may be produced.

  If one or both ROC curves are \dQuote{smooth.roc} objects with
  different smoothing specifications, the warning 
  \dQuote{Different smoothing parameters in the ROC curves. Enforcing
    the inconsistency, but unexpected results may be produced.} is issued.
  This warning can be benign, especially if ROC curves were generated
  with \code{roc(\ldots, smooth=TRUE)} with different arguments to other
  functions (such as plot), or if you really want to compare two ROC
  curves smoothed differently.

  If \code{method="delong"} and the AUC specification specifies a
  partial AUC, the warning \dQuote{Using DeLong for partial AUC is
    not supported. Using bootstrap test instead.} is issued. The
  \code{method} argument is ignored and \dQuote{bootstrap} is used instead.

  If \code{method="delong"} and the ROC
  curve is smoothed, the warning \dQuote{Using DeLong for
  smoothed ROCs is not supported. Using bootstrap instead.} is
  issued. The \code{method} argument is ignored and \dQuote{bootstrap}
  is used instead.

  DeLong ignores the direction of the ROC curve so that if two
  ROC curves have a different \code{direction}, the warning
  \dQuote{"DeLong should not be applied to ROC curves with a different
  direction."} is printed. However, the spurious computation is enforced.

  If \code{boot.stratified=FALSE} and the sample has a large imbalance between
  cases and controls, it could happen that one or more of the replicates
  contains no case or control observation, or that there are not enough
  points for smoothing, producing a \code{NA} area.
  The warning \dQuote{NA value(s) produced during bootstrap were ignored.}
  will be issued and the observation will be ignored. If you have a large
  imbalance in your sample, it could be safer to keep
  \code{boot.stratified=TRUE}.
  
  When both ROC curves have an \code{\link{auc}} of 1 (or 100\%), their covariance will always be null.
  This is true for both \dQuote{delong} and \dQuote{bootstrap} and methods. This result is misleading,
  as the covariance is of course not null.
  A \code{\link{warning}} will be displayed to inform of this condition, and of the misleading output.
}

\section{Messages}{
  The covariance can only be computed on paired data. This
  assumption is enforced by \code{\link{are.paired}}. If the ROC curves
  are not paired, the covariance is \code{0} and the message \dQuote{ROC
  curves are unpaired.} is printed. If your ROC curves are paired, make
  sure they fit \code{\link{are.paired}} criteria.
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
  \doi{10.1016/S1076-6332(97)80161-4}.

  Nancy A. Obuchowski, Donna K. McClish (1997). ``Sample size
  determination for diagnostic accurary studies involving binormal ROC
  curve indices''. \emph{Statistics in Medicine}, \bold{16}(13),
  1529--1542. DOI: \doi{10.1002/(SICI)1097-0258(19970715)16:13<1529::AID-SIM565>3.0.CO;2-H}.
  
  Xu Sun and Weichao Xu (2014) ``Fast Implementation of DeLongs Algorithm for Comparing
  the Areas Under Correlated Receiver Operating Characteristic Curves''. \emph{IEEE Signal
  Processing Letters}, \bold{21}, 1389--1393. 
  DOI: \doi{10.1109/LSP.2014.2337313}.
  
  Hadley Wickham (2011) ``The Split-Apply-Combine Strategy for Data Analysis''. \emph{Journal of Statistical Software}, \bold{40}, 1--29.
  URL: \href{https://www.jstatsoft.org/v40/i01}{www.jstatsoft.org/v40/i01}.
}

\seealso{
  \code{\link{roc}}, \code{\link{var.roc}}
  
  CRAN package \pkg{plyr}, employed in this function.
}

\examples{
data(aSAH)

# Basic example with 2 roc objects
roc1 <- roc(aSAH$outcome, aSAH$s100b)
roc2 <- roc(aSAH$outcome, aSAH$wfns)
cov(roc1, roc2)

\dontrun{
# The latter used Delong. To use bootstrap:
cov(roc1, roc2, method="bootstrap")
# Decrease boot.n for a faster execution:
cov(roc1, roc2, method="bootstrap", boot.n=1000)
}

# To use Obuchowski:
cov(roc1, roc2, method="obuchowski")

\dontrun{
# Comparison can be done on smoothed ROCs
# Smoothing is re-done at each iteration, and execution is slow
cov(smooth(roc1), smooth(roc2))
}
# or from an AUC (no smoothing)
cov(auc(roc1), roc2)

\dontrun{
# With bootstrap and return.values, one can compute the variances of the
# ROC curves in one single bootstrap run:
cov.rocs <- cov(roc1, roc2, method="bootstrap", boot.return=TRUE)
# var(roc1):
var(attr(cov.rocs, "resampled.values")[,1])
# var(roc2):
var(attr(cov.rocs, "resampled.values")[,2])
}

\dontrun{
# Covariance of partial AUC:
roc3 <- roc(aSAH$outcome, aSAH$s100b, partial.auc=c(1, 0.8), partial.auc.focus="se")
roc4 <- roc(aSAH$outcome, aSAH$wfns, partial.auc=c(1, 0.8), partial.auc.focus="se")
cov(roc3, roc4)
# This is strictly equivalent to:
cov(roc3, roc4, method="bootstrap")

# Alternatively, we could re-use roc1 and roc2 to get the same result:
cov(roc1, roc2, reuse.auc=FALSE, partial.auc=c(1, 0.8), partial.auc.focus="se")
}

# Spurious use of DeLong's test with different direction:
roc5 <- roc(aSAH$outcome, aSAH$s100b, direction="<")
roc6 <- roc(aSAH$outcome, aSAH$s100b, direction=">")
cov(roc5, roc6, method="delong")

## Test data from Hanley and Hajian-Tilaki, 1997
disease.present <- c("Yes", "No", "Yes", "No", "No", "Yes", "Yes", "No",
                     "No", "Yes", "No", "No", "Yes", "No", "No")
field.strength.1 <- c(1, 2, 5, 1, 1, 1, 2, 1, 2, 2, 1, 1, 5, 1, 1)
field.strength.2 <- c(1, 1, 5, 1, 1, 1, 4, 1, 2, 2, 1, 1, 5, 1, 1)
roc7 <- roc(disease.present, field.strength.1)
roc8 <- roc(disease.present, field.strength.2)
# Assess the covariance:
cov(roc7, roc8)

\dontrun{
# With bootstrap:
cov(roc7, roc8, method="bootstrap")
}

}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
