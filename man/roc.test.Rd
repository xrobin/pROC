\encoding{UTF-8}
\name{roc.test}
\alias{roc.test}
\alias{roc.test.default}
\alias{roc.test.roc}
\alias{roc.test.formula}
\alias{roc.test.auc}
\alias{roc.test.smooth.roc}

\title{
  Compare the AUC of two ROC curves
}
\description{
  This function compares the AUC or partial AUC of two correlated (or
  paired) or uncorrelated (unpaired) ROC curves. Several syntaxes are
  available: two object of class roc (which can be AUC or smoothed ROC),
  or either three vectors (response, predictor1, predictor2) or a
  response vector and a matrix or data.frame with two columns
  (predictors).
}
\usage{
# roc.test(...)
\S3method{roc.test}{roc}(roc1, roc2, method=c("delong", "bootstrap",
"venkatraman", "sensitivity", "specificity"), sensitivity = NULL,
specificity = NULL, alternative = c("two.sided", "less", "greater"),
paired=NULL, reuse.auc=TRUE, boot.n=2000, boot.stratified=TRUE,
ties.method="first", progress=getOption("pROCProgress")$name,
parallel=FALSE, ...)
\S3method{roc.test}{auc}(roc1, roc2, ...)
\S3method{roc.test}{smooth.roc}(roc1, roc2, ...)
\S3method{roc.test}{formula}(formula, data, ...)
\S3method{roc.test}{default}(response, predictor1, predictor2=NULL,
na.rm=TRUE, method=NULL, ...)
}

\arguments{
  \item{roc1, roc2}{the two ROC curves to compare. Either
	\dQuote{\link{roc}}, \dQuote{\link{auc}} or
	\dQuote{\link{smooth.roc}} objects (types can be mixed).
  }
  \item{response}{a vector or factor, as for the \link{roc} function.}
  \item{predictor1}{a numeric or ordered vector as for the \link{roc} function, or
  a matrix or data.frame with predictors two colums.}
  \item{predictor2}{only if predictor1 was a vector, the second
    predictor as a numeric vector.
  }
  \item{formula}{a formula of the type response~predictor1+predictor2.
  	Additional arguments  \code{data}, \code{subset} and \code{na.action}
	 are supported, see \code{\link{model.frame}} for more details.
  }
  \item{data}{a matrix or data.frame containing the variables in the
    formula. See \code{\link{model.frame}} for more details.}
  \item{na.rm}{if \code{TRUE}, the observations with \code{NA} values
    will be removed.
  }
  \item{method}{the method to use, either \dQuote{delong},
    \dQuote{bootstrap} or \dQuote{venkatraman}. The first letter is sufficient. If omitted, the
	appropriate method is selected as explained in details.
  }
  \item{sensitivity, specificity}{if \code{method="sensitivity"} or
    \code{method="specificity"}, the respective level where the test
    must be assessed as a numeric of length 1. }
  \item{alternative}{specifies the alternative hypothesis. Either of
    \dQuote{two.sided}, \dQuote{less} or \dQuote{greater}. The first letter is
    sufficient. Default: \dQuote{two.sided}. Only \dQuote{two.sided} is available
    with \code{method="venkatraman"}.
  }
  \item{paired}{a logical indicating whether you want a paired roc.test.
    If \code{NULL}, the paired status will be auto-detected by \code{\link{are.paired}}.
    If \code{TRUE} but the paired status cannot be assessed by \code{\link{are.paired}}
    will produce an error.
  }
  \item{reuse.auc}{if \code{TRUE} (default) and the \dQuote{roc} objects
    contain an \dQuote{auc} field, re-use these specifications for the
    test. See the \emph{AUC specification} section for more details.

  }
  \item{boot.n}{for \code{method="bootstrap"} and \code{method="venkatraman"} only: the number of
    bootstrap replicates or permutations. Default: \var{2000}.
  }
  \item{boot.stratified}{for \code{method="bootstrap"} only:
    should the bootstrap be stratified (same number 
    of cases/controls in each replicate than in the original sample) or
    not. Ignored with \code{method="venkatraman"}. Default: \var{TRUE}.
  }
  \item{ties.method}{for \code{method="venkatraman"} only: argument for
    \code{\link{rank}} specifying how ties are handled. Defaults to
    \dQuote{first} as described in the paper.
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
    Arguments for \code{\link{auc}},
    and \code{\link{txtProgressBar}} (only \code{char} and \code{style})
    if applicable.
  }
}

\details{
  This function compares two ROC curves. It is typically called with the two \link{roc} objects to
  compare. \code{roc.test.default} is provided as a convenience
  method and creates two \link{roc} objects before calling
  \code{roc.test.roc}.

  Three methods are available: \dQuote{delong}, \dQuote{bootstrap} and \dQuote{venkatraman} (see
  \dQuote{Computational details} section below). \dQuote{delong} and
  \dQuote{bootstrap} are tests over the AUC whereas \dQuote{venkatraman}
  compares the the ROC curves themselves.

  Default is to use \dQuote{delong} method except for comparison of partial AUC, smoothed
  curves and curves with different \code{direction}, where \code{bootstrap}
  is used. Using \dQuote{delong} for partial AUC and smoothed ROCs is not
  supported in pROC and result in an error.
  It is spurious to use \dQuote{delong} for \code{\link{roc}} with different
  \code{direction} (a warning is issued but the spurious comparison is
  enforced). \dQuote{venkatraman}'s test cannot be employed to compare smoothed
  ROC curves, or curves with partial AUC specifications. In addition,
  and comparison of ROC curves with different
  \code{direction} should be used with care (a warning is produced as well).

  If \code{alternative="two.sided"}, a two-sided test for difference in AUC is performed. If
  \code{alternative="less"}, the alternative is that the AUC of roc1 is
  smaller than the AUC of roc2. For \code{method="venkatraman"}, only
  \dQuote{two.sided} test is available.
  
  If the \code{paired} argument is not provided, the \code{\link{are.paired}} function is
  employed to detect the paired status of the ROC curves. It will test if the original \code{response} is
  identical between the two ROC curves (this is always the case if the call is made with 
  \code{roc.test.default}). This detection is unlikely to raise false positives, but
  this possibility cannot be excluded entierly. It would require equal sample sizes
  and \code{response} values and order in both ROC curves. If it happens to you, use \code{paired=FALSE}.
  If you know the ROC curves are paired you can pass \code{paired=TRUE}. However this is
  useless as it will be tested anyway.

  For \link[=smooth.roc]{smoothed ROC curves}, smoothing is performed again at each
  bootstrap replicate with the parameters originally provided.
  If a density smoothing was performed with user-provided
  \code{density.cases} or \code{density.controls} the bootstrap cannot
  be performed and an error is issued.
}

\section{AUC specification}{
  The comparison of the AUC of the ROC curves needs a specification of the
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

  The AUC specification is ignored in the Venkatraman test.

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
    are computed and the difference is stored.
    \item The following formula is used:
    \deqn{D=\frac{AUC1-AUC2}{s}}{D=(AUC1-AUC2)/s}
    where s is the standard deviation of
    the bootstrap differences and AUC1 and AUC2 the AUC of the two
    (original) ROC curves.
    \item \var{D} is then compared to the normal distribution,
    according to the value of \code{alternative}.
  }

  See also the Bootstrap section in
  \link[=pROC-package]{this package's documentation}.

  With \code{method="delong"}, the processing is done as described in
  DeLong \emph{et al.} (1988) for paired ROC curves, using the algorithm
  of Sun and Xu (2014). Only comparison of
  two ROC curves is implemented. The method has been extended for
  unpaired ROC curves where the p-value is computed with an unpaired
  t-test with unequal sample size and unequal variance, with 
  \deqn{
    D=\frac{V^r(\theta^r) - V^s(\theta^s) }{ \sqrt{S^r + S^s}}
  }{
    D=(V^r(\theta^r) - V^s(\theta^s)) / sqrt(S^r + S^s)
  }

  With \code{method="venkatraman"}, the processing is done as described
  in Venkatraman and Begg (1996) (for paired ROC curves) and Venkatraman
  (2000) (for unpaired ROC curves) with \code{boot.n} permutation of
  sample ranks (with ties breaking). For consistency reasons, the same argument \code{boot.n} as
  in bootstrap defines the number of permutations to execute,
  even though no bootstrap is performed.

  For \code{method="specificity"}, the test assesses if the sensitivity of
  the ROC curves are different at the level of specificity given by the
  \code{specificity} argument, which must be a numeric of length 1. Bootstrap is employed as with \code{method="bootstrap"}
  and \code{boot.n} and \code{boot.stratified} are available. This is
  identical to the test proposed by Pepe \emph{et al.} (2009).
  The \code{method="sensitivity"} is very similar, but assesses if the specificity of
  the ROC curves are different at the level of sensitivity given by the
  \code{sensitivity} argument.
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

  If \code{method="venkatraman"}, and \code{alternative} is
  \dQuote{less} or \dQuote{greater}, the warning \dQuote{Only two-sided
    tests are available for Venkatraman. Performing two-sided test instead.}
  is produced and a two tailed test is performed.

  Both DeLong and Venkatraman's test ignores the direction of the ROC curve so that if two
  ROC curves have a different differ in the value of
  \code{direction}, the warning \dQuote{(DeLong|Venkatraman)'s test should not be
    applied to ROC curves with different directions.} is
  printed. However, the spurious test is enforced.

  If \code{boot.stratified=FALSE} and the sample has a large imbalance between
  cases and controls, it could happen that one or more of the replicates
  contains no case or control observation, or that there are not enough
  points for smoothing, producing a \code{NA} area.
  The warning \dQuote{NA value(s) produced during bootstrap were ignored.}
  will be issued and the observation will be ignored. If you have a large
  imbalance in your sample, it could be safer to keep
  \code{boot.stratified=TRUE}.
  
  When both ROC curves have an \code{\link{auc}} of 1 (or 100\%), their variances and covariance will always be null,
  and therefore the p-value will always be 1. This is true for both \dQuote{delong}, \dQuote{bootstrap} and 
  \dQuote{venkatraman} methods. This result is misleading, as the variances and covariance are of course not null.
  A \code{\link{warning}} will be displayed to inform of this condition, and of the misleading output.
}

\section{Errors}{
  An error will also occur if you give a \code{predictor2} when
  \code{predictor1} is a \code{\link{matrix}} or a
  \code{\link{data.frame}}, if \code{predictor1} has more than two
  columns, or if you do not give a \code{predictor2} when
  \code{predictor1} is a vector.
  
  If \code{density.cases} and \code{density.controls} were provided
  for smoothing, the error \dQuote{Cannot compute the statistic on ROC
    curves smoothed with density.controls and density.cases.} is
  issued.

  If \code{method="venkatraman"} and one of the ROC curves is smoothed,
  the error \dQuote{Using Venkatraman's test for smoothed ROCs is not
    supported.} is produced.

  With \code{method="specificity"}, the error \dQuote{Argument
    'specificity' must be numeric of length 1 for a specificity test.}
  is given unless the specificity argument is specified as a numeric of
  length 1. The \dQuote{Argument 'sensitivity' must be numeric of length
    1 for a sensitivity test.} message is given for
  \code{method="sensitivity"} under similar conditions.
}

\value{
  A list of class "htest" with following content:
  \item{p.value}{the p-value of the test.}
  \item{statistic}{the value of the Z (\code{method="delong"}) or D
    (\code{method="bootstrap"}) statistics.
  }
  \item{alternative}{the alternative hypothesis.}
  \item{method}{the character string \dQuote{DeLong's test for two
      correlated ROC curves} (if \code{method="delong"}) or
    \dQuote{Bootstrap test for two correlated ROC curves} (if
    \code{method="bootstrap"}).
  }
  \item{null.value}{the expected value of the statistic under the null
    hypothesis, that is 0.}
  \item{estimate}{the AUC in the two ROC curves.}
  \item{data.name}{the names of the data that was used.}
  \item{parameter}{for \code{method="bootstrap"} only: the values of the
    \code{boot.n} and \code{boot.stratified} arguments.
  }
}

\section{Acknowledgements}{
  We would like to thank E. S. Venkatraman and Colin B. Begg for their
  support in the implementation of their test.
}

\references{
  Elisabeth R. DeLong, David M. DeLong and Daniel L. Clarke-Pearson
  (1988) ``Comparing the areas under two or more correlated receiver
  operating characteristic curves: a nonparametric
  approach''. \emph{Biometrics} \bold{44}, 837--845.
  
  James A. Hanley and Barbara J. McNeil (1982) ``The meaning and use of
  the area under a receiver operating characteristic (ROC)
  curve''. \emph{Radiology} \bold{143}, 29--36.

  Margaret Pepe, Gary Longton and Holly Janes (2009) ``Estimation and
  Comparison of Receiver Operating Characteristic Curves''. \emph{The
  Stata journal} \bold{9}, 1.
  
  Xavier Robin, Natacha Turck, Jean-Charles Sanchez and Markus Müller
  (2009) ``Combination of protein biomarkers''. \emph{useR! 2009}, Rennes.
  \url{https://www.r-project.org/nosvn/conferences/useR-2009/abstracts/user_author.html}
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
  
  Xu Sun and Weichao Xu (2014) ``Fast Implementation of DeLongs Algorithm for Comparing
  the Areas Under Correlated Receiver Operating Characteristic Curves''. \emph{IEEE Signal
  Processing Letters}, \bold{21}, 1389--1393. 
  DOI: \doi{10.1109/LSP.2014.2337313}.

  E. S. Venkatraman and Colin B. Begg (1996) ``A distribution-free
  procedure for comparing receiver operating characteristic curves from
  a paired experiment''. \emph{Biometrika} \bold{83}, 835--848.
  DOI: \doi{10.1093/biomet/83.4.835}.

  E. S. Venkatraman (2000) ``A Permutation Test to Compare Receiver
  Operating Characteristic Curves''. \emph{Biometrics} \bold{56},
  1134--1138. DOI: \doi{10.1111/j.0006-341X.2000.01134.x}.
  
  Hadley Wickham (2011) ``The Split-Apply-Combine Strategy for Data Analysis''. \emph{Journal of Statistical Software}, \bold{40}, 1--29.
  URL: \href{https://www.jstatsoft.org/v40/i01}{www.jstatsoft.org/v40/i01}.
}

\seealso{
  \code{\link{roc}}, \code{\link{power.roc.test}}
  
  CRAN package \pkg{plyr}, employed in this function.
}

\examples{
data(aSAH)

# Basic example with 2 roc objects
roc1 <- roc(aSAH$outcome, aSAH$s100b)
roc2 <- roc(aSAH$outcome, aSAH$wfns)
roc.test(roc1, roc2)

\dontrun{
# The latter used Delong's test. To use bootstrap test:
roc.test(roc1, roc2, method="bootstrap")
# Increase boot.n for a more precise p-value:
roc.test(roc1, roc2, method="bootstrap", boot.n=10000)
}

# Alternative syntaxes
roc.test(aSAH$outcome, aSAH$s100b, aSAH$wfns)
roc.test(aSAH$outcome, data.frame(aSAH$s100b, aSAH$wfns))

# If we had a good a priori reason to think that wfns gives a
# better classification than s100b (in other words, AUC of roc1
# should be lower than AUC of roc2):
roc.test(roc1, roc2, alternative="less")

\dontrun{
# Comparison can be done on smoothed ROCs
# Smoothing is re-done at each iteration, and execution is slow
roc.test(smooth(roc1), smooth(roc2))
# or:
roc.test(aSAH$outcome, aSAH$s100b, aSAH$wfns, smooth=TRUE, boot.n=100)
}
# or from an AUC (no smoothing)
roc.test(auc(roc1), roc2)

\dontrun{
# Comparison of partial AUC:
roc3 <- roc(aSAH$outcome, aSAH$s100b, partial.auc=c(1, 0.8), partial.auc.focus="se")
roc4 <- roc(aSAH$outcome, aSAH$wfns, partial.auc=c(1, 0.8), partial.auc.focus="se")
roc.test(roc3, roc4)
# This is strictly equivalent to:
roc.test(roc3, roc4, method="bootstrap")

# Alternatively, we could re-use roc1 and roc2 to get the same result:
roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(1, 0.8), partial.auc.focus="se")

# Comparison on specificity and sensitivity
roc.test(roc1, roc2, method="specificity", specificity=0.9)
roc.test(roc1, roc2, method="sensitivity", sensitivity=0.9)
}

# Spurious use of DeLong's test with different direction:
roc5 <- roc(aSAH$outcome, aSAH$s100b, direction="<")
roc6 <- roc(aSAH$outcome, aSAH$s100b, direction=">")
roc.test(roc5, roc6, method="delong")

\dontrun{
# Comparisons of the ROC curves
roc.test(roc1, roc2, method="venkatraman")
}

# Unpaired tests
roc7 <- roc(aSAH$outcome, aSAH$s100b)
# artificially create an roc8 unpaired with roc7
roc8 <- roc(aSAH$outcome[1:100], aSAH$s100b[1:100])
\dontrun{
roc.test(roc7, roc8, paired=FALSE, method="delong")
roc.test(roc7, roc8, paired=FALSE, method="bootstrap")
roc.test(roc7, roc8, paired=FALSE, method="venkatraman")
roc.test(roc7, roc8, paired=FALSE, method="specificity", specificity=0.9)
}
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{utilities}
\keyword{htest}
\keyword{roc}
