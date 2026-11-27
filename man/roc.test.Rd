\encoding{UTF-8}
\name{roc.test}
\alias{roc.test}
\alias{roc.test.default}
\alias{roc.test.roc}
\alias{roc.test.formula}
\alias{roc.test.auc}
\alias{roc.test.smooth.roc}
\alias{roc.test.list}

\title{
  Test one or more ROC curves
}
\description{
  This function tests one ROC curve against the absence of
  discrimination, or compares two ROC curves or a list of two or more,
  correlated (paired) or uncorrelated (unpaired). The DeLong, bootstrap
  and ML methods test the (partial) AUC: that of one curve against its
  value without discrimination, the difference between two curves, or
  contrasts of the AUCs of a list of curves. The Venkatraman method tests
  if two curves are perfectly superposed. The sensitivity and specificity
  methods test if the sensitivity (respectively specificity) of two ROC
  curves are different at the given level of specificity (respectively
  sensitivity).
}
\usage{
# roc.test(...)
\S3method{roc.test}{roc}(roc1, roc2, method=c("delong", "bootstrap",
"venkatraman", "sensitivity", "specificity", "ml"), sensitivity = NULL,
specificity = NULL, alternative = c("two.sided", "less", "greater"),
paired=NULL, reuse.auc=TRUE, boot.n=2000, boot.stratified=TRUE,
ties.method="random", progress = getOption("pROCProgress", interactive()),
parallel=FALSE, conf.level=0.95, cl = NULL, ...)
\S3method{roc.test}{auc}(roc1, roc2, ...)
\S3method{roc.test}{smooth.roc}(roc1, roc2, ...)
\S3method{roc.test}{list}(rocs, contrast=NULL, method=c("delong", "bootstrap", "ml"),
alternative = c("two.sided", "less", "greater"), paired=NULL, boot.n=2000,
boot.stratified=TRUE, progress = getOption("pROCProgress", interactive()),
conf.level=0.95, cl = NULL, ...)
\S3method{roc.test}{formula}(formula, data, ..., subset, na.action)
\S3method{roc.test}{default}(response, predictor1, predictor2=NULL,
na.rm=TRUE, method=NULL, paired=TRUE, ...)
}

\arguments{
  \item{roc1, roc2}{the ROC curve to test, or the two ROC curves to
	compare. Either \dQuote{\link{roc}}, \dQuote{\link{auc}} or
	\dQuote{\link{smooth.roc}} objects (types can be mixed).
	If \code{roc2} is omitted, \code{roc1} is tested alone (see the
	\dQuote{Test of one ROC curve} section).
  }
  \item{rocs}{a list of two or more ROC curves, as
    \dQuote{\link{roc}}, \dQuote{\link{auc}} or
    \dQuote{\link{smooth.roc}} objects. Paired curves must be built on
    the same observations.
  }
  \item{contrast}{the contrasts of the AUCs to test against 0: a numeric
    vector with one coefficient per curve in \code{rocs}, or a matrix
    with one such row per contrast. If \code{NULL} (default), the
    equality of all the AUCs is tested.
  }
  \item{response}{a vector or factor, as for the \link{roc} function.
    The default and formula methods, which build the curves, are
    superseded: build them with \code{\link{roc}} first (see
    \link{all_in_one_calls}).}
  \item{predictor1}{a numeric or ordered vector as for the \link{roc} function, or
  a matrix or data.frame with predictors two columns.}
  \item{predictor2}{only if predictor1 was a vector, the second
    predictor as a numeric vector. If omitted, the curve of
    \code{predictor1} is tested alone.
  }
  \item{formula}{a formula of the type response~predictor1+predictor2,
    or response~predictor to test one curve.
  	Additional arguments  \code{data}, \code{subset} and \code{na.action}
	 are supported, see \code{\link{model.frame}} for more details.
  }
  \item{data}{a matrix or data.frame containing the variables in the
    formula. See \code{\link{model.frame}} for more details.}
  \item{subset, na.action}{arguments for \code{\link{model.frame}},
    used with \code{formula} and \code{data} to build the \code{\link{roc}} curves.
  }
  \item{na.rm}{if \code{TRUE}, the observations with \code{NA} values
    in the response or in either predictor will be removed before the
    two curves are built (only when \code{predictor1} and
    \code{predictor2} are vectors). Setting it to \code{FALSE} does not
    keep the \code{NA}s: \code{\link{roc}} still removes them from each
    curve, and a paired test then uses the observations complete in both.
  }
  \item{method}{the method to use, either \dQuote{delong},
    \dQuote{bootstrap}, \dQuote{venkatraman}, \dQuote{sensitivity},
    \dQuote{specificity} or \dQuote{ml}. It can be abbreviated, as long as it is not
    ambiguous (\dQuote{s} is not sufficient). If omitted, the
	appropriate method is selected as explained in details.
  }
  \item{sensitivity, specificity}{if \code{method="sensitivity"} or
    \code{method="specificity"}, the respective level where the test
    must be assessed as a numeric of length 1. }
  \item{alternative}{specifies the alternative hypothesis. Either of
    \dQuote{two.sided}, \dQuote{less} or \dQuote{greater}. The first letter is
    sufficient. Default: \dQuote{two.sided}. Only \dQuote{two.sided} is available
    with \code{method="venkatraman"} and with several contrasts.
  }
  \item{paired}{a logical indicating whether you want a paired roc.test.
    If \code{NULL}, the paired status will be auto-detected by \code{\link{are.paired}}.
    \code{TRUE} requires curves that can be paired (see Details), and
    \code{FALSE} treats the curves as unpaired. The curves of a list are
    paired if they all are; if only some are, they are all tested as
    unpaired. Does not apply to the test of one curve.
  }
  \item{reuse.auc}{if \code{TRUE} (default) and the \dQuote{roc} objects
    contain an \dQuote{auc} field, re-use these specifications for the
    test. See the \emph{AUC specification} section for more details.
    Superseded: pass an \code{auc} object instead (see
    \link{all_in_one_calls}).

  }
  \item{boot.n}{for \code{method="bootstrap"}, \code{"sensitivity"},
    \code{"specificity"} and \code{"venkatraman"} only: the number of
    bootstrap replicates or permutations. Default: \var{2000}.
  }
  \item{boot.stratified}{for \code{method="bootstrap"},
    \code{"sensitivity"} and \code{"specificity"} only:
    should the bootstrap be stratified (same number 
    of cases/controls in each replicate than in the original sample) or
    not. Ignored with \code{method="venkatraman"}. Default: \var{TRUE}.
  }
  \item{ties.method}{for \code{method="venkatraman"} only: argument for
    \code{\link{rank}} specifying how ties are handled. Defaults to
    \dQuote{random}. With paired curves, the ties of both
    predictors are broken with the same random order of the observations,
    so that identical predictors give identical ranks. No random numbers
    are drawn when there are no ties. The test is designed for continuous
    predictors: with many ties, the p-value depends on the random
    tie-breaking (Venkatraman and Begg, 1996, section 5).
  }
  \item{progress}{whether to display a text progress bar during the
    bootstrap. Defaults to \code{getOption("pROCProgress", interactive())},
    that is on in an interactive session and off otherwise. The names of the
    progress bars of the pre-1.19 \pkg{plyr} interface (\code{"none"},
    \code{"text"}, \code{"win"}, \code{"tk"}) are still accepted, but only a
    text bar is drawn.}
  \item{cl}{how to run the bootstrap in parallel. \code{NULL} (the default)
    or \code{FALSE} runs it sequentially; a cluster object from
    \code{\link[parallel]{makeCluster}} runs it on that cluster, which
    \pkg{pROC} neither creates nor stops; \code{TRUE} uses the cluster
    registered with \code{\link[parallel]{setDefaultCluster}}; and a number
    creates a cluster of that many workers for the duration of the call.
    Any cluster type works, including \code{"FORK"}, \code{"PSOCK"},
    \code{"MPI"} and \code{"MIRAI"}. Results are reproducible from
    \code{\link{set.seed}} whatever the number of workers, but differ from a
    sequential run with the same seed. See the \dQuote{Parallel processing}
    section of \link{pROC-package}.}
  \item{parallel}{DEPRECATED. A value other than \code{FALSE} will produce a warning.
    The argument will be removed in a future version.}
  \item{conf.level}{the width of the confidence interval as [0,1], never
    in percent, used for any calculated CI. Default: 0.95, resulting in a
    95\% CI.}
  \item{\dots}{further arguments passed to or from other methods,
    especially arguments for \code{\link{roc}} and \code{roc.test.roc}
    when calling \code{roc.test.default} or \code{roc.test.formula}.
    Arguments for \code{\link{auc}} if applicable. Not used with a list
    of curves.
  }
}

\details{
  The curves are given as \dQuote{\link{roc}}, \dQuote{\link{auc}} or
  \dQuote{\link{smooth.roc}} objects: two of them, a list of two or
  more, or a single one (see the \dQuote{Test of one ROC curve}
  section).

  The \dQuote{delong}, \dQuote{bootstrap} and \dQuote{ml} methods test
  the AUC, \dQuote{venkatraman} the ROC curves themselves, and
  \dQuote{sensitivity} and \dQuote{specificity} one point of the curves
  (see the \dQuote{Computation details} section below). Only the first
  three are available for one curve and for a list of curves.

  If \code{method} is not provided, it is the first that applies of:
  \enumerate{
    \item \dQuote{bootstrap} for a partial AUC;
    \item \dQuote{ml} for unpaired curves all smoothed with
    \code{method="binormal"} or \code{"proper"}, or a single such
    curve;
    \item \dQuote{bootstrap} for other smoothed curves, and for curves
    with different \code{direction};
    \item \dQuote{delong}.
  }
  \dQuote{delong} and \dQuote{venkatraman} do not support partial AUCs
  nor smoothed curves, and \dQuote{ml} only supports the full AUC of
  unpaired binormal or proper smoothed curves. \dQuote{delong} and
  \dQuote{venkatraman} ignore the \code{direction} of the curves: with
  different directions the comparison is spurious (a warning is issued
  but the comparison is enforced).

  The \code{alternative} is on the AUC of a single curve, compared to
  its value without discrimination, on the difference between the AUC of
  roc1 and that of roc2, or on the \code{contrast} of a list of curves:
  with \code{alternative="less"}, that the AUC of the single curve is
  smaller than without discrimination, that the AUC of roc1 is smaller
  than the AUC of roc2, or that the \code{contrast} is smaller than 0.

  If the \code{paired} argument is not provided, the
  \code{\link{are.paired}} function is employed to detect the paired
  status of the ROC curves, and a message tells the result. The detection
  only shows that the curves can be paired: use \code{paired=FALSE} for
  independent samples that have the same \code{response}, for which a
  paired test is not valid. Curves built with \code{controls} and
  \code{cases} are not paired automatically: with the same numbers of
  controls and cases as given (missing values included), a warning is
  given and \code{paired=TRUE} pairs them by position. The curves built
  by \code{roc.test.default} and \code{roc.test.formula} share their
  response and are paired.

  Two paired curves with different missing values are re-paired as in
  \code{\link{are.paired}}; the test and the \code{estimate} refer to
  the re-paired curves, returned as \code{roc1} and \code{roc2}. The
  curves of a list are not re-paired.

  For \link[=smooth.roc]{smoothed ROC curves}, smoothing is performed again at each
  bootstrap replicate with the parameters originally provided.
  If a density smoothing was performed with user-provided
  \code{density.cases} or \code{density.controls} the bootstrap cannot
  be performed and an error is issued.
}

\section{Test of one ROC curve}{
  A single curve is tested against the absence of discrimination: with
  \eqn{\theta} its AUC and \var{S} the variance of \eqn{\theta} (see
  \code{\link[=var.roc]{var}}),
  \deqn{
    Z=\frac{\theta - \theta_0}{\sqrt{S}}
  }{
    Z=(\theta - \theta0) / sqrt(S)
  }
  is compared to the normal distribution (DeLong \emph{et al.}, 1988,
  p. 841, with \dQuote{delong}; Dorfman and Alf, 1969, p. 493--494, with
  \dQuote{ml}). It is named \var{D} with \code{method="bootstrap"}.
  \eqn{\theta_0}{\theta0} is 0.5 (50 with \code{percent=TRUE}), or for a
  partial AUC the area under the diagonal in the same region (0.5 if
  \code{partial.auc.correct=TRUE}). The contrasts of a list of curves
  are tested against 0 instead: a contrast on one curve of a list is not
  this test.

  \code{conf.int} is the confidence interval of \eqn{\theta},
  \eqn{\theta \pm z \sqrt{S}}{\theta +/- z sqrt(S)}, so that it agrees
  with the p-value. With \code{method="delong"} or \code{"ml"} and a
  two-sided alternative, it is the interval of \code{\link{ci.auc}} with
  \code{scale="auc"}, without truncation; with the bootstrap,
  \code{\link{ci.auc}} takes the percentiles instead.

  With \code{direction="auto"}, \code{\link{roc}} selects the direction
  from the data, and a one-sided test is not valid: set the
  \code{direction} explicitly for a one-sided test.

  The AUC of a curve smoothed with \code{method="proper"} cannot be
  lower than 0.5: only \code{alternative="greater"} is meaningful, and
  the normal approximation is doubtful for an AUC close to 0.5.

  An AUC of 0 or 1 has a null variance: the statistic is infinite.
}

\section{AUC specification}{
  For a partial AUC, pass an \code{auc} object from \code{\link{auc}}
  instead of a curve. With a curve, the \dQuote{auc} field of the curve is
  used, or the partial AUC arguments in \code{\dots}, depending on
  \code{reuse.auc}: this is superseded (see \link{all_in_one_calls}).

  The Venkatraman test compares the whole curves: it fails with an error if
  either curve has a partial AUC specification.
}


\section{Computation details}{
  This section describes the comparison of two curves and of a list of
  curves; see the \dQuote{Test of one ROC curve} section for a single
  curve.

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
    \item The confidence interval is \eqn{AUC1-AUC2 \pm z s}{AUC1-AUC2 +/- z s},
    with \var{z} the quantile of the normal distribution, so that it
    agrees with the p-value.
  }

  See also the Bootstrap section in
  \link[=pROC-package]{this package's documentation}.

  With \code{method="delong"}, the processing is done as described in
  DeLong \emph{et al.} (1988) for paired ROC curves, using the algorithm
  of Sun and Xu (2014). For unpaired ROC curves, the statistic
  \deqn{
    D=\frac{\theta^r - \theta^s}{\sqrt{S^r + S^s}}
  }{
    D=(\theta^r - \theta^s) / sqrt(S^r + S^s)
  }
  is compared to the normal distribution (Hanley and McNeil, 1983).

  With \code{method="ml"}, the variances \eqn{V^r}{V^r} and
  \eqn{V^s}{V^s} of the AUC of the two curves are computed as described in
  \code{\link{var}}, and the Wald statistic
  \deqn{
    Z=\frac{\theta^r - \theta^s}{\sqrt{V^r + V^s}}
  }{
    Z=(\theta^r - \theta^s) / sqrt(V^r + V^s)
  }
  is compared to the normal distribution (Dorfman and Alf, 1969, p. 494;
  Pesce and Metz, 2007, p. 822).

  For a list of \var{k} ROC curves, with \eqn{\theta} the vector
  of their AUCs, \var{S} its covariance matrix (see
  \code{\link[=cov.roc]{cov}}) and \var{L} the \code{contrast},
  DeLong \emph{et al.} (1988, p. 841) is applied. With
  \code{method="bootstrap"}, \var{S} is the covariance matrix of the
  (partial) AUCs of the \code{boot.n} bootstrap replicates, so that the
  test of two curves is the bootstrap test above. Unpaired curves are
  independent: with \code{method="delong"} and \code{"ml"} \var{S}
  only holds their variances, and with \code{method="bootstrap"} each
  curve is resampled separately. With one contrast,
  \deqn{
    Z=\frac{L\theta}{\sqrt{LSL'}}
  }{
    Z=L\theta / sqrt(LSL')
  }
  is compared to the normal distribution, and \code{conf.int} is the
  confidence interval of \eqn{L\theta}. With several contrasts,
  \deqn{
    X^2=(L\theta)'(LSL')^{-1}(L\theta)
  }{
    X^2=(L\theta)' (LSL')^-1 (L\theta)
  }
  is compared to the chi-square distribution with the rank of
  \eqn{LSL'} as degrees of freedom (Equation 5). The default
  \code{contrast} is the \eqn{k-1} differences between successive
  curves. A contrast without variance, for instance between curves that
  rank the observations identically, gives a statistic of 0 if it is 0
  and an infinite statistic otherwise, as for two curves.

  With \code{method="venkatraman"}, the processing is done as described
  in Venkatraman and Begg (1996) (for paired ROC curves) and Venkatraman
  (2000) (for unpaired ROC curves) with \code{boot.n} permutation of
  sample ranks (ties broken at random, see \code{ties.method}). For consistency reasons, the same argument \code{boot.n} as
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
  Curves with different \dQuote{auc} specifications or different
  smoothing parameters are compared as they are, with a warning.
  Smoothing arguments are compared as they were given (see
  \code{smoothing.args} in the Value section of
  \code{\link[=smooth.roc]{smooth}}): an explicit default value or an
  abbreviation (such as \code{kernel="e"} for \code{"epanechnikov"})
  counts as a difference. To smooth curves the same way, call
  \code{\link[=smooth.roc]{smooth}} on each of them with the same
  arguments rather than using \code{roc(\dots, smooth=TRUE)}, where the
  arguments meant for other functions (such as plot) also reach the
  smoothing.

  Bootstrap replicates without an area are ignored, with a warning: those
  without a case or control observation, which can happen with
  \code{boot.stratified=FALSE} and a large imbalance between cases and
  controls, and those without enough points for smoothing.

  When all the curves have an \code{\link{auc}} of 1 (or 100\%), their
  variances and covariances are estimated as null and the p-value of
  a difference is always 1. This result is misleading, and comes with a warning.
}

\value{
  A list of class "htest" with following content:
  \item{p.value}{the p-value of the test.}
  \item{statistic}{the value of the Z (\code{method="delong"}, paired
    curves or one curve, and \code{method="ml"}), D (\code{method="delong"} with unpaired curves,
    \code{method="bootstrap"}, \code{"sensitivity"} and
    \code{"specificity"}) or E (\code{method="venkatraman"}) statistics.
    For a list of curves, Z or D as above with one contrast, and
    X-squared with several.
  }
  \item{conf.int}{the confidence interval of the difference or contrast tested,
    or of the AUC of one curve
    (not returned with \code{method="venkatraman"} or with several contrasts). Has an attribute \code{conf.level} specifying the level of the test.}
  \item{alternative}{the alternative hypothesis.}
  \item{method}{a character string describing the test, such as
    \dQuote{DeLong's test for one ROC curve},
    \dQuote{DeLong's test for two correlated ROC curves} (if
    \code{method="delong"} with paired curves) or
    \dQuote{Bootstrap test for two ROC curves} (if
    \code{method="bootstrap"} with unpaired curves).
  }
  \item{null.value}{the expected value of the statistic under the null
    hypothesis: 0 for a difference or a contrast, and the AUC without
    discrimination for one curve (see the \dQuote{Test of one ROC
    curve} section). Not returned with several contrasts.}
  \item{estimate}{the AUC of the ROC curve or curves (with two paired curves,
    the curves rebuilt on the observations complete in both, see
    Details), or the specificities (\code{method="sensitivity"}) or
    sensitivities (\code{method="specificity"}) of the two curves at
    the given value. Not returned with \code{method="venkatraman"}.}
  \item{data.name}{the names of the data that was used.}
  \item{roc1, roc2}{the ROC curve that was tested (\code{roc1} only),
    or the two ROC curves that were compared (with paired
    curves, rebuilt on the observations complete in both, see Details).}
  \item{parameter}{for \code{method="bootstrap"}, \code{"sensitivity"}
    and \code{"specificity"}: the values of the \code{boot.n} and
    \code{boot.stratified} arguments. For \code{method="venkatraman"}: the
    value of \code{boot.n}. For a list of curves with several
    contrasts: also the degrees of freedom.
  }
  \item{rocs, contrast}{for a list of curves: the ROC curves that were
    compared, and the matrix of contrasts that was tested.}
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

  Donald D. Dorfman, Edward Alf, Jr. (1969) ``Maximum-likelihood
  estimation of parameters of signal-detection theory and determination
  of confidence intervals---rating-method data''. \emph{Journal of
  Mathematical Psychology} \bold{6}, 487--496.
  DOI: \doi{10.1016/0022-2496(69)90019-4}.
  
  James A. Hanley and Barbara J. McNeil (1982) ``The meaning and use of
  the area under a receiver operating characteristic (ROC)
  curve''. \emph{Radiology} \bold{143}, 29--36.

  James A. Hanley and Barbara J. McNeil (1983) ``A method of comparing
  the areas under receiver operating characteristic curves derived from
  the same cases''. \emph{Radiology} \bold{148}, 839--843.

  Margaret Pepe, Gary Longton and Holly Janes (2009) ``Estimation and
  Comparison of Receiver Operating Characteristic Curves''. \emph{The
  Stata journal} \bold{9}(1), 1--16. DOI: \doi{10.1177/1536867X0900900101}.

  Lorenzo L. Pesce, Charles E. Metz (2007) ``Reliable and computationally
  efficient maximum-likelihood estimation of "proper" binormal ROC
  curves''. \emph{Academic Radiology} \bold{14}, 814--829.
  DOI: \doi{10.1016/j.acra.2007.03.012}.
  
  Xavier Robin, Natacha Turck, Jean-Charles Sanchez and Markus Müller
  (2009) ``Combination of protein biomarkers''. \emph{useR! 2009}, Rennes.
  \url{https://www.r-project.org/nosvn/conferences/useR-2009/abstracts/user_author.html}
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{12}, 77.
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
}

\seealso{
  \code{\link{roc}}, \code{\link{power.roc.test}}
}

\examples{
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)
roc2 <- roc(aSAH$outcome, aSAH$wfns)

# One curve: is the AUC different from 0.5?
roc.test(roc1)

# Two curves: DeLong's test, the default
roc.test(roc1, roc2)

# One-sided: AUC of roc1 is lower than AUC of roc2
roc.test(roc1, roc2, alternative = "less")

# From an AUC, without smoothing
roc.test(auc(roc1), roc2)

# Force an unpaired test
roc.test(roc1, roc2, paired = FALSE)

# More than two paired curves: are all the AUCs equal?
rocs <- list(s100b = roc1, wfns = roc2, ndka = roc(aSAH$outcome, aSAH$ndka))
roc.test(rocs)
# One contrast: s100b against the mean of the other two
roc.test(rocs, contrast = c(1, -0.5, -0.5))
# Several contrasts, one per row
roc.test(rocs, contrast = rbind(c(1, -1, 0), c(1, 0, -1)))
# Partial AUCs are bootstrapped
paucs <- lapply(rocs, auc, partial.auc = c(1, 0.8))
\dontrun{
roc.test(paucs)}\dontshow{
roc.test(paucs, boot.n = 10)}

# Bootstrap
\dontrun{
roc.test(roc1, roc2, method = "bootstrap")}\dontshow{
roc.test(roc1, roc2, method = "bootstrap", boot.n = 10)}

# Unpaired binormal or proper smoothed curves: ML method
male <- aSAH$gender == "Male"
roc.test(smooth(roc(aSAH$outcome[male], aSAH$s100b[male])),
         smooth(roc(aSAH$outcome[!male], aSAH$s100b[!male])))

# Paired smoothed curves, and other smoothings, are bootstrapped, and smoothing
# is redone at each replicate
\dontrun{
roc.test(smooth(roc1), smooth(roc2))}\dontshow{
roc.test(smooth(roc1), smooth(roc2), boot.n = 10)}

# So is the partial AUC, given as auc objects
pauc1 <- auc(roc1, partial.auc = c(1, 0.8), partial.auc.focus = "se")
pauc2 <- auc(roc2, partial.auc = c(1, 0.8), partial.auc.focus = "se")
\dontrun{
roc.test(pauc1, pauc2)}\dontshow{
roc.test(pauc1, pauc2, boot.n = 10)}

# Comparison of the whole ROC curves, or at a given specificity
\dontrun{
roc.test(roc1, roc2, method = "venkatraman")
roc.test(roc1, roc2, method = "specificity", specificity = 0.9)}\dontshow{
roc.test(roc1, roc2, method = "venkatraman", boot.n = 10)
roc.test(roc1, roc2, method = "specificity", specificity = 0.9, boot.n = 10)}
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{utilities}
\keyword{htest}
\keyword{roc}
