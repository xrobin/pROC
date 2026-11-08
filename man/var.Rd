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
\S3method{var}{roc}(roc, method=c("delong", "bootstrap", "obuchowski", "ml"),
boot.n = 2000, boot.stratified = TRUE, reuse.auc=TRUE, 
progress = getOption("pROCProgress", interactive()), parallel=FALSE, cl = NULL, ...)
\S3method{var}{smooth.roc}(smooth.roc, ...)
}

\arguments{
  \item{roc, smooth.roc, auc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, a \dQuote{smooth.roc} object from the
	\code{\link[=smooth.roc]{smooth}} function or an \dQuote{auc} object from
	the \code{\link{auc}} function.
  }
  \item{method}{the method to use, either \dQuote{delong},
    \dQuote{bootstrap}, \dQuote{obuchowski} or \dQuote{ml}. The first letter is
    sufficient. If omitted, the appropriate method is selected as
    explained in details.
  }
  \item{reuse.auc}{if \code{TRUE} (default) and the \dQuote{roc} objects
    contain an \dQuote{auc} field, re-use these specifications for the
    test. See details.
    Superseded: pass an \code{auc} object instead (see
    \link{all_in_one_calls}).
  }
  \item{boot.n}{for \code{method="bootstrap"} only: the number of
    bootstrap replicates or permutations. Default: \var{2000}.
  }
  \item{boot.stratified}{for \code{method="bootstrap"} only:
    should the bootstrap be stratified (same number 
    of cases/controls in each replicate than in the original sample) or
    not. Default: \var{TRUE}.
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
  \item{\dots}{further arguments passed to or from other methods,
    especially arguments for \code{var.roc} when calling \code{var},
    \code{var.auc} and \code{var.smooth.roc}. Arguments for
    \code{\link{auc}} (if \code{reuse.auc=FALSE}).
  }
}

\details{
  The \code{var} function computes the variance of the AUC of a ROC
  curve. It is typically called with the \link{roc} object of
  interest. Four methods are available: \dQuote{delong},
  \dQuote{bootstrap}, \dQuote{obuchowski} and \dQuote{ml} (see \dQuote{Computational
  details} section below).

  The default is to use \dQuote{delong} method except for with
  partial AUC and smoothed curves. The full AUC of curves smoothed with
  \code{method="binormal"} or \code{"proper"} uses \dQuote{ml}; partial
  AUC and other smoothed curves use \dQuote{bootstrap}.
  Using \dQuote{delong} for partial AUC and smoothed ROCs is not
  supported, nor \dQuote{ml} for anything else than the full AUC of
  binormal and proper smoothed curves.

  For \link[=smooth.roc]{smoothed ROC curves}, smoothing is performed again at each
  bootstrap replicate with the parameters originally provided.
  If a density smoothing was performed with user-provided
  \code{density.cases} or \code{density.controls} the bootstrap cannot
  be performed and an error is issued.

  \code{var.default} forces the usage of the
  \code{\link[stats:cor]{var}} function in the \pkg{stats} package, so
  that other code relying on \code{var} should continue to function
  normally.
}

\section{AUC specification}{
  For a partial AUC, pass an \code{auc} object from \code{\link{auc}}
  instead of a curve. With a curve, the \dQuote{auc} field of the curve is
  used, or the partial AUC arguments in \code{\dots}, depending on
  \code{reuse.auc}: this is superseded (see \link{all_in_one_calls}).
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
  Hanley and Hajian-Tilaki (1997) using the algorithm by Sun and Xu (2014).

  With \code{method="obuchowski"}, the processing is done as described
  in Obuchowski and McClish (1997), Table 1 and Equation 4, p. 1530--1531. The
  formulas of Table 1 are used as published. The binormal parameters
  \eqn{A} and \eqn{B} are their maximum-likelihood estimates (pp. 1530--1531),
  those of pROC's binormal smoothing of the curve (see
  \code{\link[=smooth.roc]{smooth}}).

  With \code{method="ml"}, the variance is that of the AUC of the
  maximum-likelihood fit of the smoothed curve, by the delta method, from
  the inverse of the expected (Fisher) information matrix of the fit at
  its maximum (Dorfman and Alf, 1969, p. 490):
  \itemize{
    \item for \code{method="binormal"}, with Equation 6 of Obuchowski and
    McClish (1997, p. 1531) and the derivatives of the AUC with respect to
    \eqn{a} and \eqn{b} in their Table I;
    \item for \code{method="proper"}, as described by Pesce and Metz
    (2007, pp. 820--821): when the smallest eigenvalue of the information
    matrix is less than \eqn{10^{-6}}{1e-6} times the next one, its
    eigenvector is left out of the inverse (Moore-Penrose pseudo-inverse).
    At the chance line (an AUC of 0.5), the variance is unstable (Pesce
    and Metz, 2007, p. 828), and with a single point it does not include
    the arbitrary AUC of the curve (see \code{\link[=smooth.roc]{smooth}}):
    a warning is issued in both cases.
  }
}

\section{Binormality assumption}{
  The \dQuote{obuchowski} method makes the assumption that the data is binormal.
  If the data shows a deviation from this assumption, it might help to
  normalize the data first (that is, before calling \code{\link{roc}}),
  for example with quantile normalization: 
  \preformatted{
    norm.x <- qnorm(rank(x)/(length(x)+1))
    var(roc(response, norm.x, ...), ...)
  }

  \dQuote{delong} and \dQuote{bootstrap} methods make no such assumption.
}

\value{
  The numeric value of the variance.
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
  
  When the ROC curve has an \code{\link{auc}} of 1 (or 100\%), the variance will always be null.
  This is true for the \dQuote{delong}, \dQuote{bootstrap} and \dQuote{ml} methods that can
  not properly assess the variance in this case. This result is misleading, as the variance is of course not null.
  A \code{\link{warning}} will be displayed to inform of this condition, and of the misleading output.
}

\section{Errors}{
  If the ROC curve was smoothed with numeric \code{density.controls}
  and \code{density.cases}, an error is issued.

  If \code{method="delong"} and the AUC specification specifies a
  partial AUC, or the ROC curve is smoothed, an error is issued. The same
  applies to \code{method="obuchowski"} with a smoothed curve or with
  \code{partial.auc.focus="sensitivity"}, and to \code{method="ml"}
  with a partial AUC or a curve that was not smoothed with
  \code{method="binormal"} or \code{"proper"}.
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

  James A. Hanley and Karim O. Hajian-Tilaki (1997) ``Sampling
  variability of nonparametric estimates of the areas under receiver
  operating characteristic curves: An update''. \emph{Academic
  Radiology} \bold{4}, 49--58. DOI:
  \doi{10.1016/S1076-6332(97)80161-4}.

  Nancy A. Obuchowski, Donna K. McClish (1997). ``Sample size
  determination for diagnostic accuracy studies involving binormal ROC
  curve indices''. \emph{Statistics in Medicine}, \bold{16}(13),
  1529--1542. DOI: \doi{10.1002/(SICI)1097-0258(19970715)16:13<1529::AID-SIM565>3.0.CO;2-H}.

  Lorenzo L. Pesce, Charles E. Metz (2007) ``Reliable and computationally
  efficient maximum-likelihood estimation of "proper" binormal ROC
  curves''. \emph{Academic Radiology} \bold{14}, 814--829.
  DOI: \doi{10.1016/j.acra.2007.03.012}.
  
  Xu Sun and Weichao Xu (2014) ``Fast Implementation of DeLongs Algorithm for Comparing
  the Areas Under Correlated Receiver Operating Characteristic Curves''. \emph{IEEE Signal
  Processing Letters}, \bold{21}, 1389--1393. 
  DOI: \doi{10.1109/LSP.2014.2337313}.
}

\seealso{
 \code{\link{roc}}, \code{\link{cov.roc}}
}

\examples{
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)

# DeLong, the default
var(roc1)

# Obuchowski
var(roc1, method = "obuchowski")

# From an AUC, without smoothing
var(auc(roc1))

# Bootstrap
\dontrun{
var(roc1, method = "bootstrap")}\dontshow{
var(roc1, method = "bootstrap", boot.n = 10)}

# Binormal and proper smoothed curves: ML method
var(smooth(roc1))

# Other smoothed curves are bootstrapped, and smoothing is redone at each replicate
\dontrun{
var(smooth(roc1, method = "density"))}\dontshow{
var(smooth(roc1, method = "density"), boot.n = 10)}
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
