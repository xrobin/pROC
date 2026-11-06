\encoding{UTF-8}
\name{smooth}
\alias{smooth}
\alias{smooth.roc}
\alias{smooth.smooth.roc}
\alias{smooth.default}
\title{
 Smooth a ROC curve
}
\description{
  This function smoothes a ROC curve of numeric predictor. By default, a
  binormal smoothing is performed, but density or custom smoothings are
  supported.
}
\usage{
smooth(...)
\S3method{smooth}{default}(...)
\S3method{smooth}{roc}(roc,
method=c("binormal", "density", "fitdistr", "logcondens",
"logcondens.smooth", "proper"), n=512, bw = "nrd0", density=NULL,
density.controls=density, density.cases=density,
start=NULL, start.controls=start, start.cases=start, 
reuse.auc=TRUE, reuse.ci=FALSE, ...)
\S3method{smooth}{smooth.roc}(smooth.roc, ...)
}

\arguments{
  \item{roc, smooth.roc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, or a \dQuote{smooth.roc} object from the
	\code{\link[=smooth.roc]{smooth}} function.
  }
  \item{method}{\dQuote{binormal}, \dQuote{density}, \dQuote{fitdistr},
    \dQuote{logcondens}, \dQuote{"logcondens.smooth"}, \dQuote{proper}.
  }
  \item{n}{
    the number of equally spaced points where the smoothed curve will be
    calculated.
  }
  \item{bw}{
    if \code{method="density"} and \code{density.controls} and
    \code{density.cases} are not provided, \code{bw} is passed to
    \code{\link{density}} to determine the bandwidth of the
    density Can be a character string (\dQuote{nrd0}, \dQuote{nrd},
    \dQuote{ucv}, \dQuote{bcv} or \dQuote{SJ}, but any name
    \link[base:match.fun]{matching} a function prefixed with \dQuote{bw.} is
    supported) or a numeric value, as described in \code{\link{density}}.
    Defaults to \dQuote{\link[stats:bandwidth]{nrd0}}.
  }
  \item{density, density.controls, density.cases}{if
    \code{method="density"}, a numeric value of density (over the y
    axis) or a function returning a density (such as
    \code{\link{density}}. Negative density values are set to 0, as
    \code{\link{density}} does. If \code{method="fitdistr"}, a \code{densfun}
    argument for \code{\link[MASS]{fitdistr}}.
    If the value is different for control and case observations,
    \code{density.controls} and \code{density.cases} can be employed
    instead, otherwise \code{density} will be propagated to both
	\code{density.controls} and \code{density.cases}.
  }
  \item{start, start.controls, start.cases}{if
    \code{method="fitdistr"}, \code{start}
    arguments for \code{\link[MASS]{fitdistr}}. They can be omitted for
    some of the distributions named by a character \code{density}, but
    are required when \code{density} is a function (and for the
    \dQuote{beta}, \dQuote{chi-squared} and \dQuote{f} distributions).
    \code{start.controls}
    and \code{start.cases} allow to specify different start values for
    controls and cases.
  }
  \item{reuse.auc, reuse.ci}{if \code{TRUE} (default for reuse.auc) and the \dQuote{roc} objects
    contain \dQuote{auc} or \dQuote{ci} fields, re-use these
    specifications to regenerate \code{\link{auc}} or \code{\link{ci}}
    on the smoothed ROC curve with the original parameters. If
    \code{FALSE}, the object returned will not contain
    \dQuote{auc} or \dQuote{ci} fields. It is currently not possible to
    redefine \code{\link{auc}} and \code{\link{ci}} options directly: you need to call \code{\link{auc}} or
    \code{\link{ci}} later for that. Superseded by calling \code{\link{auc}}
    and \code{\link{ci}} on the smoothed curve (see \link{all_in_one_calls}).
  }
  \item{\dots}{further arguments passed to or from other methods, and
    especially to \code{\link{density}} (only \code{cut}, \code{adjust},
    and \code{kernel}, plus \code{window} for compatibility with S+) and
    \code{\link[MASS]{fitdistr}}.
  }
}

\details{
  If \code{method="binormal"}, a linear model is fitted to the quantiles of
  the sensitivities and specificities. Smoothed sensitivities and
  specificities are then generated from this model on \code{n} points.
  This simple approach was found to work well for most ROC curves, but
  it may produce hooked smooths in some situations (see in Hanley (1988)). 

  With \code{method="proper"}, the \dQuote{proper} binormal model of Metz
  and Pan (1999) is fitted by maximum likelihood. Its decision
  variable is the likelihood ratio of a pair of normal distributions, so
  that the curve is always convex: it has no hook and its AUC is at least
  0.5. The curve is parametrized by \eqn{d_a}{d_a} and \eqn{c}{c}
  (Metz and Pan's Eqs. 33 and 34). The data are first grouped into
  categories: each run of controls or of cases in the sorted predictor
  values is one category, and tied control and case values form one
  category (Pesce and Metz, 2007). The fit starts from the least-squares
  line of \code{method="binormal"} and from curves with the empirical AUC,
  and maximizes the
  likelihood over \eqn{d_a}{d_a}, \eqn{c}{c} and the category boundaries
  together, as described by Pesce and Metz (2007), with
  \code{\link{optim}}'s BFGS method.
  Its full AUC is computed with Metz and Pan's Eq. 36. It
  can smooth curves that binormal smoothing cannot, such as degenerate
  data, but it is much slower: from a fraction of a second to tens of
  seconds per curve, growing with the number of categories, which matters
  with bootstrap confidence intervals.

  With \code{method="density"}, the \code{\link{density}}
  function is employed to generate a smooth kernel
  density of the control and case observations as described by Zou
  \emph{et al.} (1997), unless
  \code{density.controls} or \code{density.cases} are provided
  directly. \code{bw} can be given to
  specify a bandwidth to use with \code{\link{density}}. It can be a
  numeric value or a character string (\dQuote{nrd0}, \dQuote{nrd},
  \dQuote{ucv}, \dQuote{bcv} or \dQuote{SJ}, but any name
  \link[base:match.fun]{matching} a function prefixed with \dQuote{bw.} is
  supported). In the case of a character
  string, the whole predictor data is employed to determine the numeric
  value to use on both controls and cases.
  Depending on your data, it might be a good idea to specify the
  \code{kernel} argument for \code{\link{density}}. By default,
  \dQuote{gaussian} is used, but \dQuote{epanechnikov},
  \dQuote{rectangular}, \dQuote{triangular}, \dQuote{biweight},
  \dQuote{cosine} and \dQuote{optcosine} are supported. As all the
  kernels are symmetrical, it might help to normalize the data first
  (that is, before calling \code{\link{roc}}), for example with quantile
  normalization:
  \preformatted{
    norm.x <- qnorm(rank(x)/(length(x)+1))
    smooth(roc(response, norm.x, ...), ...)
  }

  Additionally, \code{density} can be a function which must return
  either a numeric vector of densities over the y axis or a \link{list}
  with a \dQuote{y} item like the \code{\link{density}} function. It
  must accept the following input:
  \preformatted{
    density.fun(x, n, from, to, bw, kernel, ...)
  }
  It is important to honour \code{n}, \code{from} and \code{to} in order
  to have the densities evaluated on the same points for controls and
  cases. Failing to do so and returning densities of different length
  will produce an error. It is also a good idea to use a constant
  smoothing parameter (such as \code{bw}) especially when controls and
  cases have a different number of observations, to avoid producing
  smoother or rougher densities.

  If \code{method="fitdistr"}, the \code{\link[MASS]{fitdistr}}
  function from the \pkg{MASS} package is employed to fit parameters for
  the density function \code{density} with start parameters
  \code{start} (required when \code{density} is a function). The density function are fitted
  separately in control (\code{density.controls}, \code{start.controls})
  and case observations (\code{density.cases},
  \code{start.cases}). \code{density} can be one of the character values
  allowed by \code{\link[MASS]{fitdistr}} or a density function (such
  as \code{\link{dnorm}}, \code{\link{dweibull}}, ...).
  For the distributions named by a character value, the smoothed curve is
  the exact ROC curve of the fitted distributions, computed from their
  cumulative distribution functions at \code{n} thresholds spread over
  the quantiles of both fits. With a density function, the fitted
  densities are evaluated at \code{n}
  equally spaced points over the range of the observations only, so that
  the curve is an approximation that ignores the tails of the fitted distributions outside of this range.
  The discrete distributions (\dQuote{geometric}, \dQuote{negative
  binomial} and \dQuote{poisson}) require non-negative integer (count)
  values.

  The \code{method="logcondens"} and \code{method="logcondens.smooth"} use the
    \pkg{logcondens} package to generate a non smoothed or smoothed
    (respectively) log-concave density estimate of of the control and case
    observation with the \link[logcondens]{logConROC} function.

  \code{smooth.default} forces the usage of the
  \code{\link[stats]{smooth}} function in the \pkg{stats} package, so
  that other code relying on \code{smooth} should continue to function
  normally.

  Smoothed ROC curves can be passed to smooth again. In this case, the
  smoothing is not re-applied on the smoothed ROC curve but the
  original \dQuote{\link{roc}} object will be re-used.
  
  Note that a \code{smooth.roc} curve has no threshold.
}

\value{
  A list of class \dQuote{smooth.roc} with the following fields:
  \item{sensitivities}{the smoothed sensitivities defining the ROC curve.}
  \item{specificities}{the smoothed specificities defining the ROC curve.}
  \item{percent}{if the sensitivities, specificities and AUC are
    reported in percent, as in the original \dQuote{\link{roc}} object.
  }
  \item{direction}{the direction of the comparison, as in the original
    \dQuote{\link{roc}} object.}
  \item{call}{how the function was called. See \code{\link{match.call}} for
    more details.
  }
  \item{smoothing.args}{a list of the arguments given to \code{smooth}:
    \code{method}, \code{n}, \code{bw}, the \code{density.*} and
    \code{start.*} arguments, and everything passed in \dots. With the
    \code{\link{roc}(..., smooth=TRUE)} shortcut, this includes the
    arguments meant for the other functions (such as \code{\link{plot.roc}}
    or \code{\link{ci}}). Defaults filled in by the smoothers themselves
    are not recorded. The list serves to apply the smoothing again in
    bootstrap operations, \code{\link{roc.test}}, \code{\link{cov}} and
    \code{\link{are.paired}}: the stored values are passed again as values,
    not evaluated again. \code{\link{roc.test}} and \code{\link{cov}}
    compare these lists to detect curves smoothed differently.
  }
  \item{auc}{if the original ROC curve contained an AUC, it is computed
    again on the smoothed ROC.
  }
  \item{ci}{if the original ROC curve contained a CI, it is computed
    again on the smoothed ROC.
  }
  \item{fit.controls, fit.cases}{with  \code{method="fitdistr"} only: 
    the result of \pkg{MASS}'s
    \code{\link[MASS]{fitdistr}} function for controls and cases, with an
    additional \dQuote{densfun} item indicating the density function, if
    possible as character.
  }
  \item{logcondens}{with \code{method="logcondens"} and \code{method="logcondens.smooth"} only: 
    the result of \pkg{logcondens}'s \link[logcondens]{logConROC} function.
  }
  \item{model}{with \code{method="binormal"} only: 
    the linear model from \code{\link{lm}} used to smooth the ROC curve.
  }
  \item{proper}{with \code{method="proper"} only:
    a list with the estimates \code{da} and \code{c}, the
    \code{cutpoints} (category boundaries on Metz and Pan's decision
    variable \eqn{v}), the maximized log-likelihood \code{logLik} (Metz
    and Pan's Eq. 2, without the terms that do not depend on the
    parameters), and the
    \code{categories} (a data frame with the number of \code{controls}
    and \code{cases} in each category, from the least to the most
    \dQuote{positive}). With perfect performance (see
    \sQuote{Warnings}), \code{da} is 0, \code{c} is 1, and the
    \code{cutpoints} and \code{logLik} are \code{NA}.
  }
  \subsection{Attributes}{
    Additionally, the original \code{\link{roc}} object is stored as a
    \dQuote{roc} attribute.
  }
} 

\section{Warnings}{
  With \code{method="proper"}, no curve
  parameters can be estimated from perfectly separated controls and
  cases (Metz and Pan, 1999, p. 27), nor, more generally, from a curve
  whose points all have a sensitivity of 1 or a specificity of 1: the
  likelihood is maximal only in the limit of an AUC of 1. As PROPROC,
  the curve is then this limit, \eqn{d_a = 0} and \eqn{c = 1}, with an
  AUC of 1, and a warning is issued.
}

\section{Errors}{
  The message \dQuote{The 'density' function must return a numeric
    vector or a list with a 'y' item.} will be displayed if the
  \code{density} function did not return a valid output. The message
  \dQuote{Length of 'density.controls' and 'density.cases' differ.}
  will be displayed if the returned value differ in length.

  Binormal smoothing cannot smooth a ROC curve defined by only one
  point, nor degenerate data (Metz and Pan, 1999): a curve whose points
  with sensitivity and specificity strictly between 0 and 1 all have the
  same sensitivity, or all have the same specificity. Any such attempt
  fails with an error. With \code{method="proper"}, degenerate data can
  be smoothed, but a fit that does not converge fails with an error.

  If \code{method="fitdistr"} and \code{\link[MASS]{fitdistr}} fails to
  estimate the parameters of the distribution because the controls or
  the cases take a single distinct value (for example with a
  \dQuote{weibull} or \dQuote{logistic} density), an error is issued.
  With the \dQuote{normal}, \dQuote{lognormal}, \dQuote{exponential},
  \dQuote{poisson} and \dQuote{geometric} densities, \code{fitdistr}
  estimates all the parameters in closed form, so parameters of the
  density (such as \code{mean}) cannot be held fixed: passing them in
  \dots is an error. Use a density function with \code{start} instead.

  If \code{method="logcondens"} or \code{method="logcondens.smooth"} and
  the controls or the cases take fewer than two distinct values, an
  error is issued.

  If the smooth ROC curve was generated by \code{\link{roc}} with
  \code{density.controls} and \code{density.cases} numeric arguments, it
  cannot be smoothed and the error \dQuote{Cannot smooth a ROC curve
    generated directly with numeric 'density.controls' and
    'density.cases'.} is produced.
    
  \code{fitdistr} and \code{density} smoothing methods require a
  \link{numeric} \code{predictor}. If the ROC curve to smooth was
  generated with an ordered factor only binormal and proper smoothing
  can be applied and the message \dQuote{ROC curves of ordered
  predictors can be smoothed only with binormal or proper smoothing.} is
  displayed otherwise.

  \code{fitdistr}, \code{logcondens} and \code{logcondens.smooth} methods,
  and \code{method="proper"}, require additional packages (\pkg{MASS},
  \pkg{logcondens} and \pkg{mvtnorm}, respectively). In an interactive session, you will be offered to install
  a missing package; otherwise an error is issued.
}

\references{
  James E. Hanley (1988) ``The robustness of the ``binormal'' assumptions
  used in fitting ROC curves''. \emph{Medical Decision Making} \bold{8}, 197--203.

  Charles E. Metz, Xiaochuan Pan (1999) ``"Proper" binormal ROC curves:
  theory and maximum-likelihood estimation''. \emph{Journal of
  Mathematical Psychology} \bold{43}, 1--33.
  DOI: \doi{10.1006/jmps.1998.1218}.

  Lorenzo L. Pesce, Charles E. Metz (2007) ``Reliable and computationally
  efficient maximum-likelihood estimation of "proper" binormal ROC
  curves''. \emph{Academic Radiology} \bold{14}, 814--829.
  DOI: \doi{10.1016/j.acra.2007.03.012}.

  Lutz Duembgen, Kaspar Rufibach (2011) ``logcondens: Computations Related
  to Univariate Log-Concave Density Estimation''. \emph{Journal of Statistical
  Software}, \bold{39}, 1--28. URL: \doi{10.18637/jss.v039.i06}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{12}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.

  Kaspar Rufibach (2012) ``A Smooth ROC Curve Estimator Based on Log-Concave Density Estimates''.
  \emph{The International Journal of Biostatistics}, \bold{8}(1), 1--29. DOI: 
  \doi{10.1515/1557-4679.1378}. arXiv:
  \href{https://arxiv.org/abs/1103.1787}{arXiv:1103.1787}.


  William N. Venables, Brian D. Ripley (2002). ``Modern Applied Statistics with S''. New York, Springer.
  \href{http://books.google.ch/books?id=974c4vKurNkC}{Google books}.
  
  Kelly H. Zou, W. J. Hall and David E. Shapiro (1997) ``Smooth
  non-parametric receiver operating characteristic (ROC) curves for
  continuous diagnostic tests''. \emph{Statistics in Medicine}
  \bold{16}(19), 2143--2156. DOI: 
  \doi{10.1002/(SICI)1097-0258(19971015)16:19<2143::AID-SIM655>3.0.CO;2-3}.
}

\seealso{
 \code{\link{roc}}
 
CRAN packages \pkg{MASS} and \pkg{logcondens} employed in this function.
}

\examples{
data(aSAH)
roc1 <- roc(aSAH$outcome, aSAH$s100b)

# Binormal smoothing, the default
smooth(roc1)

# Another method
if (requireNamespace("logcondens", quietly = TRUE)) {
smooth(roc1, method = "logcondens")
}

# Fit a different distribution to controls and cases
if (requireNamespace("MASS", quietly = TRUE)) {
smooth(roc1, method = "fitdistr",
       density.controls = "normal", density.cases = "lognormal")
}

# Plot the smoothed curve over the empirical one
plot(roc1)
plot(smooth(roc1), add = TRUE, col = "blue")
legend("bottomright", legend = c("Empirical", "Smoothed"),
       col = c(par("fg"), "blue"), lwd = 2)
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
\keyword{smooth}
