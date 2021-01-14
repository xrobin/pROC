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
"logcondens.smooth"), n=512, bw = "nrd0", density=NULL,
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
    \dQuote{logcondens}, \dQuote{"logcondens.smooth"}.
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
    \code{\link{density}}. If \code{method="fitdistr"}, a \code{densfun}
    argument for \code{\link[MASS]{fitdistr}}.
    If the value is different for control and case observations,
    \code{density.controls} and \code{density.cases} can be employed
    instead, otherwise \code{density} will be propagated to both
	\code{density.controls} and \code{density.cases}.
  }
  \item{start, start.controls, start.cases}{if
    \code{method="fitdistr"}, optionnal \code{start}
    arguments for . \code{start.controls}
    and \code{start.cases} allows to specify different distributions for
    controls and cases.
  }
  \item{reuse.auc, reuse.ci}{if \code{TRUE} (default for reuse.auc) and the \dQuote{roc} objects
    contain \dQuote{auc} or \dQuote{ci} fields, re-use these
    specifications to regenerate \code{\link{auc}} or \code{\link{ci}}
    on the smoothed ROC curve with the original parameters. If
    \code{FALSE}, the object returned will not contain
    \dQuote{auc} or \dQuote{ci} fields. It is currently not possible to
    redefine \code{\link{auc}} and \code{\link{ci}} options directly: you need to call \code{\link{auc}} or
    \code{\link{ci}} later for that.
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

  With \code{method="density"}, the \code{\link{density}}
  function is employed to generate a smooth kernel
  density of the control and case observations as described by Zhou
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
  kernels are symetrical, it might help to normalize the data first
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
  the density function \code{density} with optionnal start parameters
  \code{start}. The density function are fitted
  separately in control (\code{density.controls}, \code{start.controls})
  and case observations (\code{density.cases},
  \code{start.cases}). \code{density} can be one of the character values
  allowed by \code{\link[MASS]{fitdistr}} or a density function (such
  as \code{\link{dnorm}}, \code{\link{dweibull}}, ...).

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
    reported in percent, as defined in argument.
  }
  \item{direction}{the direction of the comparison, as defined in argument.}
  \item{call}{how the function was called. See \code{\link{match.call}} for
    more details.
  }
  \item{smoothing.args}{a list of the arguments used for the
    smoothing. Will serve to apply the smoothing again in further
    bootstrap operations.
  }
  \item{auc}{if the original ROC curve contained an AUC, it is computed
    again on the smoothed ROC.
  }
  \item{ci}{if the original ROC curve contained a CI, it is computed
    again on the smoothed ROC.
  }
  \item{fit.controls, fit.cases}{with  \code{method="fitdistr"} only: 
    the result of \pkg{MASS}'s
    \code{\link{fitdistr}} function for controls and cases, with an
    additional \dQuote{densfun} item indicating the density function, if
    possible as character.
  }
  \item{logcondens}{with \code{method="logcondens"} and \code{method="logcondens.smooth"} only: 
    the result of \pkg{logcondens}'s \link[logcondens]{logConROC} function.
  }
  \item{model}{with \code{method="binormal"} only: 
    the linear model from \code{\link{lm}} used to smooth the ROC curve.
  }
  \subsection{Attributes}{
    Additionally, the original \code{\link{roc}} object is stored as a
    \dQuote{roc} attribute.
  }
} 

\section{Errors}{
  The message \dQuote{The 'density' function must return a numeric
    vector or a list with a 'y' item.} will be displayed if the
  \code{density} function did not return a valid output. The message
  \dQuote{Length of 'density.controls' and 'density.cases' differ.}
  will be displayed if the returned value differ in length.

  Binormal smoothing cannot smooth ROC curve defined by only one
  point. Any such attempt will fail with the error \dQuote{ROC curve not
    smoothable (not enough points).}.

  If the smooth ROC curve was generated by \code{\link{roc}} with
  \code{density.controls} and \code{density.cases} numeric arguments, it
  cannot be smoothed and the error \dQuote{Cannot smooth a ROC curve
    generated directly with numeric 'density.controls' and
    'density.cases'.} is produced.
    
  \code{fitdistr} and \code{density} smoothing methods require a
  \link{numeric} \code{predictor}. If the ROC curve to smooth was
  generated with an ordered factor only binormal smoothing can be
  applied and the message \dQuote{ROC curves of ordered predictors can
  be smoothed only with binormal smoothing.} is displayed otherwise.

  \code{fitdistr}, \code{logcondens} and \code{logcondens.smooth} methods
  require additional packages. If not available, the following message 
  will be displayed with the required command to install the package:
  \dQuote{Package ? not available, required with method='?'.
     Please install it with 'install.packages("?")'. 
  }
}

\references{
  James E. Hanley (1988) ``The robustness of the ``binormal'' assumptions
  used in fitting ROC curves''. \emph{Medical Decision Making} \bold{8}, 197--203.

  Lutz Duembgen, Kaspar Rufibach (2011) ``logcondens: Computations Related
  to Univariate Log-Concave Density Estimation''. \emph{Journal of Statistical
  Software}, \bold{39}, 1--28. URL: \href{https://www.jstatsoft.org/v39/i06/}{jstatsoft.org/v39/i06}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.

  Kaspar Rufibach (2011) ``A Smooth ROC Curve Estimator Based on Log-Concave Density Estimates''.
  \emph{The International Journal of Biostatistics}, \bold{8}, accepted. DOI: 
  \doi{10.1515/1557-4679.1378}. arXiv:
  \href{https://arxiv.org/abs/1103.1787}{1103.1787}.

  William N. Venables, Brian D. Ripley (2002). ``Modern Applied Statistics with S''. New York, Springer.
  \href{http://books.google.ch/books?id=974c4vKurNkC}{Google books}.
  
  Kelly H. Zou, W. J. Hall and David E. Shapiro (1997) ``Smooth
  non-parametric receiver operating characteristic (ROC) curves for
  continuous diagnostic tests''. \emph{Statistics in Medicine}
  \bold{18}, 2143--2156. DOI: 
  \doi{10.1002/(SICI)1097-0258(19971015)16:19<2143::AID-SIM655>3.0.CO;2-3}.
}

\seealso{
 \code{\link{roc}}
 
CRAN packages \pkg{MASS} and \pkg{logcondens} employed in this function.
}

\examples{
data(aSAH)

##  Basic example
rocobj <- roc(aSAH$outcome, aSAH$s100b)
smooth(rocobj)
# or directly with roc()
roc(aSAH$outcome, aSAH$s100b, smooth=TRUE)

# plotting
plot(rocobj)
rs <- smooth(rocobj, method="binormal")
plot(rs, add=TRUE, col="green")
rs2 <- smooth(rocobj, method="density")
plot(rs2, add=TRUE, col="blue")
rs3 <- smooth(rocobj, method="fitdistr", density="lognormal")
plot(rs3, add=TRUE, col="magenta")
if (requireNamespace("logcondens")) {
rs4 <- smooth(rocobj, method="logcondens")
plot(rs4, add=TRUE, col="brown")
rs5 <- smooth(rocobj, method="logcondens.smooth")
plot(rs5, add=TRUE, col="orange")
}
legend("bottomright", legend=c("Empirical", "Binormal", "Density", "Log-normal",
                               "Log-concave density", "Smoothed log-concave density"),
       col=c("black", "green", "blue", "magenta", "brown", "orange"), lwd=2)

## Advanced smoothing

# if we know the distributions are normal with sd=0.1 and an unknown mean:
smooth(rocobj, method="fitdistr", density=dnorm, start=list(mean=1), sd=.1)
# different distibutions for controls and cases:
smooth(rocobj, method="fitdistr", density.controls="normal", density.cases="lognormal")

# with densities
bw <- bw.nrd0(rocobj$predictor)
density.controls <- density(rocobj$controls, from=min(rocobj$predictor) - 3 * bw,
                            to=max(rocobj$predictor) + 3*bw, bw=bw, kernel="gaussian")
density.cases <- density(rocobj$cases, from=min(rocobj$predictor) - 3 * bw,
                            to=max(rocobj$predictor) + 3*bw, bw=bw, kernel="gaussian")
smooth(rocobj, method="density", density.controls=density.controls$y, 
       density.cases=density.cases$y)
# which is roughly what is done by a simple:
smooth(rocobj, method="density")

\dontrun{
## Smoothing artificial ROC curves

rand.unif <- runif(1000, -1, 1)
rand.exp <- rexp(1000)
rand.norm <- 
rnorm(1000)

# two normals
roc.norm <- roc(controls=rnorm(1000), cases=rnorm(1000)+1, plot=TRUE)
plot(smooth(roc.norm), col="green", lwd=1, add=TRUE)
plot(smooth(roc.norm, method="density"), col="red", lwd=1, add=TRUE)
plot(smooth(roc.norm, method="fitdistr"), col="blue", lwd=1, add=TRUE)
if (requireNamespace("logcondens")) {
plot(smooth(roc.norm, method="logcondens"), col="brown", lwd=1, add=TRUE)
plot(smooth(roc.norm, method="logcondens.smooth"), col="orange", lwd=1, add=TRUE)
}
legend("bottomright", legend=c("empirical", "binormal", "density", "fitdistr",
                               "logcondens", "logcondens.smooth"), 
       col=c(par("fg"), "green", "red", "blue", "brown", "orange"), lwd=c(2, 1, 1, 1))
       
# deviation from the normality
roc.norm.exp <- roc(controls=rnorm(1000), cases=rexp(1000), plot=TRUE)
plot(smooth(roc.norm.exp), col="green", lwd=1, add=TRUE)
plot(smooth(roc.norm.exp, method="density"), col="red", lwd=1, add=TRUE)
# Wrong fitdistr: normality assumed by default
plot(smooth(roc.norm.exp, method="fitdistr"), col="blue", lwd=1, add=TRUE)
# Correct fitdistr
plot(smooth(roc.norm.exp, method="fitdistr", density.controls="normal",
            density.cases="exponential"), col="purple", lwd=1, add=TRUE)
if (requireNamespace("logcondens")) {
plot(smooth(roc.norm.exp, method="logcondens"), col="brown", lwd=1, add=TRUE)
plot(smooth(roc.norm.exp, method="logcondens.smooth"), col="orange", lwd=1, add=TRUE)
}
legend("bottomright", legend=c("empirical", "binormal", "density",
                               "wrong fitdistr", "correct fitdistr",
                               "logcondens", "logcondens.smooth"),
       col=c(par("fg"), "green", "red", "blue", "purple", "brown", "orange"), lwd=c(2, 1, 1, 1, 1))

# large deviation from the normality
roc.unif.exp <- roc(controls=runif(1000, 2, 3), cases=rexp(1000)+2, plot=TRUE)
plot(smooth(roc.unif.exp), col="green", lwd=1, add=TRUE)
plot(smooth(roc.unif.exp, method="density"), col="red", lwd=1, add=TRUE)
plot(smooth(roc.unif.exp, method="density", bw="ucv"), col="magenta", lwd=1, add=TRUE)
# Wrong fitdistr: normality assumed by default (uniform distributions not handled)
plot(smooth(roc.unif.exp, method="fitdistr"), col="blue", lwd=1, add=TRUE)
if (requireNamespace("logcondens")) {
plot(smooth(roc.unif.exp, method="logcondens"), col="brown", lwd=1, add=TRUE)
plot(smooth(roc.unif.exp, method="logcondens.smooth"), col="orange", lwd=1, add=TRUE)
}
legend("bottomright", legend=c("empirical", "binormal", "density",
                               "density ucv", "wrong fitdistr",
                               "logcondens", "logcondens.smooth"),
       col=c(par("fg"), "green", "red", "magenta", "blue", "brown", "orange"), lwd=c(2, 1, 1, 1, 1))
}

# 2 uniform distributions with a custom density function
unif.density <- function(x, n, from, to, bw, kernel, ...) {
  smooth.x <- seq(from=from, to=to, length.out=n)
  smooth.y <- dunif(smooth.x, min=min(x), max=max(x))
  return(smooth.y)
}
roc.unif <- roc(controls=runif(1000, -1, 1), cases=runif(1000, 0, 2), plot=TRUE)
s <- smooth(roc.unif, method="density", density=unif.density)
plot(roc.unif)
plot(s, add=TRUE, col="grey")

\dontrun{
# you can bootstrap a ROC curve smoothed with a density function:
ci(s, boot.n=100)
}
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
\keyword{smooth}
