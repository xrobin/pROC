\encoding{UTF-8}
\name{power.roc.test}
\alias{power.roc.test}
\alias{power.roc.test.numeric}
\alias{power.roc.test.roc}

\title{
  Sample size and power computation for ROC curves
}
\description{
  Computes sample size, power, significance level or minimum AUC for ROC curves.
}
\usage{
power.roc.test(...)
\S3method{power.roc.test}{roc}(roc1, roc2, sig.level = 0.05, 
power = NULL, alternative = c("two.sided", "one.sided"), ...)
\S3method{power.roc.test}{numeric}(ncontrols = NULL, ncases = NULL, 
auc = NULL, sig.level = 0.05, power = NULL, kappa = 1, 
alternative = c("two.sided", "one.sided"), ...)
}

\arguments{
  \item{roc1, roc2}{one or two \dQuote{roc} object from the
	\code{\link{roc}} function.
  }
  \item{ncontrols, ncases}{
    number of controls and case observations available.
  }
  \item{auc}{
    expected AUC.
  }
  \item{sig.level}{expected significance level (probability of type I
    error).
  }
  \item{power}{expected power of the test (1 - probability of type II
    error).
  }
  \item{kappa}{
    expected balance between control and case observations. Must be
    positive. Only for sample size determination, that is to determine
    \code{ncontrols} and \code{ncases}.
  }
  \item{alternative}{whether a one or two-sided test is performed.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\section{Two ROC curves power calculation}{
  If two ROC curves are passed to \code{power.roc.test}, the function
  will compute either the required sample size (if both sig.level and
  power are supplied), the power or the significance level (if
  \code{sig.level=NULL}) of a test of a difference between to AUCs
  according to the formula by Obuchowski \emph{et al.} (formulas 2 and
  3, p. 1123).
}

\section{One ROC curve power calculation}{
  If one or no ROC curves are passed to \code{power.roc.test}, a one ROC
  curve power calculation is performed. The function expects either
  \code{power}, \code{sig.level} or \code{auc}, or both \code{ncontrols}
  and \code{ncases} to be missing, so that the parameter is determined
  from the others with the formula by Obuchowski \emph{et al.} (formulas
  2 and 3, p. 1123).

  For the sample size, \code{ncases} is computed directly from formulas
  2 and 3 and ncontrols is deduced with \code{kappa}.
  AUC is optimized by \code{\link{uniroot}} while \code{sig.level}
  and \code{power} are solved as quadratic equations.

  \code{power.roc.test} can also be passed a \code{roc} object from the \code{\link{roc}}
  function, but the empirical ROC will not be used, only the number of
  patients and the AUC.
}

\value{
  An object of class \code{power.htest} (such as that given by
  \code{\link{power.t.test}}) with the supplied and computed values.
}

\references{
  Nancy A. Obuchowski, Micharl L. Lieber, Frank H. Wians
  Jr. (2004). ``ROC Curves in Clinical Chemistry: Uses, Misuses, and
  Possible Solutions''. \emph{Clin Chem}, \bold{50}(7), 1118--1125. DOI:
  \href{http://dx.doi.org/10.1373/clinchem.2004.031823}{10.1373/clinchem.2004.031823}.
}

\seealso{
\code{\link{roc}}
}
\examples{
data(aSAH)

# Print a roc object:
rocobj <- roc(aSAH$outcome, aSAH$s100b)
# Determine power of a ROC curve:
power.roc.test(rocobj)
# Same:
power.roc.test(ncases=41, ncontrols=72, auc=0.73, sig.level=0.05)
# sig.level=0.05 is implicit and can be omitted:
power.roc.test(ncases=41, ncontrols=72, auc=0.73)
# Determine ncases & ncontrols:
power.roc.test(auc=0.73, sig.level=0.05, power=0.95, kappa=1.7)
# Determine sig.level:
power.roc.test(ncases=41, ncontrols=72, auc=0.73, power=0.95)
# Derermine detectable AUC:
power.roc.test(ncases=41, ncontrols=72, sig.level=0.05, power=0.95)
}


\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
