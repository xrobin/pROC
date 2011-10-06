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
power = NULL, alternative = c("two.sided", "one.sided"), 
reuse.auc=TRUE, ...)
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
  \item{reuse.auc}{if \code{TRUE} (default) and the \dQuote{roc} objects
    contain an \dQuote{auc} field, re-use these specifications for the
    test. See the \emph{AUC specification} section for more details.
  }
  \item{\dots}{further arguments passed to or from other methods,
  especially \code{\link{auc}}.}
}

\section{One ROC curve power calculation}{
  If one or no ROC curves are passed to \code{power.roc.test}, a one ROC
  curve power calculation is performed. The function expects either
  \code{power}, \code{sig.level} or \code{auc}, or both \code{ncontrols}
  and \code{ncases} to be missing, so that the parameter is determined
  from the others with the formula by Obuchowski \emph{et al.}, 2004 (formulas
  2 and 3, p. 1123).

  For the sample size, \code{ncases} is computed directly from formulas
  2 and 3 and ncontrols is deduced with \code{kappa}.
  AUC is optimized by \code{\link{uniroot}} while \code{sig.level}
  and \code{power} are solved as quadratic equations.

  \code{power.roc.test} can also be passed a \code{roc} object from the \code{\link{roc}}
  function, but the empirical ROC will not be used, only the number of
  patients and the AUC.
}

\section{Two paired ROC curves power calculation}{
  If two ROC curves are passed to \code{power.roc.test}, the function
  will compute either the required sample size (if \code{power} is supplied),
  the significance level (if \code{sig.level=NULL} and \code{power} is
  supplied) or the power of a test of a difference between to AUCs
  according to the formula by Obuchowski and McClish, 1997\emph{et al.}
  (formulas 2 to 5, p. 1530--1531). The null hypothesis is that both ROC
  curves have the same AUC than the ROC curve (either \code{roc1} or
  \code{roc2}) having the lowest AUC (thus the highest variance).

  For the sample size, \code{ncases} is computed directly from formula 2
  and ncontrols is deduced from the ratio observed in \code{roc1} and \code{roc2}.
  AUC is optimized by \code{\link{uniroot}} while \code{sig.level}
  and \code{power} are solved as quadratic equations.

  Power calculation for unpaired ROC curves is not implemented.
}

\section{Binormality assumption}{
  Both methods make the assumption that the data is binormal.
  If the data shows a deviation from this assumption, it might help to
  normalize the data first (that is, before calling \code{\link{roc}}),
  for example with quantile normalization: 
  \preformatted{
    norm.x <- qnorm(rank(x)/(length(x)+1))
    smooth(roc(response, norm.x, ...), ...)
  }
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

  Warning: if the roc object passed to roc.test contains an \code{auc}
  field and \code{reuse.auc=TRUE}, \link{auc} is not called and
  arguments such as \code{partial.auc} are silently ignored.
}

\value{
  An object of class \code{power.htest} (such as that given by
  \code{\link{power.t.test}}) with the supplied and computed values.
}

\references{
  Nancy A. Obuchowski, Donna K. McClish (1997). ``Sample size
  determination for diagnostic accurary studies involving binormal ROC
  curve indices''. \emph{Statistics in Medicine}, \bold{16}(13),
  1529--1542. DOI: \href{http://dx.doi.org/10.1002/(SICI)1097-0258(19970715)16:13<1529::AID-SIM565>3.0.CO;2-H}{(SICI)1097-0258(19970715)16:13<1529::AID-SIM565>3.0.CO;2-H}.

  Nancy A. Obuchowski, Micharl L. Lieber, Frank H. Wians
  Jr. (2004). ``ROC Curves in Clinical Chemistry: Uses, Misuses, and
  Possible Solutions''. \emph{Clin Chem}, \bold{50}(7), 1118--1125. DOI:
  \href{http://dx.doi.org/10.1373/clinchem.2004.031823}{10.1373/clinchem.2004.031823}.
}

\seealso{
 \code{\link{roc}}, \code{\link{roc.test}}
}
\examples{
data(aSAH)

#### One ROC curve ####

# Build a roc object:
rocobj <- roc(aSAH$outcome, aSAH$s100b)

# Determine power of one ROC curve:
power.roc.test(rocobj)
# Same as:
power.roc.test(ncases=41, ncontrols=72, auc=0.73, sig.level=0.05)
# sig.level=0.05 is implicit and can be omitted:
power.roc.test(ncases=41, ncontrols=72, auc=0.73)

# Determine ncases & ncontrols:
power.roc.test(auc=rocobj$auc, sig.level=0.05, power=0.95, kappa=1.7)
power.roc.test(auc=0.73, sig.level=0.05, power=0.95, kappa=1.7)

# Determine sig.level:
power.roc.test(ncases=41, ncontrols=72, auc=0.73, power=0.95, sig.level=NULL)

# Derermine detectable AUC:
power.roc.test(ncases=41, ncontrols=72, sig.level=0.05, power=0.95)


#### Two ROC curves ####

##  Full AUC
roc1 <- roc(aSAH$outcome, aSAH$ndka)
roc2 <- roc(aSAH$outcome, aSAH$wfns)
# Sample size
power.roc.test(roc1, roc2, power=0.9)
# Power test
power.roc.test(roc1, roc2)
# Significance level
power.roc.test(roc1, roc2, power=0.9, sig.level=NULL)

# Partial AUC
roc3 <- roc(aSAH$outcome, aSAH$ndka, partial.auc=c(1, 0.9))
roc4 <- roc(aSAH$outcome, aSAH$wfns, partial.auc=c(1, 0.9))

# Sample size
power.roc.test(roc3, roc4, power=0.9)

# Power test
power.roc.test(roc3, roc4)
# This is exactly equivalent:
power.roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(1, 0.9))

# Significance level
power.roc.test(roc3, roc4, power=0.9, sig.level=NULL)

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
