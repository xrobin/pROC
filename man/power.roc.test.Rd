\encoding{UTF-8}
\name{power.roc.test}
\alias{power.roc.test}
\alias{power.roc.test.numeric}
\alias{power.roc.test.roc}
\alias{power.roc.test.list}

\title{
  Sample size and power computation for ROC curves
}
\description{
  Computes sample size, power, significance level or minimum AUC for ROC curves.
}
\usage{
power.roc.test(...)
# One or Two ROC curves test with roc objects:
\S3method{power.roc.test}{roc}(roc1, roc2, sig.level = 0.05, 
power = NULL, kappa = NULL, alternative = c("two.sided", "one.sided"),
reuse.auc=TRUE, method = c("delong", "bootstrap", "obuchowski"), ...)
# One ROC curve with a given AUC:
\S3method{power.roc.test}{numeric}(auc = NULL, ncontrols = NULL, 
ncases = NULL, sig.level = 0.05, power = NULL, kappa = 1, 
alternative = c("two.sided", "one.sided"), ...)
# Two ROC curves with the given parameters:
\S3method{power.roc.test}{list}(parslist, ncontrols = NULL, 
ncases = NULL, sig.level = 0.05, power = NULL,  kappa = 1, 
alternative = c("two.sided", "one.sided"), ...)
}

\arguments{
  \item{roc1, roc2}{one or two \dQuote{roc} object from the
	\code{\link{roc}} function.
  }
  \item{auc}{
    expected AUC.
  }
  \item{parslist}{
    a \code{\link{list}} of parameters for the two ROC curves test with
    Obuchowski variance when no empirical ROC curve is known:
    \describe{
      \item{A1}{binormal A parameter for ROC curve 1}
      \item{B1}{binormal B parameter for ROC curve 1}
      \item{A2}{binormal A parameter for ROC curve 2}
      \item{B2}{binormal B parameter for ROC curve 2}
      \item{rn}{correlation between the variables in control patients}
      \item{ra}{correlation between the variables in case patients}
      \item{delta}{the difference of AUC between the two ROC curves}
    }
    For a partial AUC, the following additional parameters must be set:
    \describe{
      \item{FPR11}{Upper bound of FPR (1 - specificity) of ROC curve 1}
      \item{FPR12}{Lower bound of FPR (1 - specificity) of ROC curve 1}
      \item{FPR21}{Upper bound of FPR (1 - specificity) of ROC curve 2}
      \item{FPR22}{Lower bound of FPR (1 - specificity) of ROC curve 2}
    }
  }
  \item{ncontrols, ncases}{
    number of controls and case observations available.
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
  \item{method}{the method to compute \link[=var.roc]{variance} and
    \link[=cov.roc]{covariance}, either \dQuote{delong},
    \dQuote{bootstrap} or \dQuote{obuchowski}. The first letter is
    sufficient. Only for Two ROC curves power calculation. See
    \code{\link{var}} and \code{\link{cov}} documentations for more
    details. 
  }
  \item{\dots}{further arguments passed to or from other methods,
    especially \code{\link{auc}} (with \code{reuse.auc=FALSE} or no AUC in
    the ROC curve), \code{\link{cov}} and \code{\link{var}} (especially
    arguments \code{method}, \code{boot.n} and \code{boot.stratified}).
    Ignored (with a warning) with a \code{parslist}.
  }
}

\section{One ROC curve power calculation}{
  If one or no ROC curves are passed to \code{power.roc.test}, a one ROC
  curve power calculation is performed. The function expects either
  \code{power}, \code{sig.level} or \code{auc}, or both \code{ncontrols}
  and \code{ncases} to be missing, so that the parameter is determined
  from the others with the formula by Obuchowski \emph{et al.}, 2004 (formulas
  2 and 3, p. 1123).

  For the sample size, \code{ncases} is computed directly from formulas
  2 and 3 and ncontrols is deduced with \code{kappa} (defaults to the
  ratio of controls to cases).
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
  according to the formula by Obuchowski and McClish, 1997
  (formulas 2 and 3, p. 1530--1531). The null hypothesis is that the AUC
  of \code{roc1} is the same than the AUC of \code{roc2}, with
  \code{roc1} taken as the reference ROC curve.

  For the sample size, \code{ncases} is computed directly from formula 2
  and ncontrols is deduced with \code{kappa} (defaults to the
  ratio of controls to cases in \code{roc1}).
  \code{sig.level} and \code{power} are solved as quadratic equations.

  The variance and covariance of the ROC curve are computed with the
  \code{\link{var}} and \code{\link{cov}} functions. By default, DeLong
  method using the algorithm by Sun and Xu (2014) is used for full
  AUCs and the bootstrap for partial AUCs. It is
  possible to force the use of Obuchowski's variance by specifying
  \code{method="obuchowski"}.

  Alternatively when no empirical ROC curve is known, or if only one is
  available, a list can be passed to \code{power.roc.test}, with the
  contents defined in the \dQuote{Arguments} section. The variance and
  covariance are computed from Table 1 and Equation 4 and 5 of
  Obuchowski and McClish (1997), p. 1530--1531.

  Power calculation for unpaired ROC curves is not implemented.
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

\section{Acknowledgements}{
  The authors would like to thank Christophe Combescure and Anne-Sophie
  Jannot for their help with the implementation of this section of the package.
}

\references{
  Elisabeth R. DeLong, David M. DeLong and Daniel L. Clarke-Pearson
  (1988) ``Comparing the areas under two or more correlated receiver
  operating characteristic curves: a nonparametric
  approach''. \emph{Biometrics} \bold{44}, 837--845.

  Nancy A. Obuchowski, Donna K. McClish (1997). ``Sample size
  determination for diagnostic accurary studies involving binormal ROC
  curve indices''. \emph{Statistics in Medicine}, \bold{16},
  1529--1542. DOI: \doi{10.1002/(SICI)1097-0258(19970715)16:13<1529::AID-SIM565>3.0.CO;2-H}.

  Nancy A. Obuchowski, Micharl L. Lieber, Frank H. Wians
  Jr. (2004). ``ROC Curves in Clinical Chemistry: Uses, Misuses, and
  Possible Solutions''. \emph{Clinical Chemistry}, \bold{50}, 1118--1125. DOI:
  \doi{10.1373/clinchem.2004.031823}.
  
  Xu Sun and Weichao Xu (2014) ``Fast Implementation of DeLongs Algorithm for Comparing
  the Areas Under Correlated Receiver Operating Characteristic Curves''. \emph{IEEE Signal
  Processing Letters}, \bold{21}, 1389--1393. 
  DOI: \doi{10.1109/LSP.2014.2337313}.
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

###  Full AUC
roc1 <- roc(aSAH$outcome, aSAH$ndka)
roc2 <- roc(aSAH$outcome, aSAH$wfns)

## Sample size
# With DeLong variance (default)
power.roc.test(roc1, roc2, power=0.9)
# With Obuchowski variance
power.roc.test(roc1, roc2, power=0.9, method="obuchowski")

## Power test
# With DeLong variance (default)
power.roc.test(roc1, roc2)
# With Obuchowski variance
power.roc.test(roc1, roc2, method="obuchowski")

## Significance level
# With DeLong variance (default)
power.roc.test(roc1, roc2, power=0.9, sig.level=NULL)
# With Obuchowski variance
power.roc.test(roc1, roc2, power=0.9, sig.level=NULL, method="obuchowski")

### Partial AUC
roc3 <- roc(aSAH$outcome, aSAH$ndka, partial.auc=c(1, 0.9))
roc4 <- roc(aSAH$outcome, aSAH$wfns, partial.auc=c(1, 0.9))

## Sample size
# With bootstrap variance (default)
\dontrun{
power.roc.test(roc3, roc4, power=0.9)
}
# With Obuchowski variance
power.roc.test(roc3, roc4, power=0.9, method="obuchowski")

## Power test
# With bootstrap variance (default)
\dontrun{
power.roc.test(roc3, roc4)
# This is exactly equivalent:
power.roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(1, 0.9))
}
# With Obuchowski variance
power.roc.test(roc3, roc4, method="obuchowski")

## Significance level
# With bootstrap variance (default)
\dontrun{
power.roc.test(roc3, roc4, power=0.9, sig.level=NULL)
}
# With Obuchowski variance
power.roc.test(roc3, roc4, power=0.9, sig.level=NULL, method="obuchowski")

## With only binormal parameters given
# From example 2 of Obuchowski and McClish, 1997.
ob.params <- list(A1=2.6, B1=1, A2=1.9, B2=1, rn=0.6, ra=0.6, FPR11=0,
FPR12=0.2, FPR21=0, FPR22=0.2, delta=0.037) 

power.roc.test(ob.params, power=0.8, sig.level=0.05)
power.roc.test(ob.params, power=0.8, sig.level=NULL, ncases=107)
power.roc.test(ob.params, power=NULL, sig.level=0.05, ncases=107)

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
