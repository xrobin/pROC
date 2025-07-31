\encoding{UTF-8}
\name{pROC-package}
\alias{pROC-package}
\alias{pROC}
\docType{package}
\title{pROC}

\description{
  Tools for visualizing, smoothing and comparing receiver operating
  characteristic (ROC curves). (Partial) area under the curve (AUC) can
  be compared with statistical tests based on U-statistics or
  bootstrap. Confidence intervals can be computed for (p)AUC or ROC
  curves. Sample size / power computation for one or two ROC curves are available.
}

\details{
  The basic unit of the pROC package is the \code{\link{roc}} function. It
  will build a ROC curve, smooth it if requested (if \code{smooth=TRUE}),
  compute the AUC (if \code{auc=TRUE}), the confidence interval (CI) if 
  requested (if \code{ci=TRUE}) and plot the curve if requested (if
  \code{plot=TRUE}). 

  The \code{\link{roc}} function will call \code{\link[=smooth.roc]{smooth}},
  \code{\link{auc}},
  \code{\link{ci}} and \code{\link{plot}} as necessary. See these
  individual functions for the arguments that can be passed to them
  through \code{\link{roc}}. These function can be called separately.

  Two paired (that is \code{\link{roc}} objects with the same
  \code{response}) or unpaired (with different \code{response}) ROC
  curves can be compared with the \code{\link{roc.test}} function.
}

\section{Citation}{
  If you use pROC in published research, please cite the following paper:

  Xavier Robin, Natacha Turck, Alexandre Hainard, Natalia Tiberti,
  Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller (2011).
  ``pROC: an open-source package for R and S+ to analyze and compare ROC
  curves''. \emph{BMC Bioinformatics}, \bold{12}, p. 77. DOI: 
  \doi{10.1186/1471-2105-12-77}

  Type \code{citation("pROC")} for a BibTeX entry.

  The authors would be glad to hear how pROC is employed. You are kindly
  encouraged to notify Xavier Robin <pROC-cran@xavier.robin.name>
  about any work you publish.
}

\section{Abbreviations}{
  The following abbreviations are employed extensively in this package:
  \itemize{
    \item ROC: receiver operating characteristic
    \item AUC: area under the ROC curve
    \item pAUC: partial area under the ROC curve
    \item CI: confidence interval
    \item SP: specificity
    \item SE: sensitivity
  }
}

\section{Functions}{
  \tabular{ll}{
    \code{\link{roc}} \tab Build a ROC curve\cr
    \code{\link{are.paired}} \tab Dertermine if two ROC curves are paired \cr
    \code{\link{auc}} \tab Compute the area under the ROC curve \cr
    \code{\link{ci}} \tab Compute confidence intervals of a ROC curve \cr
    \code{\link{ci.auc}} \tab Compute the CI of the AUC \cr
    \code{\link{ci.coords}} \tab Compute the CI of arbitrary coordinates \cr
    \code{\link{ci.se}} \tab Compute the CI of sensitivities at given specificities \cr
    \code{\link{ci.sp}} \tab Compute the CI of specificities at given sensitivities \cr
    \code{\link{ci.thresholds}} \tab Compute the CI of specificity and sensitivity of thresholds \cr
    \code{\link{ci.coords}} \tab Compute the CI of arbitrary coordinates \cr
    \code{\link{coords}} \tab Coordinates of the ROC curve \cr
    \code{\link[=cov.roc]{cov}} \tab Covariance between two AUCs\cr
    \code{\link{ggroc}} \tab Plot a ROC curve with \pkg{ggplot2}\cr
    \code{\link{has.partial.auc}} \tab Determine if the ROC curve have a partial AUC\cr
    \code{\link{lines.roc}} \tab Add a ROC line to a ROC plot \cr
    \code{\link{plot.ci}} \tab Plot CIs \cr
    \code{\link[=plot.roc]{plot}} \tab Plot a ROC curve \cr
    \code{\link{power.roc.test}} \tab Sample size and power computation \cr
    \code{\link[=print.roc]{print}} \tab Print a ROC curve object \cr
    \code{\link{roc.test}} \tab Compare two ROC curves \cr
    \code{\link[=smooth.roc]{smooth}} \tab Smooth a ROC curve \cr
    \code{\link[=var.roc]{var}} \tab Variance of the AUC
  }
}

\section{Dataset}{
  This package comes with a dataset of 141 patients with aneurysmal
  subarachnoid hemorrhage: \code{\link{aSAH}}.
}

\section{Installing and using}{
  To install this package, make sure you are connected to the internet and issue the following command in the R prompt:
  
  \preformatted{
    install.packages("pROC")
  }
  
  To load the package in R:
 
  \preformatted{
    library(pROC)
  }
}

\section{Experimental: pipelines}{
	Since version 1.15.0, the \code{\link{roc}} function can be used in pipelines, for instance with \pkg{dplyr} or \pkg{magrittr}. This is still a highly experimental feature and will change significantly in future versions  (see \href{https://github.com/xrobin/pROC/issues/54}{issue 54}).
	The \code{\link{roc.data.frame}} method supports both standard and non-standard evaluation (NSE), and the \code{\link{roc_}}
	function supports standard evaluation only.
	
    \preformatted{
library(dplyr)
aSAH \%>\% 
    filter(gender == "Female") \%>\%
    roc(outcome, s100b)
	}
	
	By default it returns the \code{\link{roc}} object, which can then be piped to
	the \code{\link{coords}} function to extract coordinates that can be used
	in further pipelines.
	
	\preformatted{
aSAH \%>\% 
    filter(gender == "Female") \%>\%
    roc(outcome, s100b) \%>\%
    coords(transpose=FALSE) \%>\%
    filter(sensitivity > 0.6, 
           specificity > 0.6)
	}
	More details and use cases are available in the \code{\link{roc}} help page.
}

\section{Bootstrap}{

  All the bootstrap operations for \link[=roc.test]{significance testing}, \link[=ci]{confidence interval}, \link[=var]{variance} and \link[=cov]{covariance} computation are performed with non-parametric stratified or non-stratified resampling (according to the \code{stratified} argument) and with the percentile method, as described in Carpenter and Bithell (2000) sections 2.1 and 3.3.

  Stratification of bootstrap can be controlled
  with \code{boot.stratified}. In stratified bootstrap (the default), each replicate
  contains the same number of cases and controls than the original
  sample. Stratification is especially useful if one group has only
  little observations, or if groups are not balanced.

  The number of bootstrap replicates is controlled by \code{boot.n}. Higher numbers will give a more precise estimate of the significance tests and confidence intervals
  but take more time to compute. 2000 is recommanded by Carpenter and Bithell (2000) for confidence intervals. In our experience this is sufficient for a good estimation of the 
  first significant digit only, so we recommend the use of 10000 bootstrap replicates to obtain a good estimate of the second significant digit whenever possible.

	\subsection{Progress bars}{
	  Progress bars are no longer supported in pROC 1.19. The use of the
	  \code{progress} argument is deprecated, and any value other than \code{NULL}
	  will produce a warning. The argument will be removed in a future version.
	}
}

\section{Handling large datasets}{
	\subsection{Algorithms}{
	Over the years, a significant amount of time has been invested in making pROC run faster and faster.
	From the naive algorithm iterating over all thresholds implemented in the first version, we went to a
	C++ implementation (with \pkg{Rcpp}), and a different algorithm using cummulative sum of responses sorted 
	by the predictor, which scales only with the number of data points, independently on the number of thresholds (\code{algorithm = 2}).
	The curves themselves are identical, but computation time has been decreased massively.
	
	Earlier versions of \pkg{pROC} allowed choosing the optimal algorithm. Since pROC 1.19, only a single
	algorithm (previously \code{algorithm = 2}) is available. This allowed us to remove legacy code
	and dependencies, and make \pkg{pROC} more maintainable and future proof.
	}
	
    \subsection{Boostrap}{
    Bootstrap is typically slow because it involves repeatedly computing the ROC curve (or a part of it).
    
    Some bootstrap functions are faster than others. Typically, \code{\link{ci.thresholds}} is the fastest, and \code{\link{ci.coords}} the slowest. Use \code{\link{ci.coords}} only if the CI you need cannot be computed by the specialized CI functions \code{\link{ci.thresholds}}, \code{\link{ci.se}} and \code{\link{ci.sp}}. Note that \code{\link{ci.auc}} cannot be replaced anyway.
    
    A naive way to speed-up the boostrap is by removing the progress bar:
    
    \preformatted{
rocobj <- roc(response, round(predictor))
system.time(ci(rocobj))
system.time(ci(rocobj, progress = "none"))
    }
    
    It is of course possible to reduce the number of boostrap iterations. See the \code{boot.n} argument to \code{\link{ci}}. This will reduce the precision of the bootstrap estimate.

    \subsection{Parallel processing}{
	    Parallel processing is no longer supported in pROC 1.19. The use of the
	    \code{parallel} argument is deprecated, and any value other than
	    \code{FALSE} will produce a warning. The argument will be removed in a
	    future version.
  	}
  }
  
  \subsection{Using DeLong instead of boostrap}{
    DeLong is an asymptotically exact method to evaluate the uncertainty of an AUC (DeLong \emph{et al.} (1988)). Since version 1.9, pROC uses the algorithm proposed by Sun and Xu (2014) which has an O(N log N) complexity and is always faster than bootstrapping. By default, pROC will choose the DeLong method whenever possible.
  
        \preformatted{
rocobj <- roc(response, round(predictor), algorithm=3)
system.time(ci(rocobj, method="delong"))
    }
  }
}

\author{
Xavier Robin, Natacha Turck, Jean-Charles Sanchez and Markus Müller

Maintainer: Xavier Robin <pROC-cran@xavier.robin.name>
}

\references{
  James Carpenter and John Bithell (2000) ``Bootstrap condence intervals:
  when, which, what? A practical guide for medical statisticians''.
  \emph{Statistics in Medicine} \bold{19}, 1141--1164.
  DOI: \doi{10.1002/(SICI)1097-0258(20000515)19:9<1141::AID-SIM479>3.0.CO;2-F}.

  Elisabeth R. DeLong, David M. DeLong and Daniel L. Clarke-Pearson
  (1988) ``Comparing the areas under two or more correlated receiver
  operating characteristic curves: a nonparametric
  approach''. \emph{Biometrics} \bold{44}, 837--845.

  Tom Fawcett (2006) ``An introduction to ROC analysis''. \emph{Pattern
    Recognition Letters} \bold{27}, 861--874. DOI:
  \doi{10.1016/j.patrec.2005.10.010}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
  
  Xu Sun and Weichao Xu (2014) ``Fast Implementation of DeLongs Algorithm for Comparing
  the Areas Under Correlated Receiver Operating Characteristic Curves''. \emph{IEEE Signal
  Processing Letters}, \bold{21}, 1389--1393. 
  DOI: \doi{10.1109/LSP.2014.2337313}.
  
  Hadley Wickham (2011) ``The Split-Apply-Combine Strategy for Data Analysis''. \emph{Journal of Statistical Software}, \bold{40}, 1--29.
  URL: \doi{10.18637/jss.v040.i01}.
}

\seealso{
CRAN packages \pkg{ROCR}, \pkg{verification} or Bioconductor's \pkg{roc}
for ROC curves.

CRAN packages \pkg{MASS} and \pkg{logcondens} employed in this package.
}
\examples{
data(aSAH)

## Build a ROC object and compute the AUC ##
roc1 <- roc(aSAH$outcome, aSAH$s100b)
print(roc1)

# With a formula
roc(outcome ~ s100b, aSAH)
# With pipes, dplyr-style:
\dontrun{
library(dplyr)
aSAH \%>\% roc(outcome, s100b)}

# Create a few more curves for the next examples
roc2 <- roc(aSAH$outcome, aSAH$wfns)
roc3 <- roc(aSAH$outcome, aSAH$ndka)


## AUC ##
auc(roc1, partial.auc = c(1, .9))


## Smooth ROC curve ##
smooth(roc1)


## Summary statistics
var(roc1)
cov(roc1, roc3)


## Plot the curve ##
plot(roc1)

#  More plotting options, CI and plotting
# with all-in-one syntax:
roc4 <- roc(aSAH$outcome,
            aSAH$s100b, percent=TRUE,
            # arguments for auc
            partial.auc=c(100, 90), partial.auc.correct=TRUE,
            partial.auc.focus="sens",
            # arguments for ci
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            # arguments for plot
            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)

# Add to an existing plot. Beware of 'percent' specification!
roc5 <- roc(aSAH$outcome, aSAH$wfns,
            plot=TRUE, add=TRUE, percent=roc4$percent)


## With ggplot2 ##
if (require(ggplot2)) {
# Create multiple curves to plot
rocs <- roc(outcome ~ wfns + s100b + ndka, data = aSAH)
ggroc(rocs)
}


## Coordinates of the curve ##
coords(roc1, "best", ret=c("threshold", "specificity", "1-npv"))
coords(roc2, "local maximas", ret=c("threshold", "sens", "spec", "ppv", "npv"))


## Confidence intervals ##

# CI of the AUC
ci(roc2)

\dontrun{
# CI of the curve
sens.ci <- ci.se(roc1, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col="lightblue")
plot(sens.ci, type="bars")}

# need to re-add roc2 over the shape
plot(roc2, add=TRUE)

\dontrun{
# CI of thresholds
plot(ci.thresholds(roc2))}

## Comparisons ##

# Test on the whole AUC
roc.test(roc1, roc2, reuse.auc=FALSE)

\dontrun{
# Test on a portion of the whole AUC
roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(100, 90),
         partial.auc.focus="se", partial.auc.correct=TRUE)

# With modified bootstrap parameters
roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(100, 90),
         partial.auc.correct=TRUE, boot.n=1000, boot.stratified=FALSE)}


## Power & sample size ##

# Power
# 1 curve
power.roc.test(roc1)
# 2 curves
power.roc.test(roc3, roc2)

# Sample size 
# 1 curve
power.roc.test(roc3, power = 0.9)
# 2 curves
power.roc.test(roc1, roc2, power = 0.9)

# Also without ROC objects.
# For instance what AUC would be significantly different from 0.5?
power.roc.test(ncases=41, ncontrols=72, sig.level=0.05, power=0.95)
}

\keyword{package}
\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{print}
\keyword{htest}
\keyword{aplot}
\keyword{hplot}
\keyword{roc}
