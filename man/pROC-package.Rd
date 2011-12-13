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

  In S+, the package comes with a  graphical user interface.
}

\details{
  The basic unit of the pROC package is the \code{\link{roc}} function. It
  will build a ROC curve, smooth it if requested (if \code{smooth=TRUE}),
  compute the AUC (if \code{auc=TRUE}), the confidence interval (CI) if 
  requested (if \code{ci=TRUE}) and plot the curve if requested (if
  \code{plot=TRUE}). 

  The \code{\link{roc}} function will call \code{\link{smooth.roc}},
  \code{\link{auc}},
  \code{\link{ci}} and \code{\link{plot}} as necessary. See these
  individual functions for the arguments that can be passed to them
  through \code{\link{roc}}. These function can be called separately.

  Two paired (that is \code{\link{roc}} objects with the same
  \code{response}) or unpaired (with different \code{response}) ROC
  curves can be compared with the \code{\link{roc.test}} function.
  Sample size and power computations can be performed with the 
  \code{\link{power.roc.test}} function.
}

\section{Citation}{
  If you use pROC in published research, please cite the following paper:

  Xavier Robin, Natacha Turck, Alexandre Hainard, Natalia Tiberti,
  Frédérique Lisacek, Jean-Charles Sanchez and Markus Müller (2011).
  ``pROC: an open-source package for R and S+ to analyze and compare ROC
  curves''. \emph{BMC Bioinformatics}, \bold{12}, p. 77. DOI: 
  \href{http://dx.doi.org/10.1186/1471-2105-12-77}{10.1186/1471-2105-12-77}

  Type \code{citation("pROC")} for a BibTeX entry.

  The authors would be glad to hear how pROC is employed. You are kindly
  encouraged to notify Xavier Robin <Xavier.Robin@unige.ch>
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
    \code{\link{ci.se}} \tab Compute the CI of sensitivities at given specificities \cr
    \code{\link{ci.sp}} \tab Compute the CI of specificities at given sensitivities \cr
    \code{\link{ci.thresholds}} \tab Compute the CI of specificity and sensitivity of thresholds \cr
    \code{\link{coords}} \tab Coordinates of a ROC curve \cr
    \code{\link[=cov.roc]{cov}} \tab Covariance between two AUCs\cr
    \code{\link{has.partial.auc}} \tab Determine if the ROC curve have a partial AUC\cr
    \code{\link{lines.roc}} \tab Add a ROC line to a ROC plot \cr
    \code{\link{plot.ci}} \tab Plot CIs \cr
    \code{\link[=plot.roc]{plot}} \tab Plot a ROC curve \cr
    \code{\link{power.roc.test}} \tab Sample size and power computation \cr
    \code{\link[=print.roc]{print}} \tab Print a ROC curve object \cr
    \code{\link{roc.test}} \tab Compare the AUC of two ROC curves \cr
    \code{\link{smooth}} \tab Smooth a ROC curve \cr
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

\section{Bootstrap}{

  All the bootstrap operations for \link[=roc.test]{significance testing} and \link[=ci]{confidence interval} computation are performed with non-parametric stratified or non-stratified resampling (according to the \code{stratified} argument) and with the percentile method, as described in Carpenter and Bithell (2000) sections 2.1 and 3.3.

  Stratification of bootstrap can be controlled
  with \code{boot.stratified}. In stratified bootstrap (the default), each replicate
  contains the same number of cases and controls than the original
  sample. Stratification is especially useful if one group has only
  little observations, or if groups are not balanced.

  The number of bootstrap replicates is controlled by \code{boot.n}. Higher numbers will give a more precise estimate of the significance tests and confidence intervals
  but take more time to compute. 2000 is recommanded by Carpenter and Bithell (2000) for confidence intervals. In our experience this is sufficient for a good estimation of the 
  first significant digit only, so we recommend the use of 10000 bootstrap replicates to obtain a good estimate of the second significant digit whenever possible.

	\subsection{Progress bars}{
	  A progressbar shows the progress of bootstrap operations. It is handled by the \pkg{plyr} package (Wickham, 2011),
	  and is created by the \code{progress_*} family of functions.
	  Sensible defaults are guessed during the package loading:
	  \itemize{
	    \item In non-\link{interactive} mode, no progressbar is displayed.
	    \item In embedded GNU Emacs \dQuote{ESS}, a \code{\link{txtProgressBar}}
#ifdef windows
	    \item In Windows, a \code{\link{winProgressBar}} bar.
#endif
#ifdef unix 
	    \item In Windows, a \code{winProgressBar} bar.
#endif
	    \item In other systems with a graphical display, a \code{\link{tkProgressBar}}.
	    \item In systems without a graphical display, a \code{\link{txtProgressBar}}.
	  }

	  The default can be changed with the option \dQuote{pROCProgress}. The option must be a list with
	  a \code{name} item setting the type of progress bar (\dQuote{none}, \dQuote{win}, \dQuote{tk}
	  or \dQuote{text}). Optional items of the list are \dQuote{width}, \dQuote{char} and \dQuote{style},
	  corresponding to the arguments to the underlying progressbar functions.
	  For example, to force a text progress bar:
	  
	  \preformatted{options(pROCProgress = list(name = "text", width = NA, char = "=", style = 3)}

	  To inhibit the progress bars completely:
	  
	  \preformatted{options(pROCProgress = list(name = "none"))}
	}
}

\author{
Xavier Robin, Natacha Turck, Jean-Charles Sanchez and Markus Müller

Maintainer: Xavier Robin <Xavier.Robin@unige.ch>
}

\references{
  James Carpenter and John Bithell (2000) ``Bootstrap condence intervals:
  when, which, what? A practical guide for medical statisticians''.
  \emph{Statistics in Medicine} \bold{19}, 1141--1164.
  DOI: \href{http://dx.doi.org/10.1002/(SICI)1097-0258(20000515)19:9<1141::AID-SIM479>3.0.CO;2-F}{10.1002/(SICI)1097-0258(20000515)19:9<1141::AID-SIM479>3.0.CO;2-F}.

  Tom Fawcett (2006) ``An introduction to ROC analysis''. \emph{Pattern
    Recognition Letters} \bold{27}, 861--874. DOI:
  \href{http://dx.doi.org/10.1016/j.patrec.2005.10.010}{10.1016/j.patrec.2005.10.010}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \href{http://dx.doi.org/10.1186/1471-2105-12-77}{10.1186/1471-2105-12-77}.
  
  Hadley Wickham (2011) ``The Split-Apply-Combine Strategy for Data Analysis''. \emph{Journal of Statistical Software}, \bold{40}, 1--29.
  URL: \href{http://www.jstatsoft.org/v40/i01}{www.jstatsoft.org/v40/i01}.
}

\seealso{
CRAN packages \pkg{ROCR}, \pkg{verification} or Bioconductor's \pkg{roc}
for ROC curves.

CRAN packages \pkg{plyr}, \pkg{MASS} and \pkg{logcondens} employed in this package.
}
\examples{
data(aSAH)

# Build a ROC object and compute the AUC
roc(aSAH$outcome, aSAH$s100b)
roc(outcome ~ s100b, aSAH)

# Smooth ROC curve
roc(outcome ~ s100b, aSAH, smooth=TRUE)

# more options, CI and plotting
roc1 <- roc(aSAH$outcome,
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
roc2 <- roc(aSAH$outcome, aSAH$wfns,
            plot=TRUE, add=TRUE, percent=roc1$percent)

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
plot(sens.ci, type="bars")
}

# need to re-add roc2 over the shape
plot(roc2, add=TRUE)

\dontrun{
# CI of thresholds
plot(ci.thresholds(roc2))
}


## Comparisons ##

# Test on the whole AUC
roc.test(roc1, roc2, reuse.auc=FALSE)

\dontrun{
# Test on a portion of the whole AUC
roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(100, 90),
         partial.auc.focus="se", partial.auc.correct=TRUE)

# With modified bootstrap parameters
roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(100, 90),
         partial.auc.correct=TRUE, boot.n=1000, boot.stratified=FALSE)
}


## Sample size ##

# Two ROC curves
power.roc.test(roc1, roc2, reuse.auc=FALSE)
power.roc.test(roc1, roc2, power=0.9, reuse.auc=FALSE)

# One ROC curve
power.roc.test(auc=0.8, ncases=41, ncontrols=72)
power.roc.test(auc=0.8, power=0.9)
power.roc.test(auc=0.8, ncases=41, ncontrols=72, sig.level=0.01)
power.roc.test(ncases=41, ncontrols=72, power=0.9)
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
