\encoding{UTF-8}
\name{coords}
\alias{coords}
\alias{coords.roc}
\alias{coords.smooth.roc}

\title{
  Coordinates of a ROC curve
}
\description{
  This function returns the coordinates of the ROC curve at the
  specified point.
}
\usage{
coords(...)
\S3method{coords}{roc}(roc, x, input="threshold", ret=c("threshold",
"specificity", "sensitivity"),
as.list=FALSE, drop=TRUE, best.method=c("youden", "closest.topleft"),
best.weights=c(1, 0.5), transpose = FALSE, as.matrix=FALSE, ...)
\S3method{coords}{smooth.roc}(smooth.roc, x, input, ret=c("specificity",
"sensitivity"), as.list=FALSE, drop=TRUE, best.method=c("youden",
"closest.topleft"), best.weights=c(1, 0.5), transpose = FALSE,
as.matrix=FALSE, ...)
}
		   
\arguments{
  \item{roc, smooth.roc}{a \dQuote{roc} object from the
	\code{\link{roc}} function, or a \dQuote{smooth.roc} object from the
	\code{\link[=smooth.roc]{smooth}} function.
  }
  \item{x}{
		the coordinates to look for. Numeric (if so, their meaning is
		defined by the \code{input} argument) or one of \dQuote{all} (all
		the points of the ROC curve), \dQuote{local maximas} (the local
		maximas of the ROC curve) or \dQuote{best} (see \code{best.method}
		argument). If missing or \code{NULL}, defaults to \dQuote{all}.
  }
  \item{input}{
  	If \code{x} is numeric, the kind of input coordinate (x).
  	Typically one of \dQuote{threshold}, \dQuote{specificity} or 
  	\dQuote{sensitivity}, but can be any of the monotone coordinate available,
  	see the \dQuote{Valid input} column under \dQuote{Available coordinates}.
  	Can be shortened like \code{ret}. Defaults to \dQuote{threshold}. Note
  	that \dQuote{threshold} is not allowed in \code{coords.smooth.roc} and that
  	the argument is ignored when \code{x} is a character.
  }
  \item{ret}{The coordinates to return. See \dQuote{Available coordinates}
  	section below. Alternatively, the single value \dQuote{all} can be used to return
    every coordinate available.
  }
  \item{as.list}{DEPRECATED. If the returned object must be a list.
   Will be removed in a future version.
  }
  \item{drop}{If \code{TRUE} the result is coerced to the lowest
    possible dimension, as per \link{Extract}. By default only drops
    if \code{transpose = TRUE} and either \code{ret} or \code{x} is
    of length 1.
  }
  \item{best.method}{if \code{x="best"}, the method to determine the
    best threshold. See details in the \sQuote{Best thresholds} section.
  }
  \item{best.weights}{if \code{x="best"}, the weights to determine the
    best threshold. See details in the \sQuote{Best thresholds} section.
  }
  \item{transpose}{whether
  	to return the thresholds in columns (\code{TRUE}) or rows (\code{FALSE}).
  	Since pROC 1.16 the default value is \code{FALSE}.
  	See \link{coords_transpose} for more details the change.
  }
  \item{as.matrix}{if \code{transpose} is \code{FALSE}, whether to return
    a \code{\link{matrix}} (\code{TRUE}) or a \code{\link{data.frame}}
    (\code{FALSE}, the default). A \code{data.frame} is more convenient
    and flexible to use, but incurs a slight speed penalty. Consider
    setting this argument to \code{TRUE} if you are calling the function
    repeatedly.
  }
  \item{\dots}{further arguments passed from other methods. Ignored.}
}

\details{
  This function takes a \dQuote{roc}  or \dQuote{smooth.roc} object as
  first argument, on which the coordinates will be determined. The
  coordinates are defined by the \code{x} and \code{input}
  arguments. \dQuote{threshold} coordinates cannot be determined in a
  smoothed ROC.

  If \code{input="threshold"}, the coordinates for the threshold
  are reported, even if the exact threshold do not define the ROC
  curve. The following convenience characters are allowed: \dQuote{all},
  \dQuote{local maximas} and \dQuote{best}. They will return all the
  thresholds, only the thresholds defining local maximas (upper angles of the
  ROC curve), or only the threshold(s) corresponding to the best sum of
  sensitivity + specificity respectively. Note that \dQuote{best} can
  return more than one threshold. If \code{x} is a character, the
  coordinates are limited to the thresholds within the partial AUC if it
  has been defined, and not necessarily to the whole curve.
  
  For \code{input="specificity"} and \code{input="sensitivity"},
  the function checks if the specificity or sensitivity is one of the
  points of the ROC curve (in \code{roc$sensitivities} or
  \code{roc$specificities}). More than one point may match (in
  \emph{step} curves), then only the upper-left-most point coordinates
  are returned. Otherwise,
  the specificity and specificity of the point is interpolated and
  \code{NA} is returned as threshold.

  The coords function in this package is a generic, but it might be
  superseded by functions in other packages such as
  \pkg{colorspace} or \pkg{spatstat} if they are loaded after
  \pkg{pROC}. In this case, call the \code{pROC::coords} explicitly.

\subsection{Best thresholds}{
  If \code{x="best"}, the \code{best.method} argument controls how the
  optimal threshold is determined.
  \describe{
    \item{\dQuote{youden}}{
      Youden's J statistic (Youden, 1950) is employed. The optimal cut-off is the threshold
      that maximizes the distance to the identity (diagonal) line. Can
      be shortened to \dQuote{y}.

      The optimality criterion is:
      \deqn{max(sensitivities + specificities)}{max(sensitivities + specificities)}
    }
    \item{\dQuote{closest.topleft}}{
      The optimal threshold is the point closest to the top-left part of
      the plot with perfect sensitivity or specificity. Can be shortened
      to \dQuote{c} or \dQuote{t}.

      The optimality criterion is:
      \deqn{min((1 - sensitivities)^2 + (1- specificities)^2)}{min((1 - sensitivities)^2 + (1- specificities)^2)}
    }
  }

  In addition, weights can be supplied if false positive and false
  negative predictions are not equivalent: a numeric vector of length 2
  to the \code{best.weights} argument. The elements define
  \enumerate{
    \item the relative cost of of a false negative classification (as compared with a false positive classification)
    \item the prevalence, or the proportion of cases in the population (\eqn{\frac{n_{cases}}{n_{controls}+n_{cases}}}{n.cases/(n.controls+n.cases)}). 
  }

  The optimality criteria are modified as proposed by Perkins and Schisterman:
  \describe{
    \item{\dQuote{youden}}{
      \deqn{max(sensitivities + r * specificities)}{max(sensitivities + r \times specificities)}
    }
    \item{\dQuote{closest.topleft}}{
      \deqn{min((1 - sensitivities)^2 + r * (1- specificities)^2)}{min((1 - sensitivities)^2 + r \times (1- specificities)^2)}
    }
  }

  with

  \deqn{r = \frac{1 - prevalence}{cost * prevalence}}{r = (1 - prevalence) / (cost * prevalence)}

  By default, prevalence is 0.5 and cost is 1 so that no weight is
  applied in effect.

  Note that several thresholds might be equally optimal.
}

\subsection{Available coordinates}{

	The following table lists the coordinates that are available in the \code{ret}
	and \code{input} arguments.

	\tabular{rllll}{
		Value \tab Description \tab Formula \tab Synonyms \tab Valid input \cr
		\code{threshold} \tab The threshold value \tab - \tab - \tab Yes \cr
		\code{tn} \tab True negative count \tab - \tab - \tab Yes \cr
		\code{tp} \tab True positive count \tab - \tab - \tab Yes \cr
		\code{fn} \tab False negative count \tab - \tab - \tab Yes \cr
		\code{fp} \tab False positive count \tab - \tab - \tab Yes \cr
		\code{specificity} \tab Specificity  \tab tn / (tn + fp) \tab tnr \tab Yes \cr
		\code{sensitivity} \tab Sensitivity \tab tp / (tp + fn)  \tab recall, tpr \tab Yes \cr
		\code{accuracy} \tab Accuracy \tab  (tp + tn) / N \tab - \tab No \cr
		\code{npv} \tab Negative Predictive Value \tab tn / (tn + fn) \tab - \tab No \cr
		\code{ppv} \tab Positive Predictive Value \tab tp / (tp + fp) \tab precision \tab No \cr
		\code{precision} \tab Precision \tab tp / (tp + fp) \tab ppv \tab No \cr
		\code{recall} \tab Recall \tab tp / (tp + fn) \tab sensitivity, tpr \tab Yes \cr
		\code{tpr} \tab True Positive Rate \tab tp / (tp + fn) \tab sensitivity, recall \tab Yes \cr
		\code{fpr} \tab False Positive Rate \tab fp / (tn + fp) \tab 1-specificity \tab Yes \cr
		\code{tnr} \tab True Negative Rate \tab tn / (tn + fp) \tab specificity \tab Yes \cr
		\code{fnr} \tab False Negative Rate \tab fn / (tp + fn) \tab 1-sensitivity \tab Yes \cr
		\code{fdr} \tab False Discovery Rate \tab fp / (tp + fp) \tab 1-ppv \tab No \cr
		\code{youden} \tab Youden Index \tab 
			se + r * sp
			\tab - \tab No \cr
		\code{closest.topleft} \tab Distance to the top left corner of the ROC space \tab 
			- ((1 - se)^2 + r * (1 - sp)^2)
			\tab - \tab No \cr
	}
	
	The value \dQuote{threshold} is not allowed in \code{coords.smooth.roc}.
	
	Values can be shortenend (for example to \dQuote{thr}, \dQuote{sens} and \dQuote{spec}, or even to
	\dQuote{se}, \dQuote{sp} or \dQuote{1-np}). In addition, some values can be prefixed with
	\code{1-} to get their complement:
	\code{1-specificity}, \code{1-sensitivity}, \code{1-accuracy}, \code{1-npv}, \code{1-ppv}.
	
	The values \code{npe} and \code{ppe} are automatically replaced with
	\code{1-npv} and \code{1-ppv}, respectively (and will therefore not appear
	as is in the output, but as \code{1-npv} and \code{1-ppv} instead).
	These must be used verbatim in ROC curves with \code{percent=TRUE} 
	(ie. \dQuote{100-ppv} is never accepted).
	
	The \dQuote{youden} and \dQuote{closest.topleft} are weighted with \code{r}, 
	according to the value of the \code{best.weights} argument. See the 
	\dQuote{Best thresholds} section above for more details.
	
	For \code{ret}, the single value \dQuote{all} can be used to return
	every coordinate available.
}
} % details

\value{
  Depending on the length of \code{x} and \code{as.list} argument.

  
  \tabular{lll}{
	\tab
	length(x) == 1 or length(ret) == 1 \tab
	length(x) > 1 or length(ret) > 1 or drop == FALSE
	\cr
	
    \code{as.list=TRUE} \tab
	a list of the length of, in the order of, and named after, \code{ret}. \tab
	a list of the length of, and named after, \code{x}. Each element of this list is a list of the length of, in the order of, and named after, \code{ret}. \cr

    \code{as.list=FALSE}  \tab
	a numeric vector of the length of, in the order of, and named after, \code{ret} (if \code{length(x) == 1})
	or a numeric vector of the length of, in the order of, and named after, \code{x} (if \code{length(ret) == 1}.\tab
	a numeric matrix with one row for each \code{ret} and one column for each \code{x}\cr
  }
  
  In all cases if \code{input="specificity"} or \code{input="sensitivity"}
  and interpolation was required, threshold is returned as \code{NA}.

  Note that if giving a character as \code{x} (\dQuote{all},
  \dQuote{local maximas} or \dQuote{best}), you cannot predict the
  dimension of the return value unless \code{drop=FALSE}. Even
  \dQuote{best} may return more than one value (for example if the ROC
  curve is below the identity line, both extreme points). 

  \code{coords} may also return \code{NULL} when there a partial area is
  defined but no point of the ROC curve falls within the region.
}

\references{
  Neil J. Perkins, Enrique F. Schisterman (2006) ``The Inconsistency of "Optimal" Cutpoints
  Obtained using Two Criteria based on the Receiver Operating
  Characteristic Curve''. \emph{American Journal of Epidemiology}
  \bold{163}(7), 670--675. DOI: \doi{10.1093/aje/kwj063}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.

  W. J. Youden (1950) ``Index for rating diagnostic tests''. \emph{Cancer}, 
  \bold{3}, 32--35. DOI: 
  \doi{10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3}.
}

\seealso{
\code{\link{roc}}, \code{\link{ci.coords}}
}
\examples{

# Create a ROC curve:
data(aSAH)
roc.s100b <- roc(aSAH$outcome, aSAH$s100b, percent = TRUE)

# Get the coordinates of S100B threshold 0.55
coords(roc.s100b, 0.55, transpose = FALSE)

# Get the coordinates at 50\% sensitivity
coords(roc=roc.s100b, x=50, input="sensitivity", transpose = FALSE)
# Can be abbreviated:
coords(roc.s100b, 50, "se", transpose = FALSE)

# Works with smoothed ROC curves
coords(smooth(roc.s100b), 90, "specificity", transpose = FALSE)

# Get the sensitivities for all thresholds
cc <- coords(roc.s100b, "all", ret="sensitivity", transpose = FALSE)
print(cc$sensitivity)

# Get the best threshold
coords(roc.s100b, "best", ret="threshold", transpose = FALSE)

# Get the best threshold according to different methods
roc.ndka <- roc(aSAH$outcome, aSAH$ndka, percent=TRUE)
coords(roc.ndka, "best", ret="threshold", transpose = FALSE, 
       best.method="youden") # default
coords(roc.ndka, "best", ret="threshold", transpose = FALSE, 
       best.method="closest.topleft")

# and with different weights
coords(roc.ndka, "best", ret="threshold", transpose = FALSE, 
       best.method="youden", best.weights=c(50, 0.2))
coords(roc.ndka, "best", ret="threshold", transpose = FALSE, 
       best.method="closest.topleft", best.weights=c(5, 0.2))
       
# This is available with the plot.roc function too:
plot(roc.ndka, print.thres="best", print.thres.best.method="youden",
                                 print.thres.best.weights=c(50, 0.2)) 

# Return more values:
coords(roc.s100b, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                           "precision", "recall"), transpose = FALSE)

# Return all values
coords(roc.s100b, "best", ret = "all", transpose = FALSE)
                           
# You can use coords to plot for instance a sensitivity + specificity vs. cut-off diagram
plot(specificity + sensitivity ~ threshold, 
     coords(roc.ndka, "all", transpose = FALSE), 
     type = "l", log="x", 
     subset = is.finite(threshold))

# Plot the Precision-Recall curve
plot(precision ~ recall, 
     coords(roc.ndka, "all", ret = c("recall", "precision"), transpose = FALSE),
     type="l", ylim = c(0, 100))

# Alternatively plot the curve with TPR and FPR instead of SE/SP 
# (identical curve, only the axis change)
plot(tpr ~ fpr, 
     coords(roc.ndka, "all", ret = c("tpr", "fpr"), transpose = FALSE),
     type="l")
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
