\encoding{UTF-8}
\name{roc}
\alias{roc}
\alias{roc.formula}
\alias{roc.default}
\title{
 Build a ROC curve
}
\description{
  This is the main function of the pROC package. It builds a ROC
  curve and returns a \dQuote{roc} object, a list of class
  \dQuote{roc}. This object can be \code{prin}ted, \code{plot}ted, or
  passed to the functions \code{\link{auc}}, \code{\link{ci}},
  \code{\link{smooth.roc}} and \code{\link{coords}}. Additionally, two
  \code{roc} objects can be compared with \code{\link{roc.test}}. 
}
\usage{
roc(...)
\S3method{roc}{formula}(formula, ...)
\S3method{roc}{default}(response, predictor, controls, cases,
density.controls, density.cases,
levels=base::levels(as.factor(response)), percent=FALSE, na.rm=TRUE,
direction=c("auto", "<", ">"), algorithm = 1, smooth=FALSE, auc=TRUE, ci=FALSE,
plot=FALSE, smooth.method="binormal", ci.method=NULL, density=NULL, ...)

}

\arguments{
  \item{response}{a factor, numeric or character vector of
    responses, typically encoded with 0 (controls) and 1 (cases).
    Only two classes can be used in a ROC curve. If the vector
    contains more than two unique values, or if their order could be
    ambiguous, use \code{levels} to specify which values must be used as
    control and case value.
  }
  \item{predictor}{a numeric vector of the same length than \code{response},
    containing the predicted value of each
    observation. An ordered factor is coerced to a numeric.
  }
  \item{controls, cases}{instead of \code{response}, \code{predictor},
    the data can be supplied as two numeric or ordered vectors containing the predictor
    values for control and case observations.
  }
  \item{density.controls, density.cases}{a smoothed ROC curve can be
    built directly from two densities on identical \code{x} points, as in
    \code{\link[=smooth.roc]{smooth}}.
  }
  \item{formula}{a formula of the type \code{response~predictor}. If mulitple predictors 
    are passed, a named list of \code{roc} objects will be returned. Additional arguments
    \code{data}, \code{subset} and \code{na.action} are supported, see 
    \code{\link{model.frame}} for more details.
  }
  \item{levels}{the value of the response for controls and cases
    respectively. By default, the first two values of
    \code{levels(as.factor(response))} are taken, and the remaining levels are ignored.
    It usually captures two-class factor data correctly, but will
    frequently fail for other data types (response factor with more than 2 levels,
    or for example if your response is coded \dQuote{controls} and \dQuote{cases},
    the levels will be inverted) and must then be specified here.
    If your data is coded as \code{0} and \code{1} with \code{0}
    being the controls, you can safely omit this argument.
  }
  \item{percent}{if the sensitivities, specificities and AUC must be
    given in percent (\code{TRUE}) or in fraction (\code{FALSE}, default).
  }
  \item{na.rm}{if \code{TRUE}, the \code{NA} values will be removed
    (ignored by \code{roc.formula}).
  }
  \item{direction}{in which direction to make the comparison?
    \dQuote{auto} (default): automatically define in which group the
    median is higher and take the direction accordingly.
    \dQuote{>}: if the predictor values for the control group are
    higher than the values of the case group (controls > t >= cases).
    \dQuote{<}: if the predictor values for the control group are lower
    or equal than the values of the case group (controls < t <= cases).
  }
  \item{algorithm}{the method used to compute sensitivity and specificity, 
    an integer of length 1 between \code{0} and \code{4}.
    \code{1} (default): a safe, well-tested, pure-\R code that is efficient when the 
    number of thresholds is low. It goes with O(T*N).
    \code{2}: an alternative pure-\R algorithm that goes in
    O(N). Typically faster than \code{1} when the number of thresholds of
    the ROC curve is above 1000. Less tested than \code{1}.
    \code{3}: a C++ 
    implementation of \code{1}, about 3-5x faster. Typically the fastest with
    ROC curves with less than 3000-5000 thresholds. 
    \code{4} (debug only, slow):
    runs all 3 algorithms and makes sure they return the same values.
    \code{0}: use \pkg{microbenchmark} to choose between \code{2} and \code{3}.
  }
  \item{smooth}{if TRUE, the ROC curve is passed to \code{\link{smooth}}
    to be smoothed.
  }
  \item{auc}{compute the area under the curve (AUC)? If \code{TRUE}
    (default), additional arguments can be passed to \code{\link{auc}}.
  }
  \item{ci}{compute the confidence interval (CI)? If \code{TRUE}
    (default), additional arguments can be passed to \code{\link{ci}}.
  }
  \item{plot}{plot the ROC curve? If \code{TRUE}, additional
    arguments can be passed to \code{\link{plot.roc}}.
  }
  \item{smooth.method, ci.method}{in \code{roc.formula} and
	\code{roc.default}, the \code{method} arguments to
	\code{\link[=smooth.roc]{smooth}} (if \code{smooth=TRUE}) and
	\code{of="auc"}) must be passed as
	\code{smooth.method} and \code{ci.method} to avoid confusions.
  }
  \item{density}{\code{density} argument passed to \code{\link[=smooth.roc]{smooth}}.}
  \item{\dots}{further arguments passed to or from other methods, and
    especially:
    \itemize{
      \item \code{\link{auc}}: \code{partial.auc}, \code{partial.auc.focus}, \code{partial.auc.correct}.
      \item \code{\link{ci}}: \code{of}, \code{conf.level}, \code{boot.n}, \code{boot.stratified}, \code{progress}
      \item \code{\link{ci.auc}}:, \code{reuse.auc}, \code{method}
      \item \code{\link{ci.thresholds}}: \code{thresholds}
      \item \code{\link{ci.sp}}: \code{sensitivities}
      \item \code{\link{ci.se}}: \code{specificities}
      \item \code{\link{plot.roc}}: \code{add}, \code{col} and most
        other arguments to the \code{\link{plot.roc}} function. See
	\code{\link{plot.roc}} directly for more details.
      \item \code{\link{smooth}}: \code{method}, \code{n}, and all other
        arguments. See \code{\link{smooth}} for more details.
    }
  }
}
\details{
  This function's main job is to build a ROC object. See the
  \dQuote{Value} section to this page for more details. Before
  returning, it will call (in this order) the \code{\link[=smooth.roc]{smooth}},
  \code{\link{auc}}, \code{\link{ci}} and \code{\link{plot.roc}}
  functions if \code{smooth} \code{auc}, \code{ci} and \code{plot.roc}
  (respectively) arguments are set to TRUE. By default, only \code{auc}
  is called.

  Data can be provided as \code{response, predictor}, where the
  predictor is the numeric (or ordered) level of the evaluated signal, and
  the response encodes the observation class (control or case). The
  \code{level} argument specifies which response level must be taken as
  controls (first value of \code{level}) or cases (second). It can
  safely be ignored when the response is encoded as \code{0} and
  \code{1}, but it will frequently fail otherwise. By default, the first
  two values of \code{levels(as.factor(response))} are taken, and the
  remaining levels are ignored. This means that if your response is
  coded \dQuote{control} and \dQuote{case}, the levels will be
  inverted.
	
  In some cases, it is more convenient to pass the data as
  \code{controls, cases}, but both arguments are ignored if
  \code{response, predictor} was specified to non-\code{NULL} values.
  It is also possible to pass density data with \code{density.controls,
  density.cases}, which will result in a smoothed ROC curve even if
  \code{smooth=FALSE}, but are ignored if \code{response, predictor} or
  \code{controls, cases} are provided.

  Specifications for \code{\link{auc}}, \code{\link{ci}} and
  \code{\link{plot.roc}} are not kept if \code{auc}, \code{ci} or \code{plot} are set to
  \code{FALSE}. Especially, in the following case:
  
  \preformatted{
    myRoc <- roc(..., auc.polygon=TRUE, grid=TRUE, plot=FALSE)
    plot(myRoc)
  }

  the plot will not have the AUC polygon nor the grid. Similarly, when
  comparing \dQuote{roc} objects, the following is not possible:

  \preformatted{
    roc1 <- roc(..., partial.auc=c(1, 0.8), auc=FALSE)
    roc2 <- roc(..., partial.auc=c(1, 0.8), auc=FALSE)
    roc.test(roc1, roc2)
  }

  This will produce a test on the full AUC, not the partial AUC. To make
  a comparison on the partial AUC, you must repeat the specifications
  when calling \code{\link{roc.test}}:

  \preformatted{
    roc.test(roc1, roc2, partial.auc=c(1, 0.8))
  }

  Note that if \code{roc} was called with \code{auc=TRUE}, the latter syntax will not
  allow redefining the AUC specifications. You must use \code{reuse.auc=FALSE} for that.
  
}
\value{
  If the data contained any \code{NA} value and \code{na.rm=FALSE}, \code{NA} is
  returned. Otherwise, if \code{smooth=FALSE}, a list of class
  \dQuote{roc} with the following fields: 
  \item{auc}{if called with \code{auc=TRUE}, a numeric of class \dQuote{auc} as
    defined in \code{\link{auc}}.
  }
  \item{ci}{if called with \code{ci=TRUE}, a numeric of class \dQuote{ci} as
    defined in \code{\link{ci}}.
  }
  \item{response}{the response vector. Patients whose response is not
    \code{\link{\%in\%}} \code{levels} are discarded. If \code{NA} values
    were removed, a \code{na.action} attribute similar
    to \code{\link{na.omit}} stores the row numbers. 
  }
  \item{predictor}{the predictor vector converted to numeric as used to build the ROC
    curve. Patients whose response is not \code{\link{\%in\%}} \code{levels} are discarded. If
    \code{NA} values were removed, a \code{na.action} attribute similar
    to \code{\link{na.omit}} stores the row numbers.
  }
  \item{original.predictor, original.response}{the response and predictor vectors as passed in argument.}
  \item{levels}{the levels of the response as defined in argument.}
  \item{controls}{the predictor values for the control observations.}
  \item{cases}{the predictor values for the cases.}
  \item{percent}{if the sensitivities, specificities and AUC are
    reported in percent, as defined in argument.
  }
  \item{direction}{the direction of the comparison, as defined in argument.}
  \item{fun.sesp}{the function used to compute sensitivities and specificities.
    Will be re-used in bootstrap operations.}
  \item{sensitivities}{the sensitivities defining the ROC curve.}
  \item{specificities}{the specificities defining the ROC curve.}
  \item{thresholds}{the thresholds at which the sensitivities and
    specificities were computed.
  }
  \item{call}{how the function was called. See \code{\link{match.call}} for
    more details.
  }

  If \code{smooth=TRUE} a list of class \dQuote{smooth.roc} as returned
  by \code{\link{smooth}}, with or without additional elements
  \code{auc} and \code{ci} (according to the call).
}

\section{Errors}{
  If no control or case observation exist for the given levels of
  response, no ROC curve can be built and an error is triggered with
  message \dQuote{No control observation} or \dQuote{No case
    observation}.

  If the predictor is not a numeric or ordered, as defined by
  \code{\link{as.numeric}} or \code{\link{as.ordered}}, the message
  \dQuote{Predictor must be numeric or ordered} is returned.

  The message \dQuote{No valid data provided} is issued when the data
  wasn't properly passed. Remember you need both \code{response} and
  \code{predictor} of the same (not null) length, or bot \code{controls}
  and \code{cases}. Combinations such as \code{predictor} and
  \code{cases} are not valid and will trigger this error.
}

\references{
  Tom Fawcett (2006) ``An introduction to ROC analysis''. \emph{Pattern
    Recognition Letters} \bold{27}, 861--874. DOI:
  \href{http://dx.doi.org/10.1016/j.patrec.2005.10.010}{10.1016/j.patrec.2005.10.010}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \href{http://dx.doi.org/10.1186/1471-2105-12-77}{10.1186/1471-2105-12-77}.
}

\seealso{
 \code{\link{auc}}, \code{\link{ci}}, \code{\link{plot.roc}}, \code{\link{print.roc}}, \code{\link{roc.test}}
}

\examples{
data(aSAH)

# Basic example
roc(aSAH$outcome, aSAH$s100b,
    levels=c("Good", "Poor"))
# As levels aSAH$outcome == c("Good", "Poor"),
# this is equivalent to:
roc(aSAH$outcome, aSAH$s100b)
# In some cases, ignoring levels could lead to unexpected results
# Equivalent syntaxes:
roc(outcome ~ s100b, aSAH)
roc(aSAH$outcome ~ aSAH$s100b)
with(aSAH, roc(outcome, s100b))
with(aSAH, roc(outcome ~ s100b))

# With a formula:
roc(outcome ~ s100b, data=aSAH)

# With numeric controls/cases
roc(controls=aSAH$s100b[aSAH$outcome=="Good"], cases=aSAH$s100b[aSAH$outcome=="Poor"])
# With ordered controls/cases
roc(controls=aSAH$wfns[aSAH$outcome=="Good"], cases=aSAH$wfns[aSAH$outcome=="Poor"])

# Inverted the levels: "Poor" are now controls and "Good" cases:
roc(aSAH$outcome, aSAH$s100b,
    levels=c("Poor", "Good"))

# The result was exactly the same because of direction="auto".
# The following will give an AUC < 0.5:
roc(aSAH$outcome, aSAH$s100b,
    levels=c("Poor", "Good"), direction="<")

# If we prefer counting in percent:
roc(aSAH$outcome, aSAH$s100b, percent=TRUE)

# Test the different algorithms:
roc(aSAH$outcome, aSAH$s100b, algorithm = 1)
roc(aSAH$outcome, aSAH$s100b, algorithm = 2)
roc(aSAH$outcome, aSAH$s100b, algorithm = 3)
if (require(microbenchmark)) {
  roc(aSAH$outcome, aSAH$s100b, algorithm = 0)
}

# Plot and CI (see plot.roc and ci for more options):
roc(aSAH$outcome, aSAH$s100b,
    percent=TRUE, plot=TRUE, ci=TRUE)

# Smoothed ROC curve
roc(aSAH$outcome, aSAH$s100b, smooth=TRUE)
# this is not identical to
smooth(roc(aSAH$outcome, aSAH$s100b))
# because in the latter case, the returned object contains no AUC
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
