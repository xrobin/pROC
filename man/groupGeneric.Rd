\encoding{UTF-8}
\name{groupGeneric}
\alias{groupGeneric}
\alias{groupGeneric.pROC}
\alias{groupGeneric.auc}
\alias{groupGeneric.ci.se}
\alias{groupGeneric.ci.sp}
\alias{Ops}
\alias{Math}
\alias{Ops.auc}
\alias{Math.auc}
\alias{Ops.ci.se}
\alias{Math.ci.se}
\alias{Ops.ci.sp}
\alias{Math.ci.sp}
\alias{Ops.ci.auc}
\alias{Math.ci.auc}


\title{
  pROC Group Generic Functions
}
\description{
  Redefine \pkg{base} groupGeneric functions to handle \code{\link{auc}} and \code{\link{ci}} objects (of class \dQuote{ci.auc}, \dQuote{ci.se} and \dQuote{ci.sp}) properly on operations and mathematical operations.
  Attributes are dropped so that the AUC/CI behaves as a numeric value/matrix, respectively.
  In the case of AUC, all attributes are dropped, while in CI only the CI-specific attributes are, keeping those necessary for the matrices.
}
\usage{
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
}
		   
\arguments{
  \item{x, e1, e2}{\code{\link{auc}}, \code{\link{ci.auc}}, \code{\link{ci.se}} or
    \code{\link{ci.sp}} objects, or mixed numerics and such objects.
  }
  \item{\dots}{further arguments passed to other Math methods.
  }
}


\seealso{
\code{\link{groupGeneric}}, \code{\link{auc}}
}
\examples{
data(aSAH)

# Create two auc objects:
auc1 <- auc(roc(aSAH$outcome, aSAH$s100b))
auc2 <- auc(roc(aSAH$outcome, aSAH$wfns))

# Math
sqrt(auc1)
round(auc2, digits = 1)

# Ops
auc1 * 2
2 * auc2
auc1 + auc2

# With CI
ci1 <- ci(auc1)
ci1 * 2
sqrt(ci1)
}

\keyword{methods}
