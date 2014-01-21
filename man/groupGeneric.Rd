\encoding{UTF-8}
\name{groupGeneric}
\alias{groupGeneric}
\alias{groupGeneric.pROC}
\alias{groupGeneric.auc}
\alias{Ops}
\alias{Math}
\alias{Ops.auc}
\alias{Math.auc}


\title{
  pROC Group Generic Functions
}
\description{
  Redefine \pkg{base} groupGeneric functions to handle \code{\link{auc}} objects properly on operations and mathematical operations.
  Attributes are dropped so that the AUC behaves as a numeric value.
}
\usage{
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
}
		   
\arguments{
  \item{x, e1, e2}{\code{\link{auc}} objects, or mixed numerics and \code{auc} objects.
  }
  \item{\dots}{further arguments passed to other Math methods.
  }
}


\seealso{
\code{\link{groupGeneric}}, \code{\link{auc}}
}
\examples{
data(aSAH)

# Create a roc object:
aucobj1 <- auc(roc(aSAH$outcome, aSAH$s100b))
aucobj2 <- auc(roc(aSAH$outcome, aSAH$wfns))

# Math
sqrt(aucobj1)
round(aucobj2, digits=1)

# Ops
aucobj1 * 2
2 * aucobj2
aucobj1 + aucobj2

}

\keyword{methods}
