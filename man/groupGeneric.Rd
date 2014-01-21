\encoding{UTF-8}
\name{groupGeneric}
\alias{groupGeneric}
\alias{groupGeneric.pROC}
\alias{groupGeneric.auc}
\alias{groupGeneric.ci.coords}
\alias{groupGeneric.ci.se}
\alias{groupGeneric.ci.sp}
\alias{Ops}
\alias{Math}
\alias{Ops.auc}
\alias{Math.auc}
\alias{Ops.ci.coords}
\alias{Math.ci.coords}
\alias{Ops.ci.se}
\alias{Math.ci.se}
\alias{Ops.ci.sp}
\alias{Math.ci.sp}
\alias{Ops.ci.auc}
\alias{Math.ci.auc}
\alias{Ops.ci}
\alias{Math.ci}


\title{
  pROC Group Generic Functions
}
\description{
  Redefine \pkg{base} groupGeneric functions to handle \code{\link{auc}} and \code{\link{ci}} objects properly on operations and mathematical operations.
  Attributes are dropped so that the AUC/CI behaves as a numeric value/matrix, respectively.
  In the case of AUC, all attributes are dropped, while in CI only the CI-specific attributes are, keeping those necessary for the matrices.
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

# With CI
ciaucobj <- ci(aucobj1)
ciaucobj * 2
sqrt(ciaucobj)

}

\keyword{methods}
