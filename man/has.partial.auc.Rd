\encoding{UTF-8}
\name{has.partial.auc}
\alias{has.partial.auc}
\alias{has.partial.auc.smooth.roc}
\alias{has.partial.auc.auc}
\alias{has.partial.auc.roc}

\title{
  Does the ROC curve have a partial AUC?
}
\description{
  This function determines if the ROC curve has a partial AUC.
}
\usage{
has.partial.auc(roc)
\S3method{has.partial.auc}{auc}(roc)
\S3method{has.partial.auc}{smooth.roc}(roc)
\S3method{has.partial.auc}{roc}(roc)
}
		   
\arguments{
  \item{roc}{the ROC curve to check.}
}

\value{
  \code{TRUE} if the AUC is a partial AUC, \code{FALSE} otherwise.

  If the AUC is not defined (i. e. if roc was called with \code{AUC=FALSE}), returns \code{NULL}.
}

\seealso{
\code{\link{auc}}
}

\examples{
data(aSAH)

# Full AUC
roc1 <- roc(aSAH$outcome, aSAH$s100b)
has.partial.auc(roc1)
has.partial.auc(auc(roc1))
has.partial.auc(smooth(roc1))

# Partial AUC
roc2 <- roc(aSAH$outcome, aSAH$s100b, partial.auc = c(1, 0.9))
has.partial.auc(roc2)
has.partial.auc(smooth(roc2))

# No AUC
roc3 <- roc(aSAH$outcome, aSAH$s100b, auc = FALSE)
has.partial.auc(roc3)
}

\keyword{programming}
\keyword{logic}
\keyword{roc}
