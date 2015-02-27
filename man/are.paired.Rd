\encoding{UTF-8}
\name{are.paired}
\alias{are.paired}
\alias{are.paired.smooth.roc}
\alias{are.paired.auc}
\alias{are.paired.roc}

\title{
  Are two ROC curves paired?
}
\description{
  This function determines if two ROC curves can be paired.
}
\usage{
are.paired(...)
\S3method{are.paired}{auc}(roc1, roc2, ...)
\S3method{are.paired}{smooth.roc}(roc1, roc2, ...)
\S3method{are.paired}{roc}(roc1, roc2, return.paired.rocs=FALSE,
  reuse.auc = TRUE, reuse.ci = FALSE, reuse.smooth=TRUE, ...)
}
		   
\arguments{
  \item{roc1, roc2}{the two ROC curves to compare. Either
    \dQuote{\link{roc}}, \dQuote{\link{auc}} or
    \dQuote{\link{smooth.roc}} objects (types can be mixed).
  }
  \item{return.paired.rocs}{if \code{TRUE} and the ROC curves can be
    paired, the two paired ROC curves with \code{NA}s removed will be
    returned.
  }
  \item{reuse.auc, reuse.ci, reuse.smooth}{
    if \code{return.paired.rocs=TRUE}, determines if \code{\link{auc}},
    \code{\link{ci}} and \code{\link[=smooth.roc]{smooth}} should be re-computed
    (with the same parameters than the original ROC curves)
  }
  \item{\dots}{additionnal arguments for \code{are.paired.roc}. Ignored
    in \code{are.paired.roc}
  }
}

\details{
  Two ROC curves are paired if they are built on two variables observed
  on the same sample.

  In practice, the paired status is granted if the \code{response} and \code{levels} vector
  of both ROC curves are \link{identical}. If the \code{response}s are different, this can be
  due to missing values differing between the curves. In this case, the
  function will strip all \code{NA}s in both curves and check for
  identity again.

  It can raise false positives if the responses are identical but correspond
  to different patients.
}

\value{
  \code{TRUE} if \code{roc1} and \code{roc2} are paired, \code{FALSE}
  otherwise. 

  In addition, if \code{TRUE} and \code{return.paired.rocs=TRUE}, the
  following atributes are defined:
  \item{roc1, roc2}{the two ROC curve with all \code{NA}s values removed
    in both curves.
  }
}

\seealso{
\code{\link{roc}}, \code{\link{roc.test}}
}

\examples{
data(aSAH)
aSAH.copy <- aSAH

# artificially insert NAs for demonstration purposes
aSAH.copy$outcome[42] <- NA
aSAH.copy$s100b[24] <- NA
aSAH.copy$ndka[1:10] <- NA

# Call roc() on the whole data
roc1 <- roc(aSAH.copy$outcome, aSAH.copy$s100b)
roc2 <- roc(aSAH.copy$outcome, aSAH.copy$ndka)
# are.paired can still find that the curves were paired
are.paired(roc1, roc2) # TRUE

# Removing the NAs manually before passing to roc() un-pairs the ROC curves
nas <- is.na(aSAH.copy$outcome) | is.na(aSAH.copy$ndka)
roc2b <- roc(aSAH.copy$outcome[!nas], aSAH.copy$ndka[!nas])
are.paired(roc1, roc2b) # FALSE

# Getting the two paired ROC curves with additional smoothing and ci options
roc2$ci <- ci(roc2)
paired <- are.paired(smooth(roc1), roc2, return.paired.rocs=TRUE, reuse.ci=TRUE)
paired.roc1 <- attr(paired, "roc1")
paired.roc2 <- attr(paired, "roc2")

}

\keyword{programming}
\keyword{logic}
\keyword{roc}
