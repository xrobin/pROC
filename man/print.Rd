\encoding{UTF-8}
\name{print}
\alias{print.roc}
\alias{print.smooth.roc}
\alias{print.multiclass.roc}
\alias{print.mv.multiclass.roc}
\alias{print.ci.auc}
\alias{print.ci.thresholds}
\alias{print.ci.coords}
\alias{print.ci.se}
\alias{print.ci.sp}
\alias{print.auc}
\alias{print.multiclass.auc}

\title{
  Print a ROC curve object
}
\description{
  This function prints a ROC curve, AUC or CI object and return it invisibly.
}
\usage{
\S3method{print}{roc}(x, digits=max(3, getOption("digits") - 3), call=TRUE, ...)
\S3method{print}{multiclass.roc}(x, digits=max(3, getOption("digits") -
  3), call=TRUE, ...) 
\S3method{print}{mv.multiclass.roc}(x, digits=max(3, getOption("digits") -
  3), call=TRUE, ...) 
\S3method{print}{smooth.roc}(x, digits=max(3, getOption("digits") - 3),
call=TRUE, ...)
\S3method{print}{auc}(x, digits=max(3, getOption("digits") - 3), ...)
\S3method{print}{multiclass.auc}(x, digits=max(3, getOption("digits") - 3), ...)
\S3method{print}{ci.auc}(x, digits=max(3, getOption("digits") - 3), ...)
\S3method{print}{ci.thresholds}(x, digits=max(3, getOption("digits") - 3), ...)
\S3method{print}{ci.se}(x, digits=max(3, getOption("digits") - 3), ...)
\S3method{print}{ci.sp}(x, digits=max(3, getOption("digits") - 3), ...)
\S3method{print}{ci.coords}(x, digits=max(3, getOption("digits") - 3), ...)
}

		   
\arguments{
  \item{x}{a roc, auc or ci object, from the \link{roc}, \link{auc} or
    \link{ci} functions respectively.
  }
  \item{call}{if the call is printed.}
  \item{digits}{the number of significant figures to print. See
    \link{signif} for more details.
  }
  \item{\dots}{further arguments passed to or from other methods. In
    particular, \code{print.roc} calls \code{print.auc} and the
    \code{print.ci} variants internally, and a \code{digits} argument is
    propagated. Not used in print.auc and print.ci variants.
  }
}

\value{
  These functions return the object they were passed invisibly.
}

\seealso{
\code{\link{roc}}, \code{\link{auc}}, \code{\link{ci}}, \code{\link{coords}}
}
\examples{
data(aSAH)

# Print a roc object:
rocobj <- roc(aSAH$outcome, aSAH$s100b)
print(rocobj)

# Print a smoothed roc object
print(smooth(rocobj))

# implicit printing
 roc(aSAH$outcome, aSAH$s100b)

# Print an auc and a ci object, from the ROC object or calling
# the dedicated function:
print(rocobj$auc)
print(ci(rocobj))
}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{print}
\keyword{roc}
