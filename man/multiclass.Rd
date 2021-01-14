\encoding{UTF-8}
\name{multiclass.roc}
\alias{multiclass.roc}
\alias{multiclass.roc.default}
\alias{multiclass.roc.formula}
\title{
 Multi-class AUC
}
\description{
  This function builds builds multiple ROC curve to compute the
  multi-class AUC as defined by Hand and Till.
}
\usage{
multiclass.roc(...)
\S3method{multiclass.roc}{formula}(formula, data, ...)
\S3method{multiclass.roc}{default}(response, predictor,
levels=base::levels(as.factor(response)), 
percent=FALSE, direction = c("auto", "<", ">"), ...)

}

\arguments{
  \item{response}{a factor, numeric or character vector of
    responses (true class), typically encoded with 0 (controls) and 1 (cases), as in
    \code{\link{roc}}.
  }
  \item{predictor}{either a numeric vector, containing the value of each
    observation, as in \code{\link{roc}}, or, a matrix giving the decision value
    (e.g. probability) for each class.
  }
  \item{formula}{a formula of the type \code{response~predictor}.}
  \item{data}{a matrix or data.frame containing the variables in the
    formula. See \code{\link{model.frame}} for more details.}
  \item{levels}{the value of the response for controls and cases
    respectively. In contrast with \code{levels} argument to
    \code{\link{roc}}, all the levels are used and
    \link[=combn]{combined} to compute the multiclass AUC.
  }
  \item{percent}{if the sensitivities, specificities and AUC must be
    given in percent (\code{TRUE}) or in fraction (\code{FALSE}, default).
  }
  \item{direction}{in which direction to make the comparison?
    \dQuote{auto} (default for univariate curves):
    automatically define in which group the
    median is higher and take the direction accordingly. 
    Not available for multivariate curves.
    \dQuote{>} (default for multivariate curves):
    if the predictor values for the control group are
    higher than the values of the case group (controls > t >= cases).
    \dQuote{<}: if the predictor values for the control group are lower
    or equal than the values of the case group (controls < t <= cases).
  }
  \item{...}{further arguments passed to \code{\link{roc}}.
  }
}
\details{
This function performs multiclass AUC as defined by Hand and Till
(2001). A multiclass AUC is a mean of several \code{\link{auc}} and
cannot be plotted. Only AUCs can be computed for such curves.
Confidence intervals, standard deviation, smoothing and
comparison tests are not implemented.

The \code{multiclass.roc} function can handle two types of datasets: uni- and multi-variate.
In the univariate case, a single \code{predictor} vector is passed
and all the combinations of responses are assessed.
I the multivariate case, a \code{\link{matrix}} or \code{\link{data.frame}}
is passed as \code{predictor}. The columns must be named according to the
levels of the \code{response}.

This function has been much less tested than the rest of the package and
is more subject to bugs. Please report them if you find one.
}

\value{
  If \code{predictor} is a vector, a list of class \dQuote{multiclass.roc} 
  (univariate) or \dQuote{mv.multiclass.roc} (multivariate), 
  with the following fields: 
  \item{auc}{if called with \code{auc=TRUE}, a numeric of class \dQuote{auc} as
    defined in \code{\link{auc}}. Note that this is not the standard AUC
    but the multi-class AUC as defined by Hand and Till.
  }
  \item{ci}{if called with \code{ci=TRUE}, a numeric of class \dQuote{ci} as
    defined in \code{\link{ci}}.
  }
  \item{response}{the response vector as passed in argument. If
    \code{NA} values were removed, a \code{na.action} attribute similar
    to \code{\link{na.omit}} stores the row numbers.
  }
  \item{predictor}{the predictor vector as passed in argument. If
    \code{NA} values were removed, a \code{na.action} attribute similar
    to \code{\link{na.omit}} stores the row numbers.
  }
  \item{levels}{the levels of the response as defined in argument.}
  \item{percent}{if the sensitivities, specificities and AUC are
    reported in percent, as defined in argument.
  }
  \item{call}{how the function was called. See \code{\link{match.call}} for
    more details.
  }
}

\section{Warnings}{
  If \code{response} is an ordered factor and one of the levels
  specified in \code{levels} is missing, a warning is issued and the
  level is ignored.
}

\references{
  David J. Hand and Robert J. Till (2001). A Simple Generalisation of
  the Area Under the ROC Curve for Multiple Class Classification
  Problems. \emph{Machine Learning} \bold{45}(2), p. 171--186. DOI:
  \doi{10.1023/A:1010920819831}.
}

\seealso{
 \code{\link{auc}}
}

\examples{
####
# Examples for a univariate decision value
####
data(aSAH)

# Basic example
multiclass.roc(aSAH$gos6, aSAH$s100b)
# Produces an innocuous warning because one level has no observation

# Select only 3 of the aSAH$gos6 levels:
multiclass.roc(aSAH$gos6, aSAH$s100b, levels=c(3, 4, 5))

# Give the result in percent
multiclass.roc(aSAH$gos6, aSAH$s100b, percent=TRUE)

####
# Examples for multivariate decision values (e.g. class probabilities)
####

\dontrun{
# Example with a multinomial log-linear model from nnet
# We use the iris dataset and split into a training and test set
requireNamespace("nnet")
data(iris)
iris.sample <- sample(1:150)
iris.train <- iris[iris.sample[1:75],]
iris.test <- iris[iris.sample[76:150],]
mn.net <- nnet::multinom(Species ~ ., iris.train)

# Use predict with type="prob" to get class probabilities
iris.predictions <- predict(mn.net, newdata=iris.test, type="prob")
head(iris.predictions)

# This can be used directly in multiclass.roc:
multiclass.roc(iris.test$Species, iris.predictions)
}


# Let's see an other example with an artificial dataset
n <- c(100, 80, 150)
responses <- factor(c(rep("X1", n[1]), rep("X2", n[2]), rep("X3", n[3])))
# construct prediction matrix: one column per class

preds <- lapply(n, function(x) runif(x, 0.4, 0.6))
predictor <- as.matrix(data.frame(
                "X1" = c(preds[[1]], runif(n[2] + n[3], 0, 0.7)),
                "X2" = c(runif(n[1], 0.1, 0.4), preds[[2]], runif(n[3], 0.2, 0.8)),
                "X3" = c(runif(n[1] + n[2], 0.3, 0.7), preds[[3]])
             ))
multiclass.roc(responses, predictor)

# One can change direction , partial.auc, percent, etc:
multiclass.roc(responses, predictor, direction = ">")
multiclass.roc(responses, predictor, percent = TRUE, 
	partial.auc = c(100, 90), partial.auc.focus = "se")


# Limit set of levels
multiclass.roc(responses, predictor, levels = c("X1", "X2"))
# Use with formula. Here we need a data.frame to store the responses as characters
data <- cbind(as.data.frame(predictor), "response" = responses)
multiclass.roc(response ~ X1+X3, data)

}

\keyword{univar}
\keyword{nonparametric}
\keyword{utilities}
\keyword{roc}
