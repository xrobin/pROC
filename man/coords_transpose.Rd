\encoding{UTF-8}
\name{coords_transpose}
\alias{coords_transpose}

\title{Transposing the output of \code{coords}}

\description{This help page desribes recent and upcoming changes in the 
return values of the \code{\link{coords}} function.}

\section{Background information}{

Until the release of pROC 1.16, the \code{coords} function was returning 
a matrix with thresholds in columns, and the coordinate variables in rows.

\preformatted{
data(aSAH)
rocobj <- roc(aSAH$outcome, aSAH$s100b)
coords(rocobj, c(0.05, 0.2, 0.5))
#                   0.05       0.2       0.5
# threshold   0.05000000 0.2000000 0.5000000
# specificity 0.06944444 0.8055556 0.9722222
# sensitivity 0.97560976 0.6341463 0.2926829

}

This format didn't conform to the grammar of the 
\href{https://www.tidyverse.org/}{tidyverse} which has become prevalent in modern
R language.

In addition, the dropping of dimensions by default makes it difficult to guess
what type of data \code{coords} is going to return.

\preformatted{
coords(rocobj, "best")
#   threshold specificity sensitivity 
#   0.2050000   0.8055556   0.6341463 
# A numeric vector

}

Although it is possible to pass \code{drop = FALSE}, the fact that it is not the
default makes the behaviour unintuitive.

In pROC version 1.16, this was changed and \code{coords} now returns
a \code{\link{data.frame}} with the thresholds in rows and measurement in colums by default.

\preformatted{
 coords(rocobj, c(0.05, 0.2, 0.5), transpose = FALSE)
#      threshold specificity sensitivity
# 0.05      0.05  0.06944444   0.9756098
# 0.2       0.20  0.80555556   0.6341463
# 0.5       0.50  0.97222222   0.2926829
}
}

\section{Changes in 1.15}{
\enumerate{
	\item{Addition of the \code{transpose} argument.}
	\item{Display a warning if \code{transpose} is missing. Pass \code{transpose} explicitly to silence the warning.}
	\item{Deprecation of \code{as.list}.}
}
}

\section{Changes in 1.16}{
\enumerate{
	\item{Change of the default \code{transpose} to \code{TRUE}.}
}
THIS CHANGE IS BACKWARDS INCOMPATIBLE AND IS EXPECTED TO BREAK LEGACY CODE.
}

\section{Changes in 1.17}{
\enumerate{
	\item{Dropped the warning if \code{transpose} is missing.}
}
}

\section{Changes in future versions}{
\enumerate{
	\item{Support for  the \code{as.list} argument might be dropped in the future.
	This is still under consideration.}
	\item{
	The \code{transpose} and \code{drop} arguments might be deprecated in the future, but will remain available for a few additional major versions.
	}
}
}

\section{Related changes in ci.coords}{
In version 1.16, the format of the \code{\link{ci.coords}} return value was changed from a matrix-like object with mixed \code{x} and \code{ret} in rows and 3 columns, into a list-like object which should be easier to use programatically.
}

\section{Recommendations}{
If you are writing a new script calling the \code{coords} function, set \code{transpose = FALSE} to silence the warning and benefit from the latest improvements in pROC and obtain a tidy data.
}

\section{See also}{
\href{https://github.com/xrobin/pROC/issues/54}{The GitHub issue tracking the changes described in this manual page}.
}

