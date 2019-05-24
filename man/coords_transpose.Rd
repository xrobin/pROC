\encoding{UTF-8}
\name{coords_transpose}
\alias{coords_transpose}

\title{Transposing the output of \code{coords}}

\description{This help page desribes recent and upcoming changes in the 
return values of the \code{\link{coords}} function.}

\section{Background information}{

Since the initial release of pROC, the \code{coords} has been returning 
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

This format doesn't conform to the grammar of the 
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

In an upcoming version of pROC, this will be changed and \code{coords} will return
a \code{\link{data.frame}} with the thresholds in rows and measurement in colums by default.
}

\section{Changes in 1.15}{
\enumerate{
	\item{Addition of \code{transpose} argument.}
	\item{Display a warning if \code{transpose} is missing. Pass \code{transpose} explicitly to silence the warning.}
	\item{Deprecation of \code{as.list}.}
}

With \code{transpose = FALSE}, the output is a \code{\link{data.frame}} looking like this:
\preformatted{
 coords(rocobj, c(0.05, 0.2, 0.5), transpose = FALSE)
#      threshold specificity sensitivity
# 0.05      0.05  0.06944444   0.9756098
# 0.2       0.20  0.80555556   0.6341463
# 0.5       0.50  0.97222222   0.2926829
}

It is recommended that new developments set \code{transpose = FALSE} explicitly.

These changes are neutral to the API and do not affect functionality outside of a warning. 
}

\section{Upcoming changes in future versions}{
The next version of pROC will change the default \code{transpose} to \code{FALSE}. 
This will break any script that did not previously set \code{transpose} and will
initially come with a warning to make debugging easier. The warning will eventually be removed. Scripts that set \code{transpose} explicitly will be unaffected.
}

\section{Recommendations}{
If you are writing a script calling the \code{coords} function, set \code{transpose = FALSE} to silence the warning and make sure your script keeps running smoothly once the 
default \code{transpose} is changed to \code{FALSE}.
}

\section{See also}{
\href{https://github.com/xrobin/pROC/issues/54}{The GitHub issue tracking the changes described in this manual page}.
}

