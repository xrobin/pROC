\encoding{UTF-8}
\name{aSAH}
\docType{data}
\alias{aSAH}
\title{
 Subarachnoid hemorrhage data
}
\description{
  This dataset summarizes several clinical and one laboratory variable
  of 113 patients with an aneurysmal subarachnoid hemorrhage.
}
\usage{aSAH}

\format{A data.frame containing 113 observations of 7 variables.}

\source{
  Natacha Turck, Laszlo Vutskits, Paola Sanchez-Pena, Xavier
  Robin, Alexandre Hainard, Marianne Gex-Fabry, Catherine Fouda, Hadiji
  Bassem, Markus Mueller, Frédérique Lisacek, Louis Puybasset and
  Jean-Charles Sanchez (2010) ``A multiparameter panel method for outcome
  prediction following aneurysmal subarachnoid hemorrhage''.
  \emph{Intensive Care Medicine} \bold{36}(1), 107--115. DOI:
  \doi{10.1007/s00134-009-1641-y}.
}

\seealso{
  Other examples can be found in all the documentation pages of this
  package: 
  \code{\link{roc}}, \code{\link{auc}}, \code{\link{ci}},
  \code{\link{ci.auc}}, \code{\link{ci.se}}, \code{\link{ci.sp}},
  \code{\link{ci.thresholds}}, \code{\link{coords}},
  \code{\link{plot.ci}}, \code{\link{plot.roc}},
  \code{\link{print.roc}}, \code{\link{roc.test}} and
  \code{\link{smooth}}.

  An example analysis with pROC is shown in:
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}
}

\examples{
# load the dataset
data(aSAH)

# Gender, outcome and set
with(aSAH, table(gender, outcome))

# Age
with(aSAH, by(age, outcome, mean))
with(aSAH, by(age, outcome,
     function(x) sprintf("mean: \%.1f (+/- \%.1f), median: \%.1f (\%i-\%i)",
                         mean(x), sd(x), median(x), min(x), max(x))))

# WFNS score
with(aSAH, table(wfns=ifelse(wfns<=2, "1-2", "3-4-5"), outcome))

}

\keyword{datasets}
