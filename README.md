[![Build Status](https://travis-ci.org/xrobin/pROC.svg?branch=master)](https://travis-ci.org/xrobin/pROC)
[![Codecov coverage](https://codecov.io/github/xrobin/pROC/branch/master/graphs/badge.svg)](https://codecov.io/github/xrobin/pROC) 
[![CRAN Version](http://www.r-pkg.org/badges/version/pROC)](http://cran.r-project.org/web/packages/pROC)
[![Downloads](http://cranlogs.r-pkg.org/badges/pROC)](http://cran.rstudio.com/package=pROC)

pROC
=============

An [R](http://r-project.org/) package to display and analyze ROC curves.

For more information, see:

1. Xavier Robin, Natacha Turck, Alexandre Hainard, *et al.* (2011) “pROC: an open-source package for R and S+ to analyze and compare ROC curves”. *BMC Bioinformatics*, **7**, 77. DOI: [10.1186/1471-2105-12-77](http://dx.doi.org/10.1186/1471-2105-12-77)
2. [The official web page on ExPaSy](http://www.expasy.org/tools/pROC/)
3. [The CRAN page](https://cran.r-project.org/package=pROC)
4. [My blog](http://xavier.robin.name/tag/pROC/)
5. [The FAQ](https://github.com/xrobin/pROC/wiki/FAQ---Frequently-asked-questions)

Stable
-------

The latest stable version is best installed from the CRAN:

    install.packages("pROC")

Help
-------

Once the library is loaded with `library(pROC)`, you can get help on pROC by typing `?pROC`.

Getting started
-------

If you don't want to read the manual first, try the following:

### Loading 

```R
library(pROC)
data(aSAH)
```
### Basic ROC / AUC analysis 
```R
roc(aSAH$outcome, aSAH$s100b)
roc(outcome ~ s100b, aSAH)
```
### Smoothing
```R
roc(outcome ~ s100b, aSAH, smooth=TRUE) 
```
### more options, CI and plotting
```R
roc1 <- roc(aSAH$outcome,
            aSAH$s100b, percent=TRUE,
            # arguments for auc
            partial.auc=c(100, 90), partial.auc.correct=TRUE,
            partial.auc.focus="sens",
            # arguments for ci
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            # arguments for plot
            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)

    # Add to an existing plot. Beware of 'percent' specification!
    roc2 <- roc(aSAH$outcome, aSAH$wfns,
            plot=TRUE, add=TRUE, percent=roc1$percent)        
```
### Coordinates of the curve
```R
coords(roc1, "best", ret=c("threshold", "specificity", "1-npv"))
coords(roc2, "local maximas", ret=c("threshold", "sens", "spec", "ppv", "npv"))
```
### Confidence intervals
```R
# Of the AUC
ci(roc2)

# Of the curve
sens.ci <- ci.se(roc1, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col="lightblue")
plot(sens.ci, type="bars")

# need to re-add roc2 over the shape
plot(roc2, add=TRUE)

# CI of thresholds
plot(ci.thresholds(roc2))
```
### Comparisons
```R
    # Test on the whole AUC
    roc.test(roc1, roc2, reuse.auc=FALSE)

    # Test on a portion of the whole AUC
    roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(100, 90),
             partial.auc.focus="se", partial.auc.correct=TRUE)

    # With modified bootstrap parameters
    roc.test(roc1, roc2, reuse.auc=FALSE, partial.auc=c(100, 90),
             partial.auc.correct=TRUE, boot.n=1000, boot.stratified=FALSE)
```
### Sample size
```R
    # Two ROC curves
    power.roc.test(roc1, roc2, reuse.auc=FALSE)
    power.roc.test(roc1, roc2, power=0.9, reuse.auc=FALSE)

    # One ROC curve
    power.roc.test(auc=0.8, ncases=41, ncontrols=72)
    power.roc.test(auc=0.8, power=0.9)
    power.roc.test(auc=0.8, ncases=41, ncontrols=72, sig.level=0.01)
    power.roc.test(ncases=41, ncontrols=72, power=0.9)
```


Development
-------

Download the source code from git, unzip it if necessary, and then type `R CMD INSTALL pROC`. Alternatively, you can use the [devtools](https://github.com/hadley/devtools/wiki) package by [Hadley Wickham](http://had.co.nz/) to automate the process (make sure you follow [the full instructions to get started](http://www.rstudio.com/projects/devtools/)):

```R
install.packages("devtools")
devtools::install_github("xrobin/pROC")
```
