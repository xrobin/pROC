
- [Validation of the maximum-likelihood binormal
  fit](#validation-of-the-maximum-likelihood-binormal-fit)
  - [Design](#design)
  - [Small datasets (D&B Tables 1, 3, 4 and
    5)](#small-datasets-db-tables-1-3-4-and-5)
  - [Large datasets (D&B Table 2)](#large-datasets-db-table-2)
  - [Standard error of the AUC (delta
    method)](#standard-error-of-the-auc-delta-method)
  - [Discussion](#discussion)

<!-- binormal_validation.md is generated from binormal_validation.Rmd. Please edit that file -->

# Validation of the maximum-likelihood binormal fit

This report checks `smooth(method = "binormal")`, the maximum-likelihood
fit of the binormal model $\Phi^{-1}(TPF) = a + b \Phi^{-1}(FPF)$,
against the published results of ROCFIT, the maximum-likelihood program
of Metz, and the standard error of its AUC by the delta method
(`var(method = "ml")`). It is maintainer-only: knit it from the package
root with the version of pROC to check installed:

``` r
rmarkdown::render("tools/binormal_validation.Rmd")
# or on one core (Windows):
rmarkdown::render("tools/binormal_validation.Rmd", params = list(cores = 1))
```

It takes a few minutes on 4 cores.

Sources:

- **D&A**: Dorfman DD, Alf E (1969) Maximum-likelihood estimation of
  parameters of signal-detection theory and determination of confidence
  intervals — rating-method data. *J Math Psychol* 6:487–496.
- **D&B**: Dorfman DD, Berbaum KS (1995) Degeneracy and discrete
  receiver operating characteristic rating data. *Acad Radiol*
  2:907–915.
- **O&M**: Obuchowski NA, McClish DK (1997) Sample size determination
  for diagnostic accuracy studies involving binormal ROC curve indices.
  *Stat Med* 16:1529–1542.

``` r
library(pROC)
#> Type 'citation("pROC")' for a citation.
#> 
#> Attaching package: 'pROC'
#> The following objects are masked from 'package:stats':
#> 
#>     cov, smooth, var
library(parallel)
B <- params$B
cores <- params$cores
RNGkind("L'Ecuyer-CMRG")
set.seed(20261106)
c(pROC = as.character(packageVersion("pROC")), R = R.version.string)
#>                           pROC                              R 
#>                    "1.19.9999" "R version 4.6.2 (2026-10-29)"
```

## Design

D&B (p. 909) sample five-category datasets from the population of Metz:
the binormal curve with $a = 1.5$ and $b = 0.8$, with the operating
points (FPF, TPF) = (0.00025, 0.10), (0.000328, 0.25), (0.03040, 0.50)
and (0.28114, 0.85). The FPF 0.000328 is not on the curve: a TPF of 0.25
needs an FPF of 0.00328. The cut-points are taken from the TPFs, as in
the report of the proper fit. Each dataset is fitted with
`method = "binormal"`; the least-squares line on normal-deviate axes
(the binormal smoothing of pROC 1.19, and the start of the fit) is shown
for comparison.

ROCFIT gives an exact fit with infinite parameters for degenerate
datasets and “does not compute diagnostic accuracy in such situations”
(D&B p. 907); D&B report its means on the non-degenerate datasets only.
pROC refuses the same data with an error (fewer than two points inside
the square, or points all on one horizontal or vertical line), so its
means are taken over the datasets it fits.

``` r
a <- 1.5
b <- 0.8
cuts <- sort((a - qnorm(c(0.10, 0.25, 0.50, 0.85))) / b)
true <- c(a = a, b = b, auc = pnorm(a / sqrt(1 + b^2)))
fpf.at <- c(0.05, 0.10, 0.20)
pa.bounds <- c(0, 0.25, 0.5, 0.75, 1)
# partial areas under TPF = pnorm(a + b * qnorm(FPF)) between the FPF bounds
partial.areas <- function(a, b) {
  vapply(1:4, function(i) {
    integrate(function(x) pnorm(a + b * qnorm(x)), pa.bounds[i], pa.bounds[i + 1])$value
  }, numeric(1))
}
one <- function(n.controls, n.cases) {
  r <- roc(
    controls = findInterval(rnorm(n.controls), cuts), cases = findInterval(rnorm(n.cases, a / b, 1 / b), cuts),
    direction = "<", quiet = TRUE
  )
  s <- tryCatch(smooth(r, method = "binormal"), error = function(e) NULL)
  if (is.null(s)) {
    return(rep(NA_real_, 14))
  }
  ci <- ci.auc(s)
  ml <- s$binormal
  df <- data.frame(sp = qnorm(r$specificities), se = qnorm(r$sensitivities))
  df <- df[is.finite(df$sp) & is.finite(df$se), ]
  ls <- pROC:::binormal_ls(df)
  c(
    a = ml$a, b = ml$b, auc = as.numeric(auc(s)),
    tpf = pnorm(ml$a + ml$b * qnorm(fpf.at)), pa = partial.areas(ml$a, ml$b),
    ls.a = ls[["a"]], ls.b = ls[["b"]],
    se = sqrt(var(s)), covered = ci[1] <= true[["auc"]] && true[["auc"]] <= ci[3]
  )
}
summarize <- function(x, true, published, published.rmse) {
  fitted <- !is.na(x[, 1])
  data.frame(
    true = true, mean = colMeans(x[fitted, , drop = FALSE]), published = published,
    rmse = sqrt(colMeans(sweep(x[fitted, , drop = FALSE], 2, true)^2)), published.rmse = published.rmse
  )
}
```

## Small datasets (D&B Tables 1, 3, 4 and 5)

1000 datasets of 30 negative and 40 positive cases in D&B; 2000 here.

``` r
small <- do.call(rbind, mclapply(seq_len(B), function(i) one(30, 40), mc.cores = cores))
refused <- is.na(small[, 1])
c(refused = mean(refused), rocfit.degenerate = 419 / 1000)
#>           refused rocfit.degenerate 
#>             0.403             0.419
```

pROC refuses 40.3% of the datasets; ROCFIT classified 419 of 1000 as
degenerate (D&B Table 1, p. 911).

The parameters and the area (D&B Table 3, p. 911, “Metz nondegenerate
samples (n = 581)”), the TPF at fixed FPFs (Table 4, p. 912) and the
partial areas (Table 5, p. 912), with ROCFIT’s means and RMSEs on its
581 non-degenerate datasets:

``` r
knitr::kable(summarize(small[, 1:10],
  true = c(true, pnorm(a + b * qnorm(fpf.at)), partial.areas(a, b)),
  published = c(1.631, 0.948, 0.869, 0.528, 0.656, 0.785, 0.163, 0.221, 0.238, 0.246),
  published.rmse = c(0.5474, 0.3505, 0.0500, 0.1171, 0.1020, 0.0868, 0.0226, 0.0170, 0.0111, 0.0053)
), digits = 4)
```

|      |   true |   mean | published |   rmse | published.rmse |
|:-----|-------:|-------:|----------:|-------:|---------------:|
| a    | 1.5000 | 1.6276 |     1.631 | 0.5323 |         0.5474 |
| b    | 0.8000 | 0.9362 |     0.948 | 0.3447 |         0.3505 |
| auc  | 0.8793 | 0.8703 |     0.869 | 0.0468 |         0.0500 |
| tpf1 | 0.5730 | 0.5343 |     0.528 | 0.1124 |         0.1171 |
| tpf2 | 0.6825 | 0.6602 |     0.656 | 0.0966 |         0.1020 |
| tpf3 | 0.7958 | 0.7880 |     0.785 | 0.0827 |         0.0868 |
| pa1  | 0.1692 | 0.1637 |     0.163 | 0.0216 |         0.0226 |
| pa2  | 0.2225 | 0.2215 |     0.221 | 0.0161 |         0.0170 |
| pa3  | 0.2398 | 0.2384 |     0.238 | 0.0103 |         0.0111 |
| pa4  | 0.2479 | 0.2466 |     0.246 | 0.0049 |         0.0053 |

The least-squares line on the same datasets:

``` r
knitr::kable(summarize(small[, 11:12], true = true[1:2], published = c(NA, NA), published.rmse = c(NA, NA))[, c("true", "mean", "rmse")], digits = 4)
```

|     | true |   mean |   rmse |
|:----|-----:|-------:|-------:|
| a   |  1.5 | 1.7679 | 0.7254 |
| b   |  0.8 | 1.0921 | 0.5640 |

## Large datasets (D&B Table 2)

2000 negative and 2000 positive cases (D&B Table 2, p. 911: 1000
datasets).

``` r
large <- do.call(rbind, mclapply(seq_len(B), function(i) one(2000, 2000), mc.cores = cores))
sum(is.na(large[, 1]))
#> [1] 0
knitr::kable(summarize(large[, 1:3],
  true = true, published = c(1.499, 0.799, 0.879), published.rmse = c(0.0541, 0.0398, 0.0059)
), digits = 4)
```

|     |   true |   mean | published |   rmse | published.rmse |
|:----|-------:|-------:|----------:|-------:|---------------:|
| a   | 1.5000 | 1.5002 |     1.499 | 0.0545 |         0.0541 |
| b   | 0.8000 | 0.8002 |     0.799 | 0.0386 |         0.0398 |
| auc | 0.8793 | 0.8791 |     0.879 | 0.0060 |         0.0059 |

## Standard error of the AUC (delta method)

`var()` and `ci.auc()` of the fit use the delta method: O&M Eq. 6
(p. 1531), with the variances and covariance of $a$ and $b$ from the
inverse of the expected information at the maximum (D&A p. 490), and a
Wald interval. The mean and median of the estimated SE are compared with
the empirical SD of the AUC estimates, with the coverage of the 95% CI,
on the datasets fitted above:

``` r
se.table <- function(x) {
  x <- x[!is.na(x[, 1]), , drop = FALSE]
  data.frame(
    datasets = nrow(x), sd = sd(x[, "auc"]), mean.se = mean(x[, "se"]), median.se = median(x[, "se"]),
    coverage = mean(x[, "covered"])
  )
}
knitr::kable(rbind(small = se.table(small), large = se.table(large)), digits = 4)
```

|       | datasets |     sd | mean.se | median.se | coverage |
|:------|---------:|-------:|--------:|----------:|---------:|
| small |     1194 | 0.0459 |  0.0472 |    0.0471 |   0.9439 |
| large |     2000 | 0.0060 |  0.0062 |    0.0062 |   0.9565 |

## Discussion

**Degenerate data.** pROC refuses about as many datasets as ROCFIT
classifies as degenerate (40.3% against 41.9%, 0.8 SE apart). ROCFIT’s
classes of degeneracy are only given in its code (D&B p. 907), and D&B
may have generated the counts at the published operating points rather
than at the cut-points of the curve.

**Small datasets.** On the datasets pROC fits, the means of $a$, $b$,
the area, the TPFs and the partial areas agree with ROCFIT’s on its
non-degenerate datasets within 1.6 SE (counting the simulation error of
both, and the rounding of the published means). With so few cases both
overestimate $a$ and $b$ and slightly underestimate the area. The RMSEs
of $a$, $b$ and the area are 2% to 6% lower than ROCFIT’s.

**Large datasets.** The means and RMSEs agree with ROCFIT’s.

**Standard error.** The mean SE of the delta method is within 2.8% of
the empirical SD of the AUC estimates, and the coverage of the 95% Wald
interval is 0.944 on the small datasets and 0.957 on the large ones.

**Least squares.** The least-squares line, the binormal smoothing of
pROC 1.19 and the start of the fit, overestimates $a$ and $b$ more
(means 1.768 and 1.092), with larger RMSEs (0.725 and 0.564). It is not
a maximum-likelihood estimate: D&A (p. 491) note that an “eyeball” fit
on normal-deviate axes “would not coincide with the fit based upon the
MLE”.
