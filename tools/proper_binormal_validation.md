
- [Validation of the proper binormal
  fit](#validation-of-the-proper-binormal-fit)
  - [AUC (M&P Eq. 36)](#auc-mp-eq-36)
  - [Continuous data (P&M Tables 2–4)](#continuous-data-pm-tables-24)
  - [Five-category data (M&P Tables
    1–3)](#five-category-data-mp-tables-13)
  - [Perfect performance](#perfect-performance)
  - [Comparison with PROPROC fits](#comparison-with-proproc-fits)
  - [Discussion](#discussion)

<!-- proper_binormal_validation.md is generated from proper_binormal_validation.Rmd. Please edit that file -->

# Validation of the proper binormal fit

This report checks `smooth(method = "proper")` (issue
\#152) against the published results of the “proper” binormal model
(PBM). It is maintainer-only: knit it from the package root with the
version of pROC to check installed, and mvtnorm:

``` r
rmarkdown::render("tools/proper_binormal_validation.Rmd")
# or with fewer datasets, on one core (Windows):
rmarkdown::render("tools/proper_binormal_validation.Rmd", params = list(B = 200, cores = 1))
```

With `B = 2000` datasets per simulated population it takes about 2 hours
on 4 cores.

Sources:

- **M&P**: Metz CE, Pan X (1999) “Proper” binormal ROC curves: theory
  and maximum-likelihood estimation. *J Math Psychol* 43:1–33.
- **P&M**: Pesce LL, Metz CE (2007) Reliable and computationally
  efficient maximum-likelihood estimation of “proper” binormal ROC
  curves. *Acad Radiol* 14:814–829.
- **Owen**: Owen DB (1956) Tables for computing bivariate normal
  probabilities. *Ann Math Stat* 27:1075–1090.

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
set.seed(20261103)
c(pROC = as.character(packageVersion("pROC")), mvtnorm = as.character(packageVersion("mvtnorm")), R = R.version.string)
#>                           pROC                        mvtnorm 
#>                    "1.19.9999"                        "1.2.4" 
#>                              R 
#> "R version 4.3.3 (2024-02-29)"
```

## AUC (M&P Eq. 36)

The AUC of the PBM is computed with M&P p. 12, Eq. 36. P&M Table 1
(p. 826) lists the curve parameters of 9 population curves of known AUC,
to 6 significant digits:

``` r
table1 <- data.frame(
  auc = rep(c(0.6, 0.85, 0.95), each = 3),
  shape = rep(c("high", "medium", "symmetric"), 3),
  da = c(0.1, 0.358, 0.358287, 1.24365, 1.46571, 1.46574, 0.0999998, 2.324, 2.32618),
  c = c(-0.154357, -0.046907, 0, -0.428571, -0.14, 0, -0.853738, -0.292369, 0)
)
table1$eq36 <- mapply(pROC:::proper_binormal_auc, table1$da, table1$c)
#> Loading required namespace: mvtnorm
knitr::kable(table1, digits = 8)
```

|  auc | shape     |        da |         c |      eq36 |
|-----:|:----------|----------:|----------:|----------:|
| 0.60 | high      | 0.1000000 | -0.154357 | 0.6000000 |
| 0.60 | medium    | 0.3580000 | -0.046907 | 0.6000000 |
| 0.60 | symmetric | 0.3582870 |  0.000000 | 0.6000000 |
| 0.85 | high      | 1.2436500 | -0.428571 | 0.8500000 |
| 0.85 | medium    | 1.4657100 | -0.140000 | 0.8499997 |
| 0.85 | symmetric | 1.4657400 |  0.000000 | 0.8500003 |
| 0.95 | high      | 0.0999998 | -0.853738 | 0.9500000 |
| 0.95 | medium    | 2.3240000 | -0.292369 | 0.9500000 |
| 0.95 | symmetric | 2.3261800 |  0.000000 | 0.9500004 |

The largest difference to the published AUCs is 4.2^{-7}, within the
rounding of the parameters.

At $d_a = 0$, $F(0, 0; \rho) = 1/4 + \arcsin(\rho) / (2\pi)$ (Owen
p. 1078), so the AUC is $1 + \arcsin(\rho) / \pi$ with
$\rho = -(1 - c^2) / (1 + c^2)$:

``` r
cs <- c(-0.99, -0.5, -1e-3, -1e-6, 1e-6, 1e-3, 0.5, 0.99)
closed <- 1 + asin(-(1 - cs^2) / (1 + cs^2)) / pi
max(abs(vapply(cs, pROC:::proper_binormal_auc, numeric(1), da = 0) - closed))
#> [1] 1.110223e-16
```

## Continuous data (P&M Tables 2–4)

P&M p. 822 generate continuous data from the PBM: $x$ from the
conventional binormal model with parameters $a$ and $b$ (negatives
$N(0, 1)$, positives $N(a / b, 1 / b)$), and
$y = -[(bx - a)^2 - x^2] / 2$, the decision variable of the PBM (M&P Eq.
11). $a$ and $b$ follow from $d_a$ and $c$ (M&P p. 23). At $c = 0$
($b = 1$), $y$ is linear in $x$, which is used directly.

For each of the 9 population curves of Table 1, 2000 datasets of 50
negative and 50 positive cases are fitted with `method = "proper"`, and,
for comparison, with the least-squares binormal fit (`method = "binormal"`, the
default) and the empirical AUC. Data with perfect performance (see
below) give the limit AUC = 1, which is kept in the means as an
estimate, as PROPROC returns it; `perfect` counts them.

``` r
sim.pbm <- function(n, da, c) {
  b <- (1 + c) / (1 - c)
  a <- da * sqrt(1 + c^2) / (1 - c)
  y <- if (c == 0) identity else function(x) -((b * x - a)^2 - x^2) / 2
  list(controls = y(rnorm(n)), cases = y(rnorm(n, a / b, 1 / b)))
}
fit.auc <- function(d) {
  r <- roc(controls = d$controls, cases = d$cases, direction = "<", quiet = TRUE)
  s <- tryCatch(smooth(r, method = "proper"), error = function(e) NULL)
  c(
    proper = if (is.null(s)) NA_real_ else as.numeric(auc(s)),
    # perfect performance: the limit AUC = 1, without a log-likelihood
    perfect = !is.null(s) && is.na(s$proper$logLik),
    lm = tryCatch(as.numeric(auc(smooth(r))), error = function(e) NA_real_),
    empirical = as.numeric(auc(r))
  )
}
# P&M Tables 2-4 (pp. 826-827), 50 positive and 50 negative cases: mean AUC
# and empirical SD of the estimates
published <- data.frame(
  mean = c(0.617, 0.610, 0.608, 0.860, 0.852, 0.852, 0.957, 0.951, 0.951),
  sd = c(0.049, 0.050, 0.051, 0.034, 0.038, 0.038, 0.018, 0.021, 0.019)
)
continuous <- lapply(seq_len(nrow(table1)), function(i) {
  do.call(rbind, mclapply(seq_len(B), function(j) {
    fit.auc(sim.pbm(50, table1$da[i], table1$c[i]))
  }, mc.cores = cores))
})
summary.continuous <- do.call(rbind, lapply(seq_along(continuous), function(i) {
  p <- continuous[[i]][, "proper"]
  n <- sum(!is.na(p))
  data.frame(
    auc = table1$auc[i], shape = table1$shape[i],
    mean = mean(p, na.rm = TRUE), se = sd(p, na.rm = TRUE) / sqrt(n),
    published = published$mean[i],
    sd = sd(p, na.rm = TRUE), published.sd = published$sd[i],
    perfect = sum(continuous[[i]][, "perfect"]), failures = sum(is.na(p))
  )
}))
summary.continuous$z <- (summary.continuous$mean - summary.continuous$published) / summary.continuous$se
knitr::kable(summary.continuous, digits = c(2, 0, 4, 4, 3, 4, 3, 0, 0, 1))
```

|  auc | shape     |   mean |     se | published |     sd | published.sd | perfect | failures |   z |
|-----:|:----------|-------:|-------:|----------:|-------:|-------------:|--------:|---------:|----:|
| 0.60 | high      | 0.6179 | 0.0011 |     0.617 | 0.0483 |        0.049 |       0 |        0 | 0.8 |
| 0.60 | medium    | 0.6120 | 0.0011 |     0.610 | 0.0511 |        0.050 |       0 |        0 | 1.8 |
| 0.60 | symmetric | 0.6095 | 0.0011 |     0.608 | 0.0514 |        0.051 |       0 |        0 | 1.3 |
| 0.85 | high      | 0.8618 | 0.0008 |     0.860 | 0.0336 |        0.034 |       0 |        0 | 2.4 |
| 0.85 | medium    | 0.8533 | 0.0008 |     0.852 | 0.0380 |        0.038 |       0 |        0 | 1.6 |
| 0.85 | symmetric | 0.8537 | 0.0008 |     0.852 | 0.0375 |        0.038 |       0 |        0 | 2.0 |
| 0.95 | high      | 0.9576 | 0.0004 |     0.957 | 0.0187 |        0.018 |       1 |        0 | 1.5 |
| 0.95 | medium    | 0.9515 | 0.0005 |     0.951 | 0.0213 |        0.021 |       0 |        0 | 1.0 |
| 0.95 | symmetric | 0.9512 | 0.0004 |     0.951 | 0.0195 |        0.019 |       0 |        0 | 0.4 |

### Bias, and comparison with the least-squares fit

``` r
bias <- do.call(rbind, lapply(seq_along(continuous), function(i) {
  res <- continuous[[i]]
  true <- table1$auc[i]
  rmse <- function(x) sqrt(mean((x - true)^2, na.rm = TRUE))
  data.frame(
    auc = true, shape = table1$shape[i],
    bias.proper = mean(res[, "proper"], na.rm = TRUE) - true,
    bias.published = published$mean[i] - true,
    bias.lm = mean(res[, "lm"], na.rm = TRUE) - true,
    bias.empirical = mean(res[, "empirical"]) - true,
    rmse.proper = rmse(res[, "proper"]), rmse.lm = rmse(res[, "lm"]), rmse.empirical = rmse(res[, "empirical"]),
    failures.lm = sum(is.na(res[, "lm"]))
  )
}))
knitr::kable(bias, digits = c(2, 0, 4, 3, 4, 4, 4, 4, 4, 0))
```

|  auc | shape     | bias.proper | bias.published | bias.lm | bias.empirical | rmse.proper | rmse.lm | rmse.empirical | failures.lm |
|-----:|:----------|------------:|---------------:|--------:|---------------:|------------:|--------:|---------------:|------------:|
| 0.60 | high      |      0.0179 |          0.017 |  0.0013 |         0.0006 |      0.0515 |  0.0547 |         0.0562 |           0 |
| 0.60 | medium    |      0.0120 |          0.010 |  0.0019 |         0.0030 |      0.0525 |  0.0557 |         0.0571 |           0 |
| 0.60 | symmetric |      0.0095 |          0.008 |  0.0015 |         0.0028 |      0.0522 |  0.0555 |         0.0569 |           0 |
| 0.85 | high      |      0.0118 |          0.010 | -0.0010 |         0.0017 |      0.0357 |  0.0378 |         0.0389 |           0 |
| 0.85 | medium    |      0.0033 |          0.002 | -0.0028 |         0.0017 |      0.0381 |  0.0379 |         0.0384 |           0 |
| 0.85 | symmetric |      0.0037 |          0.002 | -0.0030 |         0.0017 |      0.0377 |  0.0375 |         0.0377 |           0 |
| 0.95 | high      |      0.0076 |          0.007 | -0.0044 |         0.0001 |      0.0202 |  0.0244 |         0.0242 |          12 |
| 0.95 | medium    |      0.0015 |          0.001 | -0.0041 |         0.0006 |      0.0213 |  0.0224 |         0.0223 |           1 |
| 0.95 | symmetric |      0.0012 |          0.001 | -0.0047 |         0.0006 |      0.0195 |  0.0207 |         0.0199 |           0 |

## Five-category data (M&P Tables 1–3)

M&P sample 1000 datasets of 30 negative and 40 positive cases (p. 21)
from the conventional binormal model with $a = 1.5$ and $b = 0.8$, in 5
categories with the population operating points (p. 20) (FPF, TPF) =
(0.00025, 0.10), (0.000328, 0.25), (0.03040, 0.50) and (0.28114, 0.85).
The FPF 0.000328 is not on the curve: a TPF of 0.25 needs an FPF of
0.00328. The cut-points are taken from the TPFs.

``` r
a <- 1.5
b <- 0.8
cuts <- sort((a - qnorm(c(0.10, 0.25, 0.50, 0.85))) / b)
fpf.at <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
true.tpf <- pnorm(a + b * qnorm(fpf.at))
true.area <- pnorm(a / sqrt(1 + b^2))
categorical <- do.call(rbind, mclapply(seq_len(B), function(j) {
  r <- roc(
    controls = findInterval(rnorm(30), cuts), cases = findInterval(rnorm(40, a / b, 1 / b), cuts),
    direction = "<", quiet = TRUE
  )
  lm.area <- tryCatch(as.numeric(auc(smooth(r))), error = function(e) NA_real_)
  s <- tryCatch(smooth(r, method = "proper"), error = function(e) NULL)
  if (is.null(s)) {
    return(c(rep(NA_real_, 8), lm.area, FALSE))
  }
  c(
    approx(1 - s$specificities, s$sensitivities, xout = fpf.at, ties = max)$y, as.numeric(auc(s)), lm.area,
    is.na(s$proper$logLik)
  )
}, mc.cores = cores))
fitted <- !is.na(categorical[, 8])
knitr::kable(data.frame(
  FPF = fpf.at, true.TPF = true.tpf, mean.TPF = colMeans(categorical[fitted, 1:7]),
  published.mean = c(0.615, 0.705, 0.802, 0.857, 0.921, 0.960, 0.988),
  RMSE = sqrt(colMeans(sweep(categorical[fitted, 1:7], 2, true.tpf)^2)),
  published.RMSE = c(0.1557, 0.1150, 0.0795, 0.0639, 0.0459, 0.0297, 0.0111)
), digits = c(2, 3, 3, 3, 4, 4))
```

|  FPF | true.TPF | mean.TPF | published.mean |   RMSE | published.RMSE |
|-----:|---------:|---------:|---------------:|-------:|---------------:|
| 0.05 |    0.573 |    0.617 |          0.615 | 0.1528 |         0.1557 |
| 0.10 |    0.683 |    0.707 |          0.705 | 0.1128 |         0.1150 |
| 0.20 |    0.796 |    0.803 |          0.802 | 0.0788 |         0.0795 |
| 0.30 |    0.860 |    0.857 |          0.857 | 0.0638 |         0.0639 |
| 0.50 |    0.933 |    0.921 |          0.921 | 0.0460 |         0.0459 |
| 0.70 |    0.973 |    0.960 |          0.960 | 0.0297 |         0.0297 |
| 0.90 |    0.994 |    0.988 |          0.988 | 0.0111 |         0.0111 |

The area under the proper fits averages 0.8799 (RMSE 0.0438), for a true
area of 0.8793 (M&P Table 3, p. 23: 0.879, RMSE 0.0437); 1 of 2000
datasets had perfect performance (AUC = 1), and 0 could not be fitted.
The least-squares fit refuses 777 of them (degenerate data, or fewer
than two points), with an area of 0.8674 (RMSE 0.0463) on the others.

## Perfect performance

No curve parameters can be estimated from perfectly separated data (M&P
p. 27), nor, more generally, from data whose operating points all have
FPF = 0 or TPF = 1 (for instance ties straddling a separation): the
likelihood is maximal only in the limit of an AUC of 1. As PROPROC (see
below), pROC then returns that limit, $d_a = 0$ and $c = 1$, with a
warning. Without fitting, the probability that all the controls are
below all the cases, and, for the five-category data, the probability of
the more general case:

``` r
N <- 2e5
high <- table1[table1$auc == 0.95 & table1$shape == "high", ]
b.high <- (1 + high$c) / (1 - high$c)
a.high <- high$da * sqrt(1 + high$c^2) / (1 - high$c)
y.high <- function(x) -((b.high * x - a.high)^2 - x^2) / 2
p.high <- mean(replicate(N, max(y.high(rnorm(50))) < min(y.high(rnorm(50, a.high / b.high, 1 / b.high)))))
p.categorical <- mean(replicate(N, max(findInterval(rnorm(30), cuts)) < min(findInterval(rnorm(40, a / b, 1 / b), cuts))))
perfect <- function() {
  k <- tabulate(findInterval(rnorm(30), cuts) + 1, 5)
  l <- tabulate(findInterval(rnorm(40, a / b, 1 / b), cuts) + 1, 5)
  all(rev(cumsum(rev(k)))[-1] == 0 | cumsum(l)[-5] == 0)
}
p.perfect <- mean(replicate(N, perfect()))
c(continuous.0.95.high = p.high, five.category = p.categorical, five.category.perfect = p.perfect) * B
#>  continuous.0.95.high         five.category five.category.perfect 
#>                  0.99                  0.00                  1.15
```

These are the expected numbers of such datasets per 2000.

## Comparison with PROPROC fits

PROPROC itself is no longer distributed, but the repository of the book
*The RJafroc Book: ROC analysis* by Chakraborty
(dpc10ster/RJafrocRocBook on GitHub) keeps the output of PROPROC (OR-DBM
MRMC 2.5, version of 2014) for 236 fits, every modality and reader of 14
reader-study datasets, with the input files
(`R/compare-3-fits/MRMCRuns`). They are fitted here again with
`method = "proper"`. For each fit, the log-likelihood that PROPROC’s
$(d_a, c)$ can reach is its maximum over the cut-points (`ll.proproc`),
which is compared with pROC’s maximum (`ll`).

``` r
commit <- "98aa3ac0a629f63cbb86ca86761e509c80955721"
url <- sprintf("https://raw.githubusercontent.com/dpc10ster/RJafrocRocBook/%s/R/compare-3-fits/MRMCRuns", commit)
datasets <- c("DOB1", "DOB2", "DOB3", "FED", "FR", "FZR", "JT", "MAG", "NICO", "OPT", "PEN", "RUS", "TONY", "VD")
# .lrc files: a title line, the ID of the first reader, the modality names,
# the direction of each modality (L: larger values are more positive); then
# per reader the controls (one column per modality), "*", the cases, "*",
# and the ID of the next reader, or "#" at the end
read.lrc <- function(lines) {
  stopifnot(all(strsplit(trimws(lines[4]), "\\s+")[[1]] == "L"))
  readers <- list()
  i <- 5
  block <- function() {
    end <- i + which(trimws(lines[i:length(lines)]) == "*")[1] - 2
    values <- strsplit(trimws(lines[i:end]), "\\s+")
    i <<- end + 2
    do.call(rbind, lapply(values, function(x) as.numeric(x[!grepl("[A-Za-z]", x)])))
  }
  repeat {
    readers[[length(readers) + 1]] <- list(controls = block(), cases = block())
    if (trimws(lines[i]) == "#") break
    i <- i + 1
  }
  readers
}
proproc <- do.call(rbind, lapply(datasets, function(d) {
  data <- read.lrc(readLines(sprintf("%s/%s/%s.lrc", url, d, d)))
  out <- read.csv(sprintf("%s/%s/%sproprocareapooled.csv", url, d, d), strip.white = TRUE)
  do.call(rbind, lapply(seq_len(nrow(out)), function(j) {
    r <- roc(
      controls = data[[out$R[j]]]$controls[, out$T[j]], cases = data[[out$R[j]]]$cases[, out$T[j]],
      direction = "<", quiet = TRUE
    )
    categories <- pROC:::proper_binormal_categories(r)
    s <- smooth(r, method = "proper")
    data.frame(
      dataset = d, modality = out$T[j], reader = out$R[j], code = out$returnCode[j],
      categories.proproc = out$numCAT[j], categories = length(categories$controls),
      auc.proproc = out$area[j], auc = as.numeric(auc(s)),
      da.proproc = out$d_a[j], da = s$proper$da, c.proproc = out$c[j], c = s$proper$c,
      ll.proproc = if (out$returnCode[j] == 3) NA else
        pROC:::proper_binormal_profile(out$d_a[j], out$c[j], categories$controls, categories$cases),
      ll = s$proper$logLik
    )
  }))
}))
#> Warning: Perfect performance: the proper binormal fit is the limit of AUC = 1
#> (d_a = 0, c = 1).
proproc$ll.diff <- proproc$ll - proproc$ll.proproc
proproc$auc.diff <- proproc$auc - proproc$auc.proproc
c(
  fits = nrow(proproc), same.categories = sum(proproc$categories == proproc$categories.proproc),
  same.auc.1e.4 = sum(abs(proproc$auc.diff) < 1e-4),
  ll.lower = sum(proproc$ll.diff < -1e-4, na.rm = TRUE), ll.higher = sum(proproc$ll.diff > 1e-4, na.rm = TRUE),
  da.below.0.01.proproc = sum(proproc$da.proproc < 0.01), da.below.0.01 = sum(proproc$da < 0.01)
)
#>                  fits       same.categories         same.auc.1e.4 
#>                   236                   236                   235 
#>              ll.lower             ll.higher da.below.0.01.proproc 
#>                     1                     0                    73 
#>         da.below.0.01 
#>                    10
summary(proproc$ll.diff)
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.       NA's 
#> -0.0021358  0.0000000  0.0000000 -0.0000091  0.0000000  0.0000000          1
summary(abs(proproc$auc.diff))
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.000e+00 1.000e-08 5.000e-08 1.864e-04 1.200e-07 4.395e-02
mean(proproc$auc.diff)
#> [1] -0.000186341
```

Eq. 36 (pROC’s AUC) at PROPROC’s own estimates, against PROPROC’s AUC:

``` r
max(abs(mapply(pROC:::proper_binormal_auc, proproc$da.proproc, proproc$c.proproc) - proproc$auc.proproc))
#> [1] 4.636661e-06
```

The fits whose AUC differs by more than 1e-4, at the chance line, or
with return code 3 (perfect separation):

``` r
knitr::kable(
  proproc[which(abs(proproc$auc.diff) > 1e-4 | abs(proproc$auc - 0.5) < 1e-3 | proproc$code == 3), c(
    "dataset", "modality", "reader", "code", "auc.proproc", "auc", "da.proproc", "da",
    "c.proproc", "c", "ll.diff"
  )],
  digits = c(0, 0, 0, 0, 4, 4, 4, 4, 3, 3, 3), row.names = FALSE
)
```

| dataset | modality | reader | code | auc.proproc |    auc | da.proproc |     da | c.proproc |     c | ll.diff |
|:--------|---------:|-------:|-----:|------------:|-------:|-----------:|-------:|----------:|------:|--------:|
| DOB3    |        1 |      2 |    6 |      0.5000 | 0.5000 |     0.0000 | 0.0000 |     0.000 | 0.000 |   0.000 |
| DOB3    |        3 |      5 |    6 |      0.8467 | 0.8028 |     0.0285 | 1.2039 |     0.606 | 0.141 |  -0.002 |
| VD      |        2 |      4 |    3 |      1.0000 | 1.0000 |     0.0000 | 0.0000 |     1.000 | 1.000 |      NA |

## Discussion

**Eq. 36.** The AUCs of the 9 population curves of P&M Table 1 agree to
the rounding of their parameters, and the closed form at $d_a = 0$ to
machine precision.

**Continuous data.** All 9 means are above the published ones, by 0.0002
to 0.0020 (0.4 to 2.4 SE; combined $z = \sum z / \sqrt{9} =$ 4.3). The
published means come from 20,000 datasets each (P&M p. 826, note to
Table 2) and are rounded to 3 decimals. With their simulation error, the
combined $z$ is 4.1; adding the rounding as a uniform error of
$\pm 0.0005$, 3.8; and if all the published means were rounded down by
0.0005, 2.0.

The difference may still be systematic, but it is at most 0.05 SD of the
estimates. It is unlikely to come from the fit: on the reader-study data
below, pROC reproduces the fits of a later PROPROC, with a mean AUC
difference of -0.00019. Nor from the AUC: Eq. 36 at PROPROC’s estimates
gives its AUCs within 5e-06. Nor from the categories: those data have
few categories, but P&M’s algorithm uses the truth-run categories of
continuous data unchanged up to 400 categories (p. 817), and 50 + 50
cases give at most 100, as in pROC. The original PROPROC collapsed them
into at most 10 or 20, which P&M found “biases the estimates in a
complicated, almost unpredictable, way” (p. 817). What remains is a
detail of the data generation that P&M do not give, or chance. The SDs
agree with the published ones within 4%.

**Bias.** The proper fit overestimates the AUC by up to 0.018, most at
low AUC and high skew, like the published fits. P&M (pp. 825–826)
attribute it to the lower bound of 0.5 on the AUC of proper models and
to the bounds on the curve shape. The least-squares fit and the
empirical AUC are nearly unbiased on these data (within 0.005). The
proper fit nevertheless has the lowest RMSE in 7 of the 9 populations,
and is within 1% of the least-squares fit in the others (AUC 0.85,
symmetric and medium skew). Its advantage is largest at AUC 0.95 with
high skew, where the least-squares fit also refuses 12 datasets.

**Five-category data.** The mean TPFs agree with M&P Table 1 within
0.002, their RMSEs with Table 2 within 2%, and the area with Table 3.
The least-squares fit refuses 777 of the 2000 datasets, and
underestimates the area of the others.

**Perfect performance and failures.** The proper fit returned the limit
AUC = 1 for 1 of the 18,000 continuous datasets and for 1 of the 2000
five-category datasets, about the numbers expected above. It failed on 0
continuous and 0 five-category datasets.

**PROPROC.** pROC finds the same categories as PROPROC in 236 of the 236
reader-study fits, and the same AUC within 1e-4 in 235, including the
fit at the chance line and the limit AUC = 1 of the perfectly separated
data (return code 3, which PROPROC also uses for these data with a tie
straddling the separation). PROPROC’s estimates of $d_a$ are below 0.01
in 73 fits (mostly about 1e-4, the fits with return code 6), pROC’s in
10: there the log-likelihood is nearly flat in $d_a$ (P&M p. 820), and
both reach the same maximum. In 1 fit, PROPROC finds a higher maximum
(0.0021 in log-likelihood) that pROC’s starts miss, with an AUC 0.044
higher.

Not covered here: on harder data with many categories, the code review
found maxima that the starts miss, in 2 of 40 samples, with
log-likelihoods 0.05 above the fit.
