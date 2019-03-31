# Benchmark code adapted from the cutpointr vignette by Christian Thiele
# https://github.com/Thie1e/cutpointr/blob/master/vignettes/cutpointr.Rmd
# Licence: GPL-3

library(cutpointr)
library(ggplot2)

# Return cutpoint that maximizes the sum of sensitivity and specificiy
# ROCR package
rocr_sensspec <- function(x, class) {
	pred <- ROCR::prediction(x, class)
	perf <- ROCR::performance(pred, "sens", "spec")
	sens <- slot(perf, "y.values")[[1]]
	spec <- slot(perf, "x.values")[[1]]
	cut <- slot(perf, "alpha.values")[[1]]
	cut[which.max(sens + spec)]
}

# pROC package
proc_sensspec <- function(x, class, 
						  levels = c("no", "yes"), algo = 2) {
	r <- pROC::roc(class, x, algorithm = algo, quiet=TRUE)
	sens <- r$sensitivities
	spec <- r$specificities
	cut <- r$thresholds
	cut[which.max(sens + spec)]
}

n <- 1000
set.seed(123)
dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
x_pos <- dat$x[dat$y == 1]
x_neg <- dat$x[dat$y == 0]
bench_1000 <- microbenchmark::microbenchmark(
	cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
			  direction = ">=", metric = youden, break_ties = mean),
	rocr_sensspec(dat$x, dat$y),
	proc_sensspec(dat$x, dat$y, algo = 2)
)

n <- 10000
set.seed(123)
dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
x_pos <- dat$x[dat$y == 1]
x_neg <- dat$x[dat$y == 0]
bench_10000 <- microbenchmark::microbenchmark(
	cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
			  direction = ">=", metric = youden, break_ties = mean, silent = TRUE),
	rocr_sensspec(dat$x, dat$y),
	proc_sensspec(dat$x, dat$y, algo = 2),
	times = 20
)

n <- 1e5
set.seed(123)
dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
bench_1e5 <- microbenchmark::microbenchmark(
	cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
			  direction = ">=", metric = youden, break_ties = mean),
	rocr_sensspec(dat$x, dat$y),
	proc_sensspec(dat$x, dat$y, algo = 2),
	times = 20
)

n <- 1e6
set.seed(123)
dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
bench_1e6 <- microbenchmark::microbenchmark(
	cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
			  direction = ">=", metric = youden, break_ties = mean),
	rocr_sensspec(dat$x, dat$y),
	proc_sensspec(dat$x, dat$y, algo = 2),
	times = 10
)

n <- 1e7
set.seed(123)
dat <- data.frame(x = rnorm(n), y = sample(c(0:1), size = n, replace = TRUE))
bench_1e7 <- microbenchmark::microbenchmark(
	cutpointr(dat, x, y, pos_class = 1, neg_class = 0,
			  direction = ">=", metric = youden, break_ties = mean),
	rocr_sensspec(dat$x, dat$y),
	proc_sensspec(dat$x, dat$y, algo = 2),
	times = 10, unit = "ms"
)

results <- rbind(
	data.frame(time = summary(bench_1000)$median,
			   solution = summary(bench_1000)$expr, 
			   n = 1000),
	data.frame(time = summary(bench_10000)$median,
			   solution = summary(bench_10000)$expr, 
			   n = 10000),
	data.frame(time = summary(bench_1e5)$median,
			   solution = summary(bench_1e5)$expr, 
			   n = 1e5),
	data.frame(time = summary(bench_1e6)$median,
			   solution = summary(bench_1e6)$expr, 
			   n = 1e6),
	data.frame(time = summary(bench_1e7)$median,
			   solution = summary(bench_1e7)$expr, 
			   n = 1e7)
)
results$solution <- as.character(results$solution)
results$solution[grep(pattern = "cutpointr", x = results$solution)] <- "cutpointr"
results$solution[grep(pattern = "rocr", x = results$solution)] <- "ROCR"
results$solution[grep(pattern = "proc", x = results$solution)] <- "pROC"

ggplot(results, aes(x = n, y = time, col = solution, shape = solution)) +
	geom_point(size = 3) + geom_line() +
	scale_y_log10(breaks = c(3, 5, 10, 25, 100, 250, 1000, 5000, 1e4, 15000)) +
	scale_x_log10(breaks = c(1000, 1e4, 1e5, 1e6, 1e7)) +
	ggtitle("Benchmark results", "n = 1000, 10000, 1e5, 1e6, 1e7") +
	ylab("Median time (milliseconds, log scale)") + xlab("n (log scale)")
