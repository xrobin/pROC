library(pROC)
library(parallel)
data(aSAH)

# Parallel bootstrapping (issue #46). pROC never allocates a core on its own:
# the user supplies the cluster, or asks for a number of workers explicitly.
#
# CRAN allows at most two cores, so every cluster here has two workers and the
# tests are skipped when the check farm asks for a single core.

B <- 50

skip_if_no_cluster <- function() {
  skip_on_cran()
  # _R_CHECK_LIMIT_CORES_ is set on the CRAN check machines.
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    skip("Running under _R_CHECK_LIMIT_CORES_")
  }
}

with_cluster <- function(n, code) {
  cl <- parallel::makeCluster(n)
  on.exit(parallel::stopCluster(cl))
  code(cl)
}


test_that("a cluster gives the same answer however many workers it has", {
  skip_if_no_cluster()
  two <- with_cluster(2, function(cl) {
    set.seed(42)
    ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = cl)
  })
  # A second cluster of a different size must not change the result: each
  # replicate owns its RNG stream, so the worker that runs it is irrelevant.
  one <- with_cluster(1, function(cl) {
    set.seed(42)
    ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = cl)
  })
  expect_equal(as.numeric(two), as.numeric(one))
})


test_that("the same seed gives the same answer twice", {
  skip_if_no_cluster()
  with_cluster(2, function(cl) {
    set.seed(42)
    a <- ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = cl)
    set.seed(42)
    b <- ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = cl)
    expect_identical(a, b)
    set.seed(1)
    d <- ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = cl)
    expect_false(isTRUE(all.equal(as.numeric(a), as.numeric(d))))
  })
})


test_that("every bootstrap entry point accepts a cluster", {
  skip_if_no_cluster()
  with_cluster(2, function(cl) {
    expect_s3_class(ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = cl), "ci.auc")
    expect_s3_class(ci.se(r.s100b, specificities = c(0.1, 0.9), boot.n = B, cl = cl), "ci.se")
    expect_s3_class(ci.sp(r.s100b, sensitivities = c(0.1, 0.9), boot.n = B, cl = cl), "ci.sp")
    expect_s3_class(ci.thresholds(r.s100b, thresholds = c(0.1, 0.5), boot.n = B, cl = cl),
                    "ci.thresholds")
    expect_s3_class(ci.coords(r.s100b, x = c(0.1, 0.9), input = "specificity",
                              ret = "sensitivity", boot.n = B, cl = cl), "ci.coords")
    expect_true(is.finite(var(r.s100b, method = "bootstrap", boot.n = B, cl = cl)))
    expect_true(is.finite(cov(r.s100b, r.ndka, method = "bootstrap", boot.n = B, cl = cl)))
    expect_s3_class(roc.test(r.s100b, r.ndka, method = "bootstrap", boot.n = B, cl = cl),
                    "htest")
  })
})


test_that("smoothed bootstraps run on a cluster", {
  skip_if_no_cluster()
  skip_if_not_installed("MASS")
  s <- smooth(r.ndka)
  with_cluster(2, function(cl) {
    expect_s3_class(ci.auc(s, boot.n = B, cl = cl), "ci.auc")
    expect_s3_class(ci.se(s, specificities = c(0.1, 0.9), boot.n = B, cl = cl), "ci.se")
  })
})


test_that("an integer creates and stops a cluster", {
  skip_if_no_cluster()
  before <- length(getDefaultCluster())
  set.seed(42)
  by.number <- ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = 2)
  by.cluster <- with_cluster(2, function(cl) {
    set.seed(42)
    ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = cl)
  })
  expect_equal(as.numeric(by.number), as.numeric(by.cluster))
  # The cluster pROC made must not outlive the call.
  expect_identical(length(getDefaultCluster()), before)
})


test_that("cl = TRUE uses the registered default cluster", {
  skip_if_no_cluster()
  cl <- parallel::makeCluster(2)
  on.exit({
    parallel::setDefaultCluster(NULL)
    parallel::stopCluster(cl)
  })
  parallel::setDefaultCluster(cl)
  set.seed(42)
  by.default <- ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = TRUE)
  set.seed(42)
  explicit <- ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = cl)
  expect_equal(as.numeric(by.default), as.numeric(explicit))
})


test_that("cl = TRUE without a default cluster is an informative error", {
  parallel::setDefaultCluster(NULL)
  expect_error(ci.auc(r.s100b, method = "bootstrap", boot.n = 5, cl = TRUE),
               "needs a default cluster")
})


test_that("sequential values of cl stay sequential", {
  set.seed(42)
  plain <- ci.auc(r.s100b, method = "bootstrap", boot.n = B)
  for (value in list(NULL, FALSE, 1)) {
    set.seed(42)
    expect_identical(
      ci.auc(r.s100b, method = "bootstrap", boot.n = B, cl = value),
      plain
    )
  }
})


test_that("a nonsensical cl is rejected", {
  expect_error(ci.auc(r.s100b, method = "bootstrap", boot.n = 5, cl = "snow"), "'cl' must be")
  expect_error(ci.auc(r.s100b, method = "bootstrap", boot.n = 5, cl = list(1, 2)), "'cl' must be")
  expect_error(ci.auc(r.s100b, method = "bootstrap", boot.n = 5, cl = -1), "'cl' must be")
})


test_that("RNG streams are one per replicate, not one per worker", {
  # The distinguishing property: a replicate's stream depends only on its
  # index, so results survive a change in the number of workers. This is what
  # clusterSetRNGStream() -- one stream per worker -- cannot promise.
  set.seed(42)
  streams <- roc_utils_rng_streams(5)
  expect_length(streams, 5)
  expect_true(all(vapply(streams, function(s) s[1] == 10407L, logical(1))))
  expect_false(identical(streams[[1]], streams[[2]]))
  # Same seed, same streams.
  set.seed(42)
  expect_identical(roc_utils_rng_streams(5), streams)
})


test_that("generating streams leaves the caller's generator as it found it", {
  set.seed(42, kind = "Mersenne-Twister")
  before.kind <- RNGkind()
  invisible(roc_utils_rng_streams(10))
  expect_identical(RNGkind(), before.kind)
  # The caller's stream advances by the single draw used to seed the streams,
  # and stays usable afterwards.
  expect_silent(runif(1))
})
