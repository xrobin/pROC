library(pROC)
data(aSAH)

context("ci.formula")

test_that("bootstrap cov works with smooth and !reuse.auc", {
  skip_slow()
  if (getRversion() > "3.6.0") {
    # Restore it: leaving sample.kind set leaks into every test file that runs
    # after this one, making their results depend on the order tests ran in.
    previous.kind <- RNGkind()
    on.exit(suppressWarnings(do.call(RNGkind, as.list(previous.kind))), add = TRUE)
    suppressWarnings(RNGkind(sample.kind = "Rounding"))
  }

  for (pair in list(
    list(ci, list()),
    list(ci.se, list(boot.n = 10)),
    list(ci.sp, list(boot.n = 10)),
    list(ci.thresholds, list(boot.n = 10)),
    list(ci.coords, list(boot.n = 10, x = 0.5)),
    list(ci.auc, list())
  )) {
    fun <- pair[[1]]

    # First calculate ci with .default
    args.default <- c(
      list(
        response = aSAH$outcome,
        predictor = aSAH$s100b
      ),
      pair[[2]]
    )
    set.seed(42) # For reproducible CI
    obs.default <- do.call(fun, args.default)

    # Then with .formula
    args.formula <- c(
      list(
        formula = outcome ~ s100b,
        data = aSAH
      ),
      pair[[2]]
    )
    set.seed(42) # For reproducible CI
    obs.formula <- do.call(fun, args.formula)

    # Here we check both returned the same result
    # We ignore attributes, as we have different
    # roc objects, and unfortunately equivalent means
    # we only test near equality
    expect_equivalent(obs.default, obs.formula)
  }
})

test_that("ci formula and default methods work with smooth = TRUE", {
  s <- smooth(roc(aSAH$outcome, aSAH$s100b, quiet = TRUE))
  for (pair in list(
    list(ci, list(), "ci.auc"),
    list(ci, list(of = "se", specificities = 0.5), "ci.se"),
    list(ci, list(of = "coords", x = 0.5), "ci.coords"),
    list(ci.auc, list(), "ci.auc"),
    list(ci.se, list(specificities = 0.5), "ci.se"),
    list(ci.sp, list(sensitivities = 0.5), "ci.sp"),
    list(ci.coords, list(x = 0.5, input = "sp"), "ci.coords")
  )) {
    fun <- pair[[1]]
    args <- c(list(smooth = TRUE, boot.n = 5, quiet = TRUE), pair[[2]])
    set.seed(42)
    obs.smooth <- do.call(fun, c(list(s), args[!names(args) %in% c("smooth", "quiet")]))
    set.seed(42)
    obs.default <- do.call(fun, c(list(response = aSAH$outcome, predictor = aSAH$s100b), args))
    set.seed(42)
    obs.formula <- do.call(fun, c(list(formula = outcome ~ s100b, data = aSAH), args))
    expect_s3_class(obs.default, pair[[3]])
    expect_s3_class(obs.formula, pair[[3]])
    expect_equal(as.numeric(unlist(unclass(obs.default))), as.numeric(unlist(unclass(obs.smooth))))
    expect_equal(as.numeric(unlist(unclass(obs.formula))), as.numeric(unlist(unclass(obs.smooth))))
  }
})

test_that("auc and ci default/formula return the NA/NaN of roc()", {
  set.seed(42)
  d <- data.frame(
    y = c(0, 0, 0, 1, 1, 1),
    xna = c(1, 2, NA, 4, 5, 6),
    xinf = c(1, 2, Inf, 4, 5, 6)
  )
  for (fun in list(auc, ci, ci.auc, ci.se, ci.sp, ci.thresholds)) {
    expect_identical(fun(d$y, d$xna, na.rm = FALSE, quiet = TRUE), NA)
    expect_identical(fun(y ~ xna, d, na.rm = FALSE, quiet = TRUE), NA)
    expect_warning(r <- fun(d$y, d$xinf, quiet = TRUE), "Infinite")
    expect_identical(r, NaN)
    expect_warning(r <- fun(y ~ xinf, d, quiet = TRUE), "Infinite")
    expect_identical(r, NaN)
  }
  expect_identical(ci.coords(d$y, d$xna, x = 3, na.rm = FALSE, quiet = TRUE), NA)
  expect_identical(ci.coords(y ~ xna, d, x = 3, na.rm = FALSE, quiet = TRUE), NA)
  expect_warning(r <- ci.coords(d$y, d$xinf, x = 3, quiet = TRUE), "Infinite")
  expect_identical(r, NaN)
  expect_warning(r <- ci.coords(y ~ xinf, d, x = 3, quiet = TRUE), "Infinite")
  expect_identical(r, NaN)
})
