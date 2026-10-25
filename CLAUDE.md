# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

pROC is an R package (on CRAN) for displaying and analyzing ROC curves: (partial) AUC, confidence intervals, curve comparison tests, smoothing, and power/sample size. It has a small Rcpp component. It is maintained by Xavier Robin.

## Commands

Run from an R prompt in the repo root (devtools):

```r
devtools::test()                              # run tests (slow tests skipped)
run_slow_tests <- TRUE; devtools::test()      # include slow tests
devtools::test(filter = "coords")             # single file: tests/testthat/test-coords.R
testthat::snapshot_review()                   # review vdiffr / snapshot changes
devtools::check()
```

Full check from the shell (run from the **parent** directory of the repo):

```sh
VERSION=$(grep Version pROC/DESCRIPTION | sed "s/.\+ //")
R CMD build pROC
NOT_CRAN=1 RUN_SLOW_TESTS=true R CMD check pROC_$VERSION.tar.gz   # NOT_CRAN=1 enables vdiffr tests
R CMD check --as-cran pROC_$VERSION.tar.gz
```

After changing `src/*.cpp` exports, regenerate bindings with `Rcpp::compileAttributes()` (updates `R/RcppExports.R` and `src/RcppExports.cpp`; don't edit those by hand).

## Repository conventions

- **No roxygen.** `NAMESPACE` and `man/*.Rd` are maintained by hand. A new exported function or S3 method must be added to `NAMESPACE` manually (`export(...)` / `S3method(...)`), and documentation is edited directly in the `.Rd` files.
- User-visible changes go in `NEWS` (newest version at top; see "Docs and NEWS" below). Release steps (bump `Version`/`Date` in `DESCRIPTION` and `NEWS`, revdep check via `run_revdep.R`, tag `v$VERSION`) are listed in README.md.
- Branches: development happens on `develop` and is merged into `master` for release.
- Code was reformatted with `styler` (tidyverse style); keep that style.
- Deprecated arguments (e.g. `algorithm`, `progress`, `parallel`, `coords(transpose=...)`) are still accepted and emit a warning instead of erroring. Keep that backward compatibility.

## Architecture

**S3 classes and dispatch.** The main generics (`roc`, `auc`, `ci`, `ci.auc`, `ci.se`, `ci.sp`, `ci.thresholds`, `ci.coords`, `coords`, `smooth`, `var`, `cov`, `roc.test`, `plot`, `print`, ...) dispatch on these classes: `roc`, `smooth.roc`, `multiclass.roc`, `mv.multiclass.roc`, plus `formula`/`default` (and `data.frame` for `roc`) entry points. The `formula` and `default` methods usually just build a `roc` object (via `roc_utils_extract_formula` for formulas) and then call the `.roc` method. Result objects carry their own classes (`auc`, `ci.auc`, `ci.se`, `ci.sp`, `ci.thresholds`, `ci.coords`), and `R/groupGeneric.R` defines `Ops`/`Math` on them so that arithmetic drops the attributes.

**The `roc` object.** `roc.default` (R/roc.R) accepts one of three inputs: response+predictor, cases+controls, or density.cases+density.controls (the last one goes straight to a `smooth.roc`). It validates the input, resolves `levels` and `direction` (`"auto"` compares medians), then calls `roc_cc_nochecks`, the internal builder without argument checks. The object stores `sensitivities`, `specificities`, `thresholds`, `cases`, `controls`, `direction`, `percent`, `original.predictor/response`, `call`, and optionally `auc`/`ci`. A `smooth.roc` keeps the source curve in `attr(x, "roc")`, and an `auc` keeps its curve in `attr(auc, "roc")`.

**Argument pass-through via `...`.** `roc(..., auc=TRUE, ci=TRUE, plot=TRUE, smooth=TRUE)` forwards `...` to `auc`, `ci`, `plot.roc` and `smooth`. Partial AUC settings (`partial.auc`, `partial.auc.focus`, `partial.auc.correct`) are stored as attributes on the `auc` object, and later functions reuse them (`reuse.auc=TRUE`). Many functions re-call `auc(roc, ...)`, and every replay of a smoothing (roc.test, cov, are.paired, bootstrap) goes through `roc_utils_resmooth()`, which passes the stored `smoothing.args` as values. Changing argument names can therefore break callers far away.

**Percent handling.** `percent=TRUE` stores se/sp/AUC on a 0–100 scale. Statistical code (`ci.auc`, `var`, `cov`, `power.roc.test`, ...) converts to fractions with `roc_utils_unpercent()`, computes, then rescales the result. Converters for every result class are in `R/roc.utils.percent.R`. Partial AUC bounds must use the same scale as the curve.

**Statistical engines.**
- DeLong (`R/delong.R` + `src/delong.cpp::delongPlacementsCpp`): fast placement-value computation, used for AUC CI, var/cov, and `roc.test` (paired and unpaired).
- Bootstrap (`R/bootstrap.R`): every replicate goes through `bootstrap.replicate()` (one curve) or `bootstrap.test.replicate()` (two curves, `bootstrap.test`/`bootstrap.cov`). Both resample with `roc_utils_resample()`, return the caller's NA when a class is lost or the smoothing replay (`roc_utils_resmooth()`) fails, and `bootstrap.replicates()` checks each replicate's shape with `vapply()`. The statistics (`bootstrap.auc`, `bootstrap.se/sp`, `bootstrap.coords`, `bootstrap.thresholds`) are plain functions of the curve. Partial AUC and smoothed curves fall back to bootstrap because DeLong does not support them.
- Other tests: Venkatraman (`R/venkatraman.R`), Obuchowski (`R/obuchowski.R`, used by `power.roc.test`), Hand & Till multiclass AUC (`R/multiclass.R`).
- `R/coords.R` plus helpers in `R/roc.utils.R` (`roc_utils_calc_coords`, `roc_utils_optim_crit`, `roc_utils_match_coords_*_args`) compute curve coordinates and the "best" threshold (Youden or closest.topleft, with weights).

**Plotting.** Base graphics: `plot.roc`, `lines.roc`, `plot.ci` (R/plot.roc.R, R/plot.ci.R). ggplot2 (Suggests only): `ggroc` and `geom_roc_auc`/`geom_polygon_auc`. Optional packages are loaded with `load.suggested.package()`.

## Tests

testthat edition with helpers in `tests/testthat/helper-*.R`:
- `helper-rocs.R` builds shared fixture curves from the bundled `aSAH` dataset (`r.s100b`, `r.wfns.percent`, `r.ndka.partial1`, ...). Reuse these instead of creating new ones.
- `helper-*-expected.R` hold large hard-coded expected values. `print_output/` holds expected `print()` output. `_snaps/` holds vdiffr SVG snapshots.
- Wrap expensive tests in `skip_slow()`, which runs only if `RUN_SLOW_TESTS=true` or `run_slow_tests <- TRUE` in the global env. Use `expect_doppelganger` / `expect_ggroc_doppelganger` from `helper-vdiffr.R`; they skip when vdiffr/ggplot2 are missing.
- CI (GitHub Actions `check-standard.yaml`, AppVeyor on Windows, codecov) runs with `RUN_SLOW_TESTS=true`.
- Tests that use a Suggests package (ggplot2, vdiffr, MASS, logcondens, ...) need `skip_if_not_installed()`. The `depends-only` CI job runs without Suggests.
- Test hygiene:
  - no mocking of foreign classes or of `interactive()`;
  - no `recordPlot` internals;
  - restore any global state a test changes (e.g. `RNGkind`, via `on.exit`);
  - never use a random seed such as `sample.int()`; for seeds, see below.

### Bootstrap tests

Bootstrap tests (`ci.*`, `var`/`cov` and `roc.test` with `method = "bootstrap"`, ...) are random, and CI runs them on ~8 jobs per push. A test that fails 1% of the time will break CI regularly.
- **Assert only what holds for every resample.** Checks like "errors", "does not error" or "returns one best point" must hold by construction, not on typical resamples. A small margin in the fixture is not enough: aSAH points that tie or nearly tie flip between resamples. If the shared fixtures can't guarantee the property, build small synthetic data with a clear margin.
- **Ordered predictors with few levels (`r.wfns`) lose levels in resamples.** An empty level gives a threshold with the same se/sp as its neighbour, which creates ties for "best".
- **Measure before committing.** Run each new bootstrap expectation about 1000 times without a seed (e.g. `mean(replicate(1000, inherits(try(...), "try-error")))`) and check that the failure rate is exactly 0.
- **`set.seed(42)` is for exact expected values** (`print()` output, snapshots), not for hiding a flaky assertion: the next change to RNG use brings the failure back.

## Rules for code changes

These rules apply to every change: bug fixes, reviews, and the bug-hunt PR (#148). The item numbers (BUG-…) refer to PR #148, where each rule was set or applied.

### Criteria, in this order
1. **Consistency with base R and the R ecosystem.** Behaviour that is R's own semantics is not a pROC bug, e.g. partial matching, or R's missing-argument error.
2. **Least change from the last CRAN release.** Check what the release actually does, by running it.
3. **Less code:**
   - fix at the source, in one place;
   - delete point fixes rather than add machinery;
   - reuse existing helpers.

### Policies
- **P1:** degenerate input gives NA/NaN, or is dropped with a warning, as base R would do. No new `stop()` where the release returned a value.
  - A clear error replacing a cryptic one is fine.
  - A degenerate bootstrap replicate becomes NA, reported by the existing NA warning.
- **P2:** no `grepl()` on condition messages, and no blanket `suppressWarnings()`, in code or tests.
- **P3:** no floating-point tolerance as a fix; make values exact at the source, e.g. compute on integer counts.
- **P4:** fix in one layer; helpers over copies.
- **P5:** test hygiene, as in "Tests" and "Bootstrap tests" above.
- **P6:** classes from packages pROC doesn't depend on (e.g. bit64's integer64) are out of scope.
- **P7:** statistical changes need a source that was actually read, quoted with page and equation in the commit message. If a source is inaccessible, ask the maintainer for it.
- **P8:** a bug that only occurs through a shortcut, while the modular call works, is fixed only if it is a regression from the release or breaks a documented example. Otherwise document the modular call.
  - The shortcuts: `roc(smooth/ci/plot=TRUE)`, the curve-building default and formula methods of `auc`, `ci*`, `plot.roc`, `lines.roc` and `roc.test`, and `roc_`/`roc.data.frame` passing arguments on.

### Arguments and behaviour
- **Never silently reinterpret an argument.**
  - BUG-124: no `input` default on smoothed curves, so that `0.9` never silently switches between threshold and specificity.
  - BUG-231: a count `input` stays a count; the rate inputs (`tnr`, `tpr`, …) exist to hold a point fixed.
- **No silent changes:** no allow-lists or name lists that can silently miss something. Never silently drop values.
- **Twin functions behave the same.** When one gets an option, an override or a refusal, check the others:
  - `coords`/`ci.coords`;
  - `roc.test`/`cov`/`power.roc.test`;
  - `plot.roc(print.thres)`/`geom_roc_threshold`;
  - `power.roc.test(roc)`/`power.roc.test(auc =)`.
- **Users can override automatic detection consistently,** e.g. `paired` in all two-curve functions (BUG-232).
- **Respect explicit user choices:** no warning when the user explicitly asked for the behaviour (BUG-232: `paired = FALSE`).
- **Accept what the ecosystem accepts,** e.g. MASS density names in any case. Don't be stricter than the underlying package unless pROC already has the knowledge and the release already failed.
- **Statistics output follows base R,** e.g. a one-sided `conf.int` for a one-sided alternative, as `t.test` does (BUG-119).
- **Don't raise the R requirement for a test:** skip the test on old R instead (BUG-110).
- **CRAN policy:** write files only to `tempdir()`. To write anywhere else, ask first, in interactive sessions only (BUG-120).

### Docs and NEWS
- **Rd is a specification, not a changelog.** Change it only where:
  - it is wrong;
  - it misses a real, non-obvious specification;
  - a new argument needs documenting.
- **Never restate a fix in the Rd.** NEWS records fixes.
- **Don't quote error messages in the Rd,** don't state what a definition already implies, and don't repeat the same paragraph across pages.
- **NEWS:**
  - one line per user-visible change;
  - no entries for functions that are new in the release being prepared: fold them into the "Added" entry;
  - in-PR fix-ups get no entry.

### Commits
- **One commit holds code, tests, Rd and NEWS together.** Fold in-PR fix-ups into the commit that introduced the problem.
- **Full commit messages:** reprex, before/after output, and the sources.
- **When an instruction can't be applied as written, don't improvise:** skip it and report why.

### Working with the maintainer
- **Every analysis of a bug starts with "The bug":**
  - the call;
  - what you get;
  - what you should get, and why;
  - who is affected (in the release, or new in this development cycle);
  - how bad it is (silent wrong numbers, an error, cosmetic, inconsistency).

  Say plainly when it is not a bug.
- **One item at a time.** Short, direct answers to the question asked.
- **Only treat an answer as approval when it clearly answers the question asked.** Never claim the maintainer approved something without a record of it.
- **Verify by running:** the release and the current code, before reporting numbers. Check delegated work yourself.
- **Never push to `develop` or `master`, merge, or post on GitHub unless the maintainer asked for it.** Pushing to your own working branch is fine.
