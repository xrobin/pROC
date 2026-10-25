# pROC: rules for code changes

pROC is maintained by Xavier Robin. These rules apply to every change: bug fixes, reviews, and the bug-hunt PR (#148). The item numbers (BUG-…) refer to PR #148, where each rule was set or applied.

## Criteria, in this order
1. **Consistency with base R and the R ecosystem.** Behaviour that is R's own semantics is not a pROC bug, e.g. partial matching, or R's missing-argument error.
2. **Least change from the last CRAN release.** Check what the release actually does, by running it.
3. **Less code:**
   - fix at the source, in one place;
   - delete point fixes rather than add machinery;
   - reuse existing helpers.

## Policies
- **P1:** degenerate input gives NA/NaN, or is dropped with a warning, as base R would do. No new `stop()` where the release returned a value.
  - A clear error replacing a cryptic one is fine.
  - A degenerate bootstrap replicate becomes NA, reported by the existing NA warning.
- **P2:** no `grepl()` on condition messages, and no blanket `suppressWarnings()`, in code or tests.
- **P3:** no floating-point tolerance as a fix; make values exact at the source, e.g. compute on integer counts.
- **P4:** fix in one layer; helpers over copies.
- **P5:** test hygiene:
  - fixed seeds;
  - no mocking of foreign classes or of `interactive()`;
  - no `recordPlot` internals;
  - restore any global state you change (e.g. `RNGkind`, via `on.exit`).
- **P6:** classes from packages pROC doesn't depend on (e.g. bit64's integer64) are out of scope.
- **P7:** statistical changes need a source that was actually read, quoted with page and equation in the commit message. If a source is inaccessible, ask the maintainer for it.
- **P8:** a bug that only occurs through a shortcut, while the modular call works, is fixed only if it is a regression from the release or breaks a documented example. Otherwise document the modular call.
  - The shortcuts: `roc(smooth/ci/plot=TRUE)`, the curve-building default and formula methods of `auc`, `ci*`, `plot.roc`, `lines.roc` and `roc.test`, and `roc_`/`roc.data.frame` passing arguments on.

## Arguments and behaviour
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

## Docs and NEWS
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

## Commits
- **One commit holds code, tests, Rd and NEWS together.** Fold in-PR fix-ups into the commit that introduced the problem.
- **Full commit messages:** reprex, before/after output, and the sources.
- **When an instruction can't be applied as written, don't improvise:** skip it and report why.

## Working with the maintainer
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
- **Never post to GitHub, push, or merge without explicit approval.**
