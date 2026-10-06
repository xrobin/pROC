# pROC reverse-dependency checks

Checks every CRAN package that depends on pROC twice — once against the CRAN
release, once against this working tree — and reports the ones that get worse.

## Running it

```sh
tools/revdep/run_revdep.sh          # all four steps
tools/revdep/run_revdep.sh 3 4      # re-check and re-report only
```

Run it on a node with internet access and several cores. **Do not submit it to
Slurm**: the check downloads a few hundred packages as it goes. Steps 2 and 3
take hours. Everything lands in `$REVDEP_WORK`, by default
`/scratch/$USER/pROC_revdeps`.

| variable | default | meaning |
|---|---|---|
| `REVDEP_WORK` | `/scratch/$USER/pROC_revdeps` | work directory |
| `NCPUS` | 8 | parallel installs/checks (`make` stays at `-j1`) |
| `REVDEP_BIOC` | 3.23 | Bioconductor release to match the R version |
| `REVDEP_SIDES` | `old,new` | run one side only, e.g. `new` |
| `REVDEP_TIMEOUT` | 45m | per-package check timeout |
| `REVDEP_ONLY` | unset | check only these packages, comma-separated |

Smoke-test the whole pipeline on two packages before committing to a run of
several hours:

```sh
REVDEP_ONLY=aplore3,alternativeROC tools/revdep/run_revdep.sh 3 4
```

## The steps

**1 — `01_requirements.R`** works out what to check and what to install. The
set of reverse dependencies is recomputed from live repository metadata every
run, never reused: packages are added to and archived from CRAN continually.
It also separates dependencies that are on CRAN/Bioconductor, in the CRAN
Archive, in a repository declared via `Additional_repositories`, or nowhere.

**2 — `02_library.R`** brings the shared library up to that specification:
installs what is missing, refreshes what is outdated, and — importantly —
**load-tests every package in a fresh process**. A package can be present in a
library yet fail to load, and `R CMD check` then blames the package that
depends on it. Anything unloadable is reinstalled. The pass repeats until it
stops improving.

**3 — `03_check.R`** builds both sides in a fresh timestamped run directory and
checks all reverse dependencies against each.

**4 — `04_report.R`** compares the two sides and writes `revdep_report.md` and
`revdep_summary.csv` into the run directory. The report lists every package
whose status got worse, with the specific check items that degraded.

## Things that will silently corrupt a run

These are not hypothetical; each one invalidated an earlier attempt.

**`reverse=` must be a character vector.** `tools::check_packages_in_dir()`
accepts either an explicit vector of package names or a `list(which=...)`
specification. With the list form it computes the reverse dependencies of
*every tarball in the directory*. That is fine on the first pass, when the
directory holds only pROC — but on a restart the reverse dependencies are
already sitting there, and 206 packages becomes 2138. Step 3 always passes the
explicit vector.

**pROC must be visible only on the side under test.** Two other copies exist on
this system: the EasyBuild `R-bundle-CRAN` ships one, and the shared library
holds one after step 2. If either is on the library path, the baseline side
quietly checks against the wrong version and the whole comparison is
meaningless. `00_env.sh` sets `R_LIBS_SITE=/dev/null` for the first, and step 3
moves pROC out of the shared library for the second. Each side then installs
its own from its own tarball, because a tarball in the check directory
outranks the repository copy.

**`_R_CHECK_FORCE_SUGGESTS_` must be false.** Under `--as-cran` it defaults to
true, and a suggested package that is not on CRAN then aborts the check with an
ERROR at "checking package dependencies" — before a single test runs. Nine of
the reverse dependencies suggest a package that is archived, GitHub-only or
commercial. CRAN's own machines set this to false for the same reason.

**Do not disable `\donttest` examples.** It is tempting, because a few Shiny and
OpenCPU examples never return. But `\donttest` is often exactly where a package
calls `roc()`, so switching it off hides the regressions worth finding. Step 3
keeps them and bounds them with `_R_CHECK_*_ELAPSED_TIMEOUT_` instead.

**Never source `00_env.sh` inside a pipeline.** `source 00_env.sh | grep ...`
runs it in a subshell and every module load is discarded, leaving no R on the
`PATH` and a confusing cascade of failures.

## System libraries

Get them from EasyBuild modules — `00_env.sh` loads fourteen, which between
them covered every one of the ~850 source installs in the September 2026 run
without a single compile failure. If something new fails to build, look for a
module before anything else; `module -t avail` lists about 4000. The
`$REVDEP_WORK/local/usr` prefix of headers extracted from Ubuntu packages is a
legacy fallback for the few libraries with no GCCcore-15.2.0 module (jq,
poppler, rsvg, OpenCL). One of its pkg-config files (`libjq.pc`) was missing
the mandatory `Description:` field, which made pkg-config reject it silently;
that is the failure mode to expect from that directory.
