#!/bin/bash
# Reverse-dependency check for pROC, end to end.
#
#   tools/revdep/run_revdep.sh              # every step
#   tools/revdep/run_revdep.sh 2 3 4        # only these steps
#
# Steps: 1 requirements | 2 library | 3 check | 4 report
#
# Must run on a node with internet access; the check downloads 200+ packages.
# Do not submit it to Slurm for that reason. Steps 2 and 3 take hours; they are
# safe to re-run, and step 2 is incremental.
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
: "${REVDEP_WORK:=/scratch/$USER/pROC_revdeps}"
: "${NCPUS:=8}"
export REVDEP_WORK NCPUS
mkdir -p "$REVDEP_WORK"
LOGDIR="$REVDEP_WORK/logs"; mkdir -p "$LOGDIR"

STEPS=("$@"); [[ ${#STEPS[@]} -eq 0 ]] && STEPS=(1 2 3 4)

# Resolve the module environment once; every step replays the snapshot.
if [[ ! -s "$REVDEP_WORK/env.snapshot" ]]; then
  echo "resolving module environment (Lmod is slow; this takes minutes) ..."
  "$HERE/env_snapshot.sh"
fi

run_step() {
  # Assign separately: bash expands every word of a `local` command before it
  # performs any of the assignments, so ${n} would still be unbound here.
  local n=$1
  local script=$2
  local log="$LOGDIR/step${n}_$(date +%Y%m%d-%H%M%S).log"
  echo "=== step $n: $script -> $log"
  "$HERE/with_env.sh" Rscript --vanilla "$HERE/$script" 2>&1 | tee "$log"
}

for s in "${STEPS[@]}"; do
  case "$s" in
    1) run_step 1 01_requirements.R ;;
    2) run_step 2 02_library.R ;;
    3) run_step 3 03_check.R ;;
    4) run_step 4 04_report.R ;;
    *) echo "unknown step: $s" >&2; exit 2 ;;
  esac
done
echo "done. work dir: $REVDEP_WORK"
