#!/bin/bash
# Run a command in the revdep environment, replaying the snapshot taken by
# env_snapshot.sh. Lmod costs ~9 minutes for these modules; the snapshot is
# instant.
#
# env -i drops everything not in the snapshot, so any variable meant to steer a
# run must be forwarded explicitly. Forget one and it silently has no effect:
# REVDEP_ONLY was dropped that way, and a two-package smoke test quietly became
# a 206-package run. Hence the prefix expansion rather than a hand-kept list.
#   usage: with_env.sh <command> [args...]
set -euo pipefail
: "${REVDEP_WORK:=/scratch/$USER/pROC_revdeps}"
SNAP="$REVDEP_WORK/env.snapshot"
[[ -s "$SNAP" ]] || { echo "no snapshot; run env_snapshot.sh first" >&2; exit 1; }
mapfile -d '' vars < "$SNAP"
over=()
set +u
for v in ${!REVDEP_@} NCPUS USER HOME; do
  [[ -n "${!v:-}" ]] && over+=("$v=${!v}")
done
set -u
exec env -i "${vars[@]}" "${over[@]}" "$@"
