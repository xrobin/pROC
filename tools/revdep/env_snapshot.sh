#!/bin/bash
# Lmod costs ~9 minutes for these modules on this filesystem. Resolve the
# environment once, snapshot it, and let with_env.sh replay it instantly.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/00_env.sh"
env -0 > "$REVDEP_WORK/env.snapshot"
echo "snapshot written to $REVDEP_WORK/env.snapshot"
