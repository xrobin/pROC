# pROC reverse-dependency checks — environment.
# Source this; do not execute it. Sourcing it inside a pipeline silently
# discards every module load, because the pipeline runs it in a subshell.
#
# Verified on sciCORE worker10-beta with R 4.6.1-gfbf-2026.1 (Sep 2026).

set +u
module purge
for m in R Pandoc Xvfb UDUNITS CMake ImageMagick MPFR NLopt protobuf \
         GSL/2.8-GCC-15.2.0 Rust PostgreSQL libclc JAGS GDAL; do
  module load "$m" || echo "MODULE FAIL: $m" >&2
done
set -u
# GDAL must come after JAGS and must not be preceded by a visible PROJ:
# it depends on a hidden PROJ and the autoswap is blocked.

: "${REVDEP_WORK:=/scratch/$USER/pROC_revdeps}"
export REVDEP_WORK
export REVDEP_LIB="$REVDEP_WORK/Library"
export TMPDIR="$REVDEP_WORK/tmp"
mkdir -p "$REVDEP_LIB" "$TMPDIR"

# Exactly one user library. R_LIBS_SITE is silenced deliberately: the
# EasyBuild R-bundle-CRAN ships its own pROC, and if it is visible the
# baseline side checks against that instead of the CRAN release.
export R_LIBS_USER="$REVDEP_LIB"
export R_LIBS="$REVDEP_LIB"
export R_LIBS_SITE=/dev/null
export R_ENVIRON_USER=/dev/null
export R_PROFILE_USER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/Rprofile.R"

# R-level parallelism only; nested make on top of it overloads the node.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MAKEFLAGS=-j1

[[ -n "${EBROOTUDUNITS:-}" ]] && \
  export UDUNITS2_XML_PATH="$EBROOTUDUNITS/share/udunits/udunits2.xml"

# Headers with no GCCcore-15.2.0 module of their own (jq, poppler, rsvg,
# OpenCL). Prefer a module whenever one exists; this prefix is the fallback.
LP="$REVDEP_WORK/local/usr"
if [[ -d "$LP/include" ]]; then
  export PKG_CONFIG_PATH="$LP/lib/x86_64-linux-gnu/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
  export CPATH="$LP/include${CPATH:+:$CPATH}"
  export LIBRARY_PATH="$LP/lib/x86_64-linux-gnu${LIBRARY_PATH:+:$LIBRARY_PATH}"
  export LD_LIBRARY_PATH="$LP/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
[[ -f "$REVDEP_WORK/Makevars" ]] && export R_MAKEVARS_USER="$REVDEP_WORK/Makevars"

BOOST=/scicore/soft/easybuild/apps/Boost/1.90.0-GCCcore-15.2.0/lib
[[ -d "$BOOST" ]] && export LD_LIBRARY_PATH="$BOOST${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

TEXROOT="$REVDEP_WORK/texroot/usr/share/texlive/texmf-dist"
if [[ -d "$TEXROOT" ]]; then
  export TEXMFHOME="$TEXROOT"
  export TEXINPUTS="$TEXROOT//:${TEXINPUTS:-}"
  export TFMFONTS="$TEXROOT/fonts/tfm//:${TFMFONTS:-}"
  export VFFONTS="$TEXROOT/fonts/vf//:${VFFONTS:-}"
  export T1FONTS="$TEXROOT/fonts/type1//:${T1FONTS:-}"
  export TEXFONTMAPS="$TEXROOT/fonts/map//:${TEXFONTMAPS:-}"
  export ENCFONTS="$TEXROOT/fonts/enc//:${ENCFONTS:-}"
fi

# Keep every per-user cache inside the work directory. R_user_dir() resolves
# through the XDG variables, and $HOME is read-only on some nodes (and inside
# the agent sandbox), which makes packages fail in ways that look unrelated:
# gdtools cannot write its font cache, so ggiraph fails to load, and ggimage,
# shadowtext, ggtree, enrichplot and clusterProfiler all fail behind it.
export XDG_CACHE_HOME="$REVDEP_WORK/cache"
export XDG_DATA_HOME="$REVDEP_WORK/share"
export XDG_CONFIG_HOME="$REVDEP_WORK/config"
export GDTOOLS_CACHE_DIR="$REVDEP_WORK/cache/gdtools"
mkdir -p "$XDG_CACHE_HOME" "$XDG_DATA_HOME" "$XDG_CONFIG_HOME" "$GDTOOLS_CACHE_DIR"

# Xvfb ships only misc-fixed; knitr plots at 120dpi ask for Helvetica 14.
[[ -d /usr/share/fonts/X11/Type1 ]] && \
  export XVFB_SERVERARGS="-screen 0 1280x1024x24 -fp /usr/share/fonts/X11/Type1"
