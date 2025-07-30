.libPaths("/home/xavier/R/library-pROC_revdeps")
.libPaths("/home/xavier/R/library")

pak::pkg_install("r-lib/revdepcheck")
library(revdepcheck)

# Don't make in parallel. Avoids running out of memory on some build tasks
Sys.setenv("MAKEFLAGS"="")
revdep_reset()
revdepcheck::revdep_check(num_workers=2, timeout = as.difftime(60, units = "mins"))

