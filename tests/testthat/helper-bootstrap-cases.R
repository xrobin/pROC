# The bootstrap code paths exercised by test-bootstrap.R, and the expectations
# in helper-bootstrap-expected.R are generated from this same list by
# tools/gen_bootstrap_expected.R. Keeping one list means the tests and the
# generated expectations can never drift apart.
#
# boot.n is deliberately tiny: these are characterisation tests, pinning the
# RNG draw order, the aggregation and the output shape, not the statistical
# properties of the bootstrap. They must stay cheap enough to run on CRAN.

bootstrap.B <- 10

# Compact, structure-aware signature of a bootstrap result.
bootstrap.sig <- function(x) {
  if (inherits(x, "htest")) {
    c(unname(x$statistic), x$p.value)
  } else if (is.list(x) && !is.data.frame(x) && !is.matrix(x)) {
    as.numeric(unlist(x))
  } else {
    as.numeric(x)
  }
}

bootstrap.shape <- function(x) {
  d <- dim(x)
  if (is.null(d)) length(x) else d
}

bootstrap.cases <- list(
  ## ---- ci.auc ---------------------------------------------------------
  ci.auc.strat       = quote(ci.auc(r.wfns, method = "bootstrap", boot.n = bootstrap.B)),
  ci.auc.nonstrat    = quote(ci.auc(r.wfns, method = "bootstrap", boot.n = bootstrap.B, boot.stratified = FALSE)),
  ci.auc.percent     = quote(ci.auc(r.wfns.percent, method = "bootstrap", boot.n = bootstrap.B)),
  ci.auc.continuous  = quote(ci.auc(r.s100b, method = "bootstrap", boot.n = bootstrap.B)),
  ci.auc.partial     = quote(ci.auc(r.s100b.partial, method = "bootstrap", boot.n = bootstrap.B)),
  ci.auc.conflevel   = quote(ci.auc(r.s100b, method = "bootstrap", boot.n = bootstrap.B, conf.level = 0.90)),
  ci.auc.smooth      = quote(ci.auc(smooth(r.ndka), boot.n = bootstrap.B)),
  ci.auc.smooth.ns   = quote(ci.auc(smooth(r.ndka), boot.n = bootstrap.B, boot.stratified = FALSE)),

  ## ---- ci.se / ci.sp --------------------------------------------------
  ci.se.strat        = quote(ci.se(r.s100b, specificities = c(0.1, 0.5, 0.9), boot.n = bootstrap.B)),
  ci.se.nonstrat     = quote(ci.se(r.s100b, specificities = c(0.1, 0.5, 0.9), boot.n = bootstrap.B, boot.stratified = FALSE)),
  ci.se.percent      = quote(ci.se(r.s100b.percent, specificities = c(10, 50, 90), boot.n = bootstrap.B)),
  ci.se.smooth       = quote(ci.se(smooth(r.ndka), specificities = c(0.1, 0.5, 0.9), boot.n = bootstrap.B)),
  ci.sp.strat        = quote(ci.sp(r.s100b, sensitivities = c(0.1, 0.5, 0.9), boot.n = bootstrap.B)),
  ci.sp.nonstrat     = quote(ci.sp(r.s100b, sensitivities = c(0.1, 0.5, 0.9), boot.n = bootstrap.B, boot.stratified = FALSE)),
  ci.sp.smooth       = quote(ci.sp(smooth(r.ndka), sensitivities = c(0.1, 0.5, 0.9), boot.n = bootstrap.B)),

  ## ---- ci.thresholds --------------------------------------------------
  ci.thr.strat       = quote(ci.thresholds(r.s100b, thresholds = c(0.1, 0.3, 0.5), boot.n = bootstrap.B)),
  ci.thr.nonstrat    = quote(ci.thresholds(r.s100b, thresholds = c(0.1, 0.3, 0.5), boot.n = bootstrap.B, boot.stratified = FALSE)),

  ## ---- ci.coords ------------------------------------------------------
  ci.coords.sp       = quote(ci.coords(r.s100b, x = c(0.1, 0.5, 0.9), input = "specificity",
                                       ret = c("sensitivity", "ppv"), boot.n = bootstrap.B)),
  ci.coords.nonstrat = quote(ci.coords(r.s100b, x = c(0.1, 0.5, 0.9), input = "specificity",
                                       ret = c("sensitivity", "ppv"), boot.n = bootstrap.B,
                                       boot.stratified = FALSE)),
  ci.coords.best     = quote(ci.coords(r.s100b, x = "best", ret = c("sensitivity", "specificity"),
                                       boot.n = bootstrap.B)),

  ## ---- var / cov ------------------------------------------------------
  var.boot           = quote(var(r.s100b, method = "bootstrap", boot.n = bootstrap.B)),
  var.boot.ns        = quote(var(r.s100b, method = "bootstrap", boot.n = bootstrap.B, boot.stratified = FALSE)),
  var.boot.smooth    = quote(var(smooth(r.ndka), boot.n = bootstrap.B)),
  cov.boot           = quote(cov(r.s100b, r.ndka, method = "bootstrap", boot.n = bootstrap.B)),
  cov.boot.ns        = quote(cov(r.s100b, r.ndka, method = "bootstrap", boot.n = bootstrap.B, boot.stratified = FALSE)),

  ## ---- roc.test -------------------------------------------------------
  roctest.boot       = quote(roc.test(r.s100b, r.ndka, method = "bootstrap", boot.n = bootstrap.B)),
  roctest.boot.ns    = quote(roc.test(r.s100b, r.ndka, method = "bootstrap", boot.n = bootstrap.B, boot.stratified = FALSE)),
  roctest.venk       = quote(roc.test(r.s100b, r.ndka, method = "venkatraman", boot.n = bootstrap.B)),

  ## ---- paths of the replicate wrapper not pinned above -----------------
  ## (cluster F decision on PR #148, generated before the refactor)
  roctest.smooth     = quote(roc.test(smooth(r.s100b), smooth(r.ndka, method = "density"), boot.n = bootstrap.B)),
  roctest.smooth.ns  = quote(roc.test(smooth(r.s100b), smooth(r.ndka), boot.n = bootstrap.B, boot.stratified = FALSE)),
  cov.smooth         = quote(cov(smooth(r.s100b), smooth(r.ndka), boot.n = bootstrap.B)),
  roctest.unpaired   = quote(roc.test(r.s100b, roc(aSAH$outcome[1:70], aSAH$wfns[1:70], quiet = TRUE),
                                      method = "bootstrap", paired = FALSE, boot.n = bootstrap.B)),
  roctest.unpaired.ns = quote(roc.test(r.s100b, roc(aSAH$outcome[1:70], aSAH$wfns[1:70], quiet = TRUE),
                                       method = "bootstrap", paired = FALSE, boot.n = bootstrap.B,
                                       boot.stratified = FALSE)),
  # Both curves lose a class in some resamples (3 cases in 30, 2 in 10)
  roctest.lost.class = quote(roc.test(
    roc(c(rep(0, 27), rep(1, 3)), c(1:27, 20.5, 25.5, 28), quiet = TRUE),
    roc(c(rep(0, 8), rep(1, 2)), c(1:8, 4.5, 9), quiet = TRUE),
    method = "bootstrap", paired = FALSE, boot.stratified = FALSE, boot.n = 50
  )),
  roctest.sp         = quote(roc.test(r.s100b, r.ndka, method = "specificity", specificity = 0.9, boot.n = bootstrap.B)),
  roctest.pauc.percent = quote(roc.test(
    roc(aSAH$outcome, aSAH$s100b, percent = TRUE, partial.auc = c(100, 80), partial.auc.correct = TRUE, quiet = TRUE),
    roc(aSAH$outcome, aSAH$ndka, percent = TRUE, partial.auc = c(100, 80), partial.auc.correct = TRUE, quiet = TRUE),
    method = "bootstrap", boot.n = bootstrap.B
  )),
  ci.coords.smooth   = quote(ci.coords(smooth(r.ndka), x = c(0.5, 0.9), input = "specificity",
                                       ret = c("sensitivity", "specificity"), boot.n = bootstrap.B)),
  ci.coords.random   = quote(ci.coords(r.s100b, x = "best", best.policy = "random",
                                       ret = c("sensitivity", "specificity"), boot.n = bootstrap.B))
)

# Run one case reproducibly. Warnings are expected on some paths (NA replicates,
# ties) and are not what these tests are about.
bootstrap.run <- function(name, seed = 42) {
  set.seed(seed)
  withCallingHandlers(eval(bootstrap.cases[[name]], envir = parent.frame()),
                      warning = function(w) invokeRestart("muffleWarning"))
}
