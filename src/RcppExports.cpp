// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RcppVersion
String RcppVersion();
RcppExport SEXP _pROC_RcppVersion() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(RcppVersion());
    return rcpp_result_gen;
END_RCPP
}
// delongPlacementsCpp
List delongPlacementsCpp(List roc);
RcppExport SEXP _pROC_delongPlacementsCpp(SEXP rocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type roc(rocSEXP);
    rcpp_result_gen = Rcpp::wrap(delongPlacementsCpp(roc));
    return rcpp_result_gen;
END_RCPP
}
// rocUtilsPerfsAllC
List rocUtilsPerfsAllC(NumericVector thresholds, NumericVector controls, NumericVector cases, std::string direction);
RcppExport SEXP _pROC_rocUtilsPerfsAllC(SEXP thresholdsSEXP, SEXP controlsSEXP, SEXP casesSEXP, SEXP directionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type controls(controlsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cases(casesSEXP);
    Rcpp::traits::input_parameter< std::string >::type direction(directionSEXP);
    rcpp_result_gen = Rcpp::wrap(rocUtilsPerfsAllC(thresholds, controls, cases, direction));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pROC_RcppVersion", (DL_FUNC) &_pROC_RcppVersion, 0},
    {"_pROC_delongPlacementsCpp", (DL_FUNC) &_pROC_delongPlacementsCpp, 1},
    {"_pROC_rocUtilsPerfsAllC", (DL_FUNC) &_pROC_rocUtilsPerfsAllC, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_pROC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
