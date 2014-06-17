#include <Rcpp.h>

#include <pROC/auc.h>

// [[Rcpp::export]]
double runit_partial_auc_correct(const double pauc, const double from, const double to) {
	return pROC::correctPartialAuc(pauc, from, to);
}
