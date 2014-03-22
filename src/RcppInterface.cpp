#include <Rcpp.h>
using Rcpp::List;
using Rcpp::NumericVector;

#include "ROC.h"

// [[Rcpp::export]]
List computeSeSpList(const NumericVector& thresholds, const NumericVector& controls, const NumericVector& cases, const char direction) {
	
	ROC<Predictor> aROC(controls, cases, direction);

	List ret;
	ret["se"] = aROC.getSensitivity();
	ret["sp"] = aROC.getSpecificity();
	
	return ret;
}