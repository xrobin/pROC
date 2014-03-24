#include <Rcpp.h>
using Rcpp::List;
using Rcpp::NumericVector;

#include "ROC.h"
#include "RcppConversions.h"

// [[Rcpp::export]]
List computeSeSpList(const NumericVector& thresholds, const NumericVector& controls, const NumericVector& cases, const char direction) {
	
	ROC<Predictor> aROC(controls, cases, direction);

	List ret;
	ret["se"] = aROC.getSensitivity();
	ret["sp"] = aROC.getSpecificity();
	
	return ret;
}

// [[Rcpp::export]]
double computeAuc(const std::vector<double>& se, const std::vector<double>& sp, const Rcpp::List& aucParamsList) {
  AucParams aucParams = as<AucParams>(aucParamsList);
  if (aucParams.partial) {
    return computePartialAuc(se, sp, aucParams);
  }
  else {
    return computeFullAuc(se, sp);
  }
}