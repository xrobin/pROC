#include <Rcpp.h>
using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::as;

#include "ROC.h"
#include "RcppConversions.h"

// [[Rcpp::export]]
List computeSeSpList(const NumericVector& thresholds, const NumericVector& controls, const NumericVector& cases, const std::string& direction) {
	
	ROC<Predictor> aROC(controls, cases, direction);

	List ret;
	ret["se"] = aROC.getSensitivity();
	ret["sp"] = aROC.getSpecificity();
	
	return ret;
}

// [[Rcpp::export]]
double computeAuc(const Rcpp::List& aROCList, const Rcpp::List& aucParamsList) {
  AucParams aucParams = as<AucParams>(aucParamsList);
  ROC<Predictor> aROC = as<ROC<Predictor>>(aROCList);
  return aROC.auc(aucParams);
}

// [[Rcpp::export]]
double computeAucSeSp(const std::vector<double>& se, const std::vector<double>& sp, const Rcpp::List& aucParamsList) {
  AucParams aucParams = as<AucParams>(aucParamsList);
  if (aucParams.partial) {
    return computePartialAuc(se, sp, aucParams);
  }
  else {
    return computeFullAuc(se, sp);
  }
}