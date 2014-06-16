/* pROC: Tools Receiver operating characteristic (ROC curves) with
   (partial) area under the curve, confidence intervals and comparison. 
   Copyright (C) 2014 Xavier Robin

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <Rcpp.h>

#include <pROC.h>
#include <RcppConversions.h>

using namespace pROC;
using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::as;

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
  return computeAuc(aROC, aucParams);
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