/* pROC: Tools Receiver operating characteristic (ROC curves) with
   (partial) area under the curve, confidence intervals and comparison. 
   Copyright (C) 2014 Xavier Robin.

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
#include <vector>
#include <algorithm>
#include "Random.h"
#include "auc.h"
#include "RcppConversions.h"
#include "ROC.h"

using std::vector;
using std::pair;
using Rcpp::NumericVector;


// [[Rcpp::export]]
std::vector<double> bootstrapAucStratified(const size_t bootN, const SEXP& aROCSEXP, const SEXP& aucParamsSEXP) {
	// Get proper ROC
	ROC<Predictor> aROC = as<ROC<Predictor>>(aROCSEXP);
	// Get proper AUC params
	AucParams aucParams = as <AucParams>(aucParamsSEXP);

	// keep all AUCs in a vector of size bootN
	vector<double> aucs;
	aucs.reserve(bootN);
	
	ROC<ResampledPredictorStratified> aResampledROC = aROC.getResampled<ResampledPredictorStratified>();
	
	for (size_t i = 0; i < bootN; i++) {
		// Resample
		if (i > 0) { // A fresh ROC<ResampledPredictor> is already resampled, no need to resample again at the first step
			aResampledROC.resample();
		}
		// Compute AUC
		aucs.push_back(aResampledROC.auc(aucParams));
	}
	
	return aucs;
}


// [[Rcpp::export]]
std::vector<double> bootstrapAucNonStratified(const size_t bootN, const SEXP& aROCSEXP, const SEXP& aucParamsSEXP) {
	// Get proper ROC
	ROC<Predictor> aROC = as<ROC<Predictor>>(aROCSEXP);
	// Get proper AUC params
	AucParams aucParams = as <AucParams>(aucParamsSEXP);

	// keep all AUCs in a vector of size bootN
	vector<double> aucs;
	aucs.reserve(bootN);
	
	ROC<ResampledPredictorNonStratified> aResampledROC = aROC.getResampled<ResampledPredictorNonStratified>();
	
	for (size_t i = 0; i < bootN; i++) {
		// Resample
		if (i > 0) { // A fresh ROC<ResampledPredictor> is already resampled, no need to resample again at the first step
			aResampledROC.resample();
		}
		// Compute AUC
		aucs.push_back(aResampledROC.auc(aucParams));
	}
	
	return aucs;
}

// // [[Rcpp::export]]
//std::vector<double> bootstrapPaucStratified(const size_t bootN, std::vector<double> controls, std::vector<double> cases) {
