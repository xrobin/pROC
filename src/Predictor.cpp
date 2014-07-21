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
#include <string>
#include <vector>

#include <pROC/Predictor.h>
#include <pROC/rocUtils.h>
#include "Random.h"

using Rcpp::NumericVector;
using std::string;
using std::vector;

double pROC::ResampledPredictor::operator[] (const int anIdx) const {
	return anIdx < resampledNControls ? Predictor::operator[] (controlsIdx[anIdx]) : Predictor::operator[] (casesIdx[anIdx - resampledNControls] + resampledNControls);
}

vector<int> pROC::Predictor::getOrder(const string& direction) const {
	return getPredictorOrder(*this, direction);
}


vector<int> pROC::ResampledPredictor::getOrder(const string& direction) const {
	return getPredictorOrder(*this, direction);
}

void pROC::ResampledPredictorStratified::resample() {
	setRandomIdx(nControls, controlsIdx);
	setRandomIdx(nCases, casesIdx);
	resampledNControls = controlsIdx.size();
	resampledNCases = casesIdx.size();
}

void pROC::ResampledPredictorNonStratified::resample() {
	setRandomNonStratifiedSample(nControls, nCases, controlsIdx, casesIdx);
	resampledNControls = controlsIdx.size();
	resampledNCases = casesIdx.size();
}
