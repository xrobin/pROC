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
#include "Predictor.h"
#include "rocUtils.h"
#include "Random.h"

using std::vector;
using std::string;
using Rcpp::NumericVector;

vector<int> Predictor::getOrder(const std::string& direction) const {
	return getPredictorOrder(*this, direction);
}


vector<int> ResampledPredictor::getOrder(const std::string& direction) const {
	return getPredictorOrder(*this, direction);
}

void ResampledPredictorStratified::resample() {
	setRandomIdx(nControls, controlsIdx);
	setRandomIdx(nCases, casesIdx);
}

void ResampledPredictorNonStratified::resample() {
	setRandomNonStratifiedSample(nControls, nCases, controlsIdx, casesIdx);
}
