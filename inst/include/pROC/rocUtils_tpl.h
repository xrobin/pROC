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
#pragma once 

#include <algorithm> // std::reverse_copy
#include <vector>

namespace pROC {
	template<typename T> 
	std::vector<T> getReversedVector(const std::vector<T>& x) {
		std::vector<T> y;
		y.assign(x.rbegin(), x.rend());
		return y;
	}

	template <typename PredictorType>
	std::vector<double> computeThresholds(const PredictorType& predictor, const std::string& aDirection) {
		std::vector<double> thresholds;
		Rcpp::NumericVector cases = predictor.getCases();
		Rcpp::NumericVector controls = predictor.getControls();
		thresholds.insert(thresholds.begin(), cases.begin(), cases.end());
		thresholds.insert(thresholds.begin(), controls.begin(), controls.begin());
		makeUniqueInPlace(thresholds);
		if (aDirection == ">") { // reverse when ">"
			std::reverse(thresholds.begin(), thresholds.end());
		}
		return thresholds;
	}
}