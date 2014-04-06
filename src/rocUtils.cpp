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
// #include <utility> // std::pair
#include <algorithm>
// #include <Rcpp.h>

#include "rocUtils.h"
using std::vector;
// using std::pair;
using Rcpp::NumericVector;


void makeUniqueInPlace(vector<double>& thresholds) {
  // Sort and remove consecutive duplicate values
  std::sort(thresholds.begin(), thresholds.end());
  vector<double>::iterator it = std::unique(thresholds.begin(), thresholds.end());
  // The vector still contains duplicate values the end of the vector. The new end is in 'it', use it to resize the vector.
  thresholds.resize(std::distance(thresholds.begin(), it));
}

Rcpp::NumericVector getResampledVector(const Rcpp::NumericVector& x, const std::vector<int>& idx) {
	Rcpp::NumericVector resampled(x.size());
	for (int i = 0; i < x.size(); ++i) {
		resampled[i] = x[idx[i]];
	}
	return resampled;
}
