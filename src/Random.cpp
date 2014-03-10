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

#include <vector>
#include <algorithm>
#include <Rcpp.h>
using std::vector;
using Rcpp::NumericVector;


void setRandomIdx(const int size, vector<int>& idxVector) {
  NumericVector idx = Rcpp::runif(size, 0, size);
  idxVector.clear();
  idxVector.insert(idxVector.begin(), idx.begin(), idx.end());
}

void setRandomNonStratifiedSample(const int controlsSize, const int casesSize,
                             vector<int>& controlsIdx, vector<int>& casesIdx) {
  // relevant sizes
  int totalSize = controlsSize + casesSize;
  // prepare return values
  controlsIdx.clear();
  casesIdx.clear();
  // take random index from R
  NumericVector idx = Rcpp::runif(totalSize, 0, totalSize);
  // sort it
  std::sort(idx.begin(), idx.end()); // avoid branch prediction fail?
  // find the first element corresponding to a case
  auto firstCaseIdx = std::lower_bound(idx.begin(), idx.end(), controlsSize);
  // Insert idx into controlsIdx and casesIdx
  controlsIdx.insert(controlsIdx.begin(), idx.begin(), firstCaseIdx);
  casesIdx.insert(casesIdx.begin(), firstCaseIdx, idx.end());
  // casesIdx has values too large by controlsSize
  for (int& caseIdxVal: casesIdx) {
    caseIdxVal -= controlsSize;
  };
}
