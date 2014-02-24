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
#include <vector>
#include "Predictor.h"
using namespace Rcpp;
using std::vector;

// [[Rcpp::export]]
List computeSeSpList(NumericVector thresholds, NumericVector controls, NumericVector cases, std::string direction) {
  // Vector sizes
  const size_t ncontrols(controls.size());
  const size_t ncases(cases.size());
  const size_t npredictors = ncontrols + ncases;
  const size_t nthresholds = thresholds.size();
  
  const Predictor predictor(controls, cases);
  
  // Compute an index vector for controls and cases as R's order() function
  vector<size_t> index = predictor.getOrder(direction);
  
  // Cummulative sum
  // no need for tp/fp, compute se/sp directly during the cummulative sum
  vector<double> se(npredictors), sp(npredictors);
  // And store the cummulative sums (tp, fp) in two variables
  size_t currentTpSum = 0, currentFpSum = 0;

  size_t validPositions = 0;
  for (size_t i = 0; i < npredictors; ++i) {
    // Compute Se/Sp
    if (predictor.isControl(index[i])) { // we have one control
      ++currentFpSum;
    }
    else { // we have one case
      ++currentTpSum;
    }
    
    // Determine if if is a duplicate
    bool currentDupPred = false;
    // Is predictor[i] the same as predictor[i+1]?
    if (i < npredictors - 1) {
      currentDupPred = predictor[index[i + 1]] == predictor[index[i]];
    }
    // If different, add the Se/Sp as a valid position
    if (!currentDupPred) {
      se[validPositions] = static_cast<double>(currentTpSum) / ncases;
      sp[validPositions] = (ncontrols - static_cast<double>(currentFpSum)) / ncontrols;
      ++validPositions;
    }
  }
  
  // Se/sp were over-allocated - resize to the valid positions only
  se.resize(validPositions + 1);
  sp.resize(validPositions + 1);
  
  // Reverse - can we find a way to do it right from the start?
  std::reverse(se.begin(), se.end() - 1);
  std::reverse(sp.begin(), sp.end() - 1);
  
  // Anchor the last position to 1/0
  se[validPositions] = 0;
  sp[validPositions] = 1;

  // Ensure the sum of !duplicated == thresholds.size()
  // Todo: maybe remove this in a future version of pROC...
  if (validPositions != nthresholds - 1) {
    stop("Bug in pROC: fast algorithm (C++ version) computed an incorrect number of sensitivities and specificities.");
  }

  List ret;
  ret["se"] = se;
  ret["sp"] = sp;
  
  return ret;
}
