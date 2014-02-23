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
List rocUtilsPerfsAllC(NumericVector thresholds, NumericVector controls, NumericVector cases, std::string direction) {
  NumericVector se(thresholds.size());
  NumericVector sp(thresholds.size());
  long tp, tn;
  long i; // iterator over cases & controls
  
  if (direction == ">") {
    for (long t = 0; t < thresholds.size(); t++) {
      double threshold = thresholds(t);
        tp = 0;
        for (i = 0; i < cases.size(); i++) {
          if (cases(i) <= threshold) {
            tp++;
        }
      }
      se(t) = (double)tp / cases.size();
      tn = 0;
      for (i = 0; i < controls.size(); i++) {
        if (controls(i) > threshold) {
          tn++;
        }
      }
      sp(t) = (double)tn / controls.size();
    }
  }
 else {
    for (long t = 0; t < thresholds.size(); t++) {
      double threshold = thresholds(t);
      tp = 0;
      for (i = 0; i < cases.size(); i++) {
        if (cases(i) >= threshold) {
          tp++;
        }
      }
      se(t) = (double)tp / cases.size();
      long tn = 0;
      for (i = 0; i < controls.size(); i++) {
        if (controls(i) < threshold) {
          tn++;
        }
      }
      sp(t) = (double)tn / controls.size();
    }
  }
  List ret;
  ret["se"] = se;
  ret["sp"] = sp;
  return(ret);
}


// [[Rcpp::export]]
List rocUtilsPerfsCumsumC(NumericVector thresholds, NumericVector controls, NumericVector cases, std::string direction) {
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
  // Assess the duplication stage in the same loop
  vector<bool> duplicated(npredictors); 
  
  for (size_t i = 0; i < npredictors; ++i) {
    if (index[i] < ncontrols) { // we have one control
      sp[i] = (ncontrols - static_cast<double>(++currentFpSum)) / ncontrols;
      se[i] = static_cast<double>(currentTpSum) / ncases;
    }
    else { // we have one case
      sp[i] = (ncontrols - static_cast<double>(currentFpSum)) / ncontrols;
      se[i] = static_cast<double>(++currentTpSum) / ncases;
    }
    // We can compute the duplicates in the same loop
    bool currentDupPred = false;
    bool currentDupSeSp = false;
    // Is predictor[i] the same as predictor[i+1]?
    if (i < npredictors - 1) {
      size_t nextIdx = index[i + 1];
      size_t currIdx = index[i];
      currentDupPred = predictor[nextIdx] == predictor[currIdx];
    }
    // Are SE[i] & SP[i] the same as SE[i-1] & SP[i-1]
    if (i > 0) {
      currentDupSeSp = se[i] == se[i - 1] && sp[i] == sp[i - 1];
    }
    duplicated[i] = currentDupPred || currentDupSeSp;
  }
  
  // Ensure the sum of !duplicated == thresholds.size()
  // Todo: maybe remove this in a future version of pROC...
  size_t nNonDuplicated = 0;
  for (size_t i = 0; i < npredictors; ++i) {
    if (! duplicated[i]) ++nNonDuplicated;
  }
  if (nNonDuplicated != nthresholds - 1) {
    stop("Bug in pROC: fast algorithm (C++ version) computed an incorrect number of sensitivities and specificities.");
  }
  
  // Push the non-duplicated SE / SP into the final vectors
  vector<double> finalSe(nthresholds), finalSp(nthresholds);
  finalSe[nthresholds - 1] = 0;
  finalSp[nthresholds - 1] = 1;
  size_t nextFinalIdx = 0, nextIdx = npredictors - 1;
  while (nextFinalIdx < nthresholds - 1) {
    if (!duplicated[nextIdx]) {
      finalSe[nextFinalIdx] = se[nextIdx];
      finalSp[nextFinalIdx] = sp[nextIdx];
      ++nextFinalIdx;
    }
    --nextIdx;
  }

  List ret;
  ret["se"] = finalSe;
  ret["sp"] = finalSp;
  
  return ret;
}
