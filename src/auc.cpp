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
#include <vector> // std::vector
#include <utility> // std::pair
#include <algorithm> // std::swap
#include <string> // std::string
#include "auc.h"
#include "rocUtils.h"

using std::vector;
using std::pair;
using std::string;
using Rcpp::List;
using Rcpp::stop;
using Rcpp::NumericVector;


/*double aucCC(const NumericVector& controls, const NumericVector& cases, const AucParams& aucParams) {
    // Compute SE/SP of sample
    vector<double> thresholds = computeThresholds(controls, cases);
    pair<vector<double>, vector<double>> sesp = computeSeSpPair(thresholds, controls, cases);
    return computeAuc(sesp, aucParams);
}

double aucCC(const NumericVector& controls, const NumericVector& cases, 
             const vector<int>& controlsIdx, const vector<int>& casesIdx,
             const AucParams& aucParams) {
    // Compute SE/SP of sample
    vector<double> thresholds = computeThresholds(controls, cases, controlsIdx, casesIdx);
    pair<vector<double>, vector<double>> sesp = computeSeSpPair(thresholds, controls, cases, controlsIdx, casesIdx);
    return computeAuc(sesp, aucParams);
}*/

double computeFullAuc(const vector<double>& se, const vector<double>& sp) {
  double auc = 0;
  size_t lastElement = se.size() - 1;
  
  // Handle first element separately
  auc += (sp[0] - 0) * (se[0] + 1) / 2;

  for (size_t i = 0; i < lastElement; ++i) {
    auc += (sp[i + 1] - sp[i]) * (se[i + 1] + se[i]) / 2;
  }
  
  // Handle last element separately
  auc += (1 - sp[lastElement]) * (se[lastElement]) / 2;
  
  return auc;
}


double computePartialAuc(const vector<double>& se, const vector<double>& sp, 
                         double from, double to, string focus, bool correct) {
	return computePartialAuc(se, sp, AucParams(from, to, focus, correct));
}

double computePartialAuc(const vector<double>& se, const vector<double>& sp, 
                         double from, double to, bool focusOnSp, bool correct) {
	return computePartialAuc(se, sp, AucParams(from, to, focusOnSp, correct));
}

double computePartialAuc(const vector<double>& se, const vector<double>& sp, const AucParams& aucParams) {
  double auc = 0;
  
  // Copy se/sp into y/x
  // If focus == "sp" we can just get a reference, otherwise we must make a copy
  // TODO: figure out how to avoid this copy ()
  const vector<double>& x = aucParams.focusOnSp ? sp : getReversedVector(se);
  const vector<double>& y = aucParams.focusOnSp ? se : getReversedVector(sp);
  
  // iterate backwards until we reach the start of the range
  size_t i = x.size() - 1;
  while (x[i] > aucParams.from) {
    --i;
  }
  // i is now the first element of the range
  
  // special case: we're already beyond the range
  // We need to perform partial interpolation between 2 out-of-range bounds
  if (x[i] < aucParams.to) {
    double diff_horiz = aucParams.from - aucParams.to;
    double proportion_before = (x[i + 1] - aucParams.from) / (x[i + 1] - x[i]);
    double proportion_after  = (aucParams.to - x[i]) / (x[i + 1] - x[i]);
    double y_interpolated_before = y[i + 1] + proportion_before * (y[i] - y[i + 1]);
    double y_interpolated_after = y[i] - proportion_after * (y[i] - y[i + 1]);
    return diff_horiz * (y_interpolated_before + y_interpolated_after) / 2;
  }
  // Ok now we're really in range.
  // Do we need to add the part of AUC before the range?
  // Only if we're NOT the last element AND not exactly "from"
  if (x[i] < 1 && x[i] < aucParams.from) {
          double proportion = (aucParams.from - x[i + 1]) / (x[i] - x[i + 1]);
          double y_interpolated = y[i + 1] + proportion * (y[i] - y[i + 1]);
          auc += (aucParams.from - x[i]) * (y[i] + y_interpolated) / 2; 
  }
  
  // Now we're going down the SP and adding AUCs normally
  // TODO special case where x[i] == from?
  while (i > 0 && x[i - 1] >= aucParams.to) {
    auc += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2;
    --i;
  }
  // We're now either at i = 0 (nothing more to add)
  // Or just at the border of the range (then it was added @ previous step because x[i - 1] >= aucParams.to was true so now x[i] == aucParams.to)
  // or just before being out of range (because if i-1 at the previous step was, it didn't decrement!
  if ( i > 0 && x[i] > aucParams.to) {
    double proportion = (x[i] - aucParams.to) / (x[i] - x[i - 1]);
    double y_interpolated = y[i] + proportion * (y[i - 1] - y[i]);
    auc += (x[i] - aucParams.to) * (y[i] + y_interpolated) / 2;
  }
  
  if (aucParams.correct) {
    stop("Corrected pAUC not implemented in C++");
  }

  return auc;
}

double computeAuc(const pair<vector<double>, vector<double>>& sesp, const AucParams& aucParams = AucParams()) {
  const vector<double> se = sesp.first;
  const vector<double> sp = sesp.second;
  
  if (aucParams.partial) {
    return computePartialAuc(se, sp, aucParams);
  }
  else {
    return computeFullAuc(se, sp);
  }
}

