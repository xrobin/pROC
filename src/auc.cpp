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

AucParams::AucParams(const List& l) {
  try {
    bool test = l["partial.auc"];
    if (! test) {
      partial = false;
    }
    else {
      stop("Reached a line that should be unreachable (partial.auc=TRUE). Please report this bug to the maintainer of pROC. Please type packageDescription(\"pROC\", fields=\"Maintainer\") to obtain this information.");
    }
  }
  catch (const Rcpp::not_compatible&) {
    NumericVector test = l["partial.auc"];
    partial = true;
    from = test[0];
    to = test[1];
    if (to > from) {
      std::swap(from, to);
    }
    correct = l["partial.auc.correct"];
    string paucFocus = l["partial.auc.focus"];
    if (paucFocus == "specificity") {
      focusOnSp = true;
    }
    else if (paucFocus == "sensitivity") {
      focusOnSp = false;
    }
    else {
      string errMsg = string("Invalid partial.auc.focus: ") + paucFocus + ". This probably denotes a bug in pROC. If so, please type packageDescription(\"pROC\", fields=\"Maintainer\") to and report it to the maintainer.";
      stop(errMsg);
    }
  }
}



double aucCC(const vector<double>& controls, const vector<double>& cases, const AucParams& aucParams) {
    // Compute SE/SP of sample
    vector<double> thresholds = computeThresholds(controls, cases);
    pair<vector<double>, vector<double>> sesp = computeSeSp(thresholds, controls, cases);
    return computeAuc(sesp, aucParams);
}

double aucCC(const vector<double>& controls, const vector<double>& cases, 
             const vector<size_t>& controlsIdx, const vector<size_t>& casesIdx,
             const AucParams& aucParams) {
    // Compute SE/SP of sample
    vector<double> thresholds = computeThresholds(controls, cases, controlsIdx, casesIdx);
    pair<vector<double>, vector<double>> sesp = computeSeSp(thresholds, controls, cases, controlsIdx, casesIdx);
    return computeAuc(sesp, aucParams);
}

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

double computePartialAuc(const vector<double>& se, const vector<double>& sp, const AucParams& aucParams) {
  double auc = 0;
  
  // Copy se/sp into y/x
  // If focus == "sp" we can just get a reference, otherwise we must make a copy
  // TODO: figure out how to avoid this copy ()
  const vector<double>& x = aucParams.focusOnSp ? sp : getReversedVector(se);
  const vector<double>& y = aucParams.focusOnSp ? se : getReversedVector(sp);
  
  bool beforeRange = true;
  for (size_t i = 0; i < x.size() - 1; ++i) {
    if (x[i] > aucParams.from) { // We're before the range: nothing to do
      continue;
    }
    else {
      // We are in range
      if (beforeRange) { // but the previous wasn't
        // Are we still in range?
        if (x[i] < aucParams.to) {
          // special case: partial interpolation between 2 out-of-range bounds
          double diff_horiz = aucParams.from - aucParams.to;
          double proportion_before = (x[i-1] - aucParams.from) / (x[i-1] - x[i]);
          double proportion_after  = (aucParams.to - x[i]) / (x[i-1] - x[i]);
          double y_interpolated_before = y[i-1] + proportion_before * (y[i] - y[i-1]);
          double y_interpolated_after = y[i] - proportion_after * (y[i] - y[i-1]);
          std::cout << "adding partial: " << diff_horiz * (y_interpolated_before + y_interpolated_after)/2 << std::endl;
          auc += diff_horiz * (y_interpolated_before + y_interpolated_after) / 2;
          break; // We're done - we're already past the range
        }
        // add previous partial 
        if (i > 0) { // if first element is in range, just do normal processing
          if (x[i] - x[i-1] == 0 || x[i] == aucParams.from) { // if no horizontal span from last point, or we're exactly at the from point, auc += 0, just go to next iteration
            continue;
          }
          double proportion = (aucParams.from - x[i - 1]) / (x[i] - x[i-1]);
          double y_interpolated = y[i-1] + proportion * (y[i] - y[i-1]);
          std::cout << "adding previous partial: " << (aucParams.from - x[i]) * (y[i] + y_interpolated) / 2 << std::endl;
          auc += (aucParams.from - x[i]) * (y[i] + y_interpolated) / 2; 
        }
        beforeRange = false;
      }
      
      if (x[i+1] < aucParams.to) {
        //  add last partial and stop
        if (x[i] - x[i+1] != 0 && x[i] != aucParams.to) {
          // x[i] - x[i+1] != 0: no horizontal span: auc += 0 anyway
          // x[i] != aucParams.to: we're already exactly at the "to" (last) point: nothing to add
          double proportion = (x[i] - aucParams.to) / (x[i] - x[i+1]);
          double y_interpolated = y[i] + proportion * (y[i+1] - y[i]);
          std::cout << "adding last partial: " << (x[i] - aucParams.to) * (y[i] + y_interpolated) / 2 << std::endl;
          auc += (x[i] - aucParams.to) * (y[i] + y_interpolated) / 2;
        }
        break;
      }
      else {
        // Normal add everythhing betwen i and i + 1
        std::cout << "adding full range: " << (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2 << std::endl;
        auc += (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2;
      }
    }
  }

  return auc;
}

double computeAuc(const pair<vector<double>, vector<double>>& sesp, const AucParams& aucParams) {
  const vector<double> se = sesp.first;
  const vector<double> sp = sesp.second;
  
  if (aucParams.partial) {
    return computePartialAuc(se, sp, aucParams);
  }
  else {
    return computeFullAuc(se, sp);
  }
}

// [[Rcpp::export]]
double computeAuc(const std::vector<double>& se, const std::vector<double>& sp, const Rcpp::List& aucParamsList) {
  AucParams aucParams(aucParamsList);
  if (aucParams.partial) {
    return computePartialAuc(se, sp, aucParams);
  }
  else {
    return computeFullAuc(se, sp);
  }
}
