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
#include <utility>
#include <algorithm>
#include <string>
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



double aucCC(const std::vector<double>& controls, const std::vector<double>& cases, const AucParams& aucParams) {
    // Compute SE/SP of sample
    vector<double> thresholds = computeThresholds(controls, cases);
    pair<vector<double>, vector<double>> sesp = computeSeSp(thresholds, controls, cases);
    return computeAuc(sesp, aucParams);
}

double aucCC(const std::vector<double>& controls, const std::vector<double>& cases, 
             const std::vector<size_t>& controlsIdx, const std::vector<size_t>& casesIdx,
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

template<typename T> vector<T> getReversedVector(const vector<T>& x) {
  vector<T> y(x);
  std::reverse(y.begin(), y.end());
  return y;
}

double computePartialAuc(const vector<double>& se, const vector<double>& sp, const AucParams& aucParams) {
  double auc = 0;
  
  // Copy se/sp into y/x
  // If focus == "sp" we can just get a reference, otherwise we must make a copy
  // TODO: figure out how to avoid this copy ()
  const vector<double>& x = aucParams.focusOnSp ? sp : getReversedVector(se);
  const vector<double>& y = aucParams.focusOnSp ? se : getReversedVector(sp);
  
  //auto first_in = aucParams.from < 1 std::lower_bound(idx.begin(), idx.end(), controlsSize) : ;
  //auto first_out = std::lower_bound(idx.begin(), idx.end(), controlsSize);
  
  bool beforeRange = true;
  for (size_t i = 0; i < x.size() - 1; ++i) {
    if (x[i] > aucParams.from) {
      continue;
    }
    else {
      // We are in range
      if (beforeRange) { // but the previous wasn't
        // TODO: add previous partial
        ...;
        beforeRange = false;
      }
      
      // TODO: check if we must handle i == 0 differently. Probably in if (beforeRange)
      if (x[i+1] < aucParams.to) {
        // TODO: add last partial and stop
        ...;
        break;
      }
      else {
        auc += (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2;
      }
    }
  }


/*
    
    // focus == "sp" > x == "sp" y == "se"
    
    # find the SEs and SPs in the interval
    x.inc <- x[x <= partial.auc[1] & x >= partial.auc[2]]
    y.inc <- y[x <= partial.auc[1] & x >= partial.auc[2]]
    # compute the AUC strictly in the interval
    diffs.x <- x.inc[-1] - x.inc[-length(x.inc)]
    means.vert <- (y.inc[-1] + y.inc[-length(y.inc)])/2
    auc <- sum(means.vert * diffs.x)
    # add the borders:
    if (length(x.inc) == 0) { # special case: the whole AUC is between 2 se/sp points. Need to interpolate from both
      diff.horiz <- partial.auc[1] - partial.auc[2]
      # determine indices
      idx.hi <- match(FALSE, x < partial.auc[1])
      idx.lo <- idx.hi - 1
      # proportions
      proportion.hi <- (x[idx.hi] - partial.auc[1]) / (x[idx.hi] - x[idx.lo])
      proportion.lo <- (partial.auc[2] - x[idx.lo]) / (x[idx.hi] - x[idx.lo])
      # interpolated y's
      y.hi <- y[idx.hi] + proportion.hi * (y[idx.lo] - y[idx.hi])
      y.lo <- y[idx.lo] - proportion.lo * (y[idx.lo] - y[idx.hi])
      # compute AUC
      mean.vert <- (y.hi + y.lo)/2
      auc <- mean.vert*diff.horiz
    }
    else { # if the upper limit is not exactly present in SPs, interpolate
      if (!(partial.auc[1] %in% x.inc)) {
        # find the limit indices
        idx.out <- match(FALSE, x < partial.auc[1])
        idx.in <- idx.out - 1
        # interpolate y
        proportion <- (partial.auc[1] - x[idx.out]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.out] + proportion * (y[idx.in] - y[idx.out])
        # add to AUC
        auc <- auc + (partial.auc[1] - x[idx.in]) * (y[idx.in] + y.interpolated)/2
      }
      if (!(partial.auc[2] %in% x.inc)) { # if the lower limit is not exactly present in SPs, interpolate
        # find the limit indices in and out
        #idx.out <- length(x) - match(TRUE, rev(x) < partial.auc[2]) + 1
        idx.out <- match(TRUE, x > partial.auc[2]) - 1
        idx.in <- idx.out + 1
        # interpolate y
        proportion <- (x[idx.in] - partial.auc[2]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.in] + proportion * (y[idx.out] - y[idx.in])
        # add to AUC
        auc <- auc + (x[idx.in] - partial.auc[2]) * (y[idx.in] + y.interpolated)/2
      }
    }
    */
  return 0;
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
