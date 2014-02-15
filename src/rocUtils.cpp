#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>
#include <Rcpp.h>
#include "rocUtils.h"
using std::vector;
using std::pair;

/*void setRandomSample(const vector<double>& from, vector<double>& to) {
  Rcpp::NumericVector idx = Rcpp::runif(from.size(), 0, from.size());
  for (size_t i = 0; i < from.size(); ++i) {
    to[i] = from[idx(i)];
  }
}*/

void setRandomIdx(const size_t size, vector<size_t>& idxVector) {
  Rcpp::NumericVector idx = Rcpp::runif(size, 0, size);
  idxVector.clear();
  idxVector.insert(idxVector.begin(), idx.begin(), idx.end());
}

void setRandomNonStratifiedSample(const size_t controlsSize, const size_t casesSize,
                             vector<size_t>& controlsIdx, vector<size_t>& casesIdx) {
  // relevant sizes
  size_t totalSize = controlsSize + casesSize;
  // prepare return values
  controlsIdx.clear();
  casesIdx.clear();
  // take random index from R
  Rcpp::NumericVector idx = Rcpp::runif(totalSize, 0, totalSize);
  // sort it
  std::sort(idx.begin(), idx.end()); // avoid branch prediction fail?
  // find the first element corresponding to a case
  auto firstCaseIdx = std::lower_bound(idx.begin(), idx.end(), controlsSize);
  // Insert idx into controlsIdx and casesIdx
  controlsIdx.insert(controlsIdx.begin(), idx.begin(), firstCaseIdx);
  casesIdx.insert(casesIdx.begin(), firstCaseIdx, idx.end());
  // casesIdx has values too large by controlsSize
  for (size_t& caseIdxVal: casesIdx) {
    caseIdxVal -= controlsSize;
  };
}

void makeUniqueInPlace(vector<double>& thresholds) {
  // Sort and remove consecutive duplicate values
  std::sort(thresholds.begin(), thresholds.end());
  std::vector<double>::iterator it = std::unique(thresholds.begin(), thresholds.end());
  // The vector still contains duplicate values the end of the vector. The new end is in 'it', use it to resize the vector.
  thresholds.resize(std::distance(thresholds.begin(), it));
}

vector<double> computeThresholds(const vector<double>& predictor) {
  vector<double> thresholds = predictor; // Hold thresholds in a copy of predictor
  makeUniqueInPlace(thresholds);
  return thresholds;
}

vector<double> computeThresholds(const vector<double>& controls, const vector<double>& cases) {
  vector<double> thresholds;
  thresholds.insert(thresholds.begin(), cases.begin(), cases.end());
  thresholds.insert(thresholds.begin(), controls.begin(), controls.end());
  makeUniqueInPlace(thresholds);
  return thresholds;
}

vector<double> computeThresholds(const vector<double>& controls, const vector<double>& cases, const std::vector<size_t>& controlsIdx, const std::vector<size_t>& casesIdx) {
  vector<double> thresholds;
  thresholds.reserve(controlsIdx.size() + casesIdx.size());
  // Insert cases
  for (size_t idx: casesIdx) {
    thresholds.push_back(cases[idx]);
  }
  // Insert controls
  for (size_t idx: controlsIdx) {
    thresholds.push_back(controls[idx]);
  }
  makeUniqueInPlace(thresholds);
  return thresholds;
}

pair<vector<double>, vector<double>> computeSeSp(const vector<double>& thresholds, const vector<double>& controls, const vector<double>& cases) {
  vector<double> se(thresholds.size()),
                 sp(thresholds.size());
  long tp, tn;
  
  for (size_t t = 0; t < thresholds.size(); ++t) {
      double threshold = thresholds[t];
      tp = 0;
      for (size_t i = 0; i < cases.size(); ++i) {
        if (cases[i] >= threshold) {
          tp++;
        }
      }
      se[t] = (double) tp / (double) cases.size();

      tn = 0;
      for (size_t j = 0; j < controls.size(); ++j) {
        if (controls[j] < threshold) {
          tn++;
        }
      }
      sp[t] = (double) tn / (double) controls.size();
    }
  
  return std::make_pair(se, sp);
}



pair<vector<double>, vector<double>> computeSeSp(const vector<double>& thresholds, const vector<double>& controls, const vector<double>& cases,
                                                 const std::vector<size_t>& controlsIdx, const std::vector<size_t>& casesIdx) {
  vector<double> se(thresholds.size()), 
                 sp(thresholds.size());
  long tp, tn;
  
  for (size_t t = 0; t < thresholds.size(); ++t) {
      double threshold = thresholds[t];
      tp = 0;
      for (size_t idx: casesIdx) {
        if (cases[idx] >= threshold) {
          tp++;
        }
      }
      se[t] = (double) tp / (double) cases.size();

      tn = 0;
      for (size_t idx: controlsIdx) {
        if (controls[idx] < threshold) {
          tn++;
        }
      }
      sp[t] = (double) tn / (double) controls.size();
    }
  
  return std::make_pair(se, sp);
}


