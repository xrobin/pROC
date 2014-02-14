#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>
#include <Rcpp.h>
using std::vector;
using std::pair;

void setRandomSample(const vector<double>& from, vector<double>& to) {
  Rcpp::NumericVector idx = Rcpp::runif(from.size(), 0, from.size());
  for (size_t i = 0; i < from.size(); ++i) {
    to[i] = from[idx(i)];
  }
}

void setRandomUnpairedSample(const vector<double>& fromControls, const vector<double>& fromCases,
                             vector<double>& toControls, vector<double>& toCases) {
  // relevant sizes
  size_t totalSize = fromControls.size() + fromCases.size();
  size_t controlsSize = fromControls.size();
  // prepare return values
  toControls.clear();
  toControls.reserve(controlsSize); // best guess on the size.
  toCases.clear();
  toCases.reserve(fromCases.size());
  // take random index from R
  Rcpp::NumericVector idx = Rcpp::runif(totalSize, 0, totalSize);
  std::sort(idx.begin(), idx.end()); // avoid branch prediction fail?
  for (size_t id: idx) { // cast to size_t
    if (id < controlsSize) {
      toControls.push_back(fromControls[id]);
    }
    else {
      toCases.push_back(fromCases[id - controlsSize]);
    }
  }
  // remove potential unused reserved elements before returning
  toControls.shrink_to_fit();
  toCases.shrink_to_fit();
}

vector<double> computeThresholds(const vector<double>& predictor) {
  // Hold thresholds in a new vector
  vector<double> thresholds = predictor;
  
  // Sort and remove consecutive duplicate values
  std::sort(thresholds.begin(), thresholds.end());
  std::vector<double>::iterator it = std::unique(thresholds.begin(), thresholds.end());
  // The vector still contains duplicate values the end of the vector. The new end is in 'it', use it to resize the vector.
  thresholds.resize(std::distance(thresholds.begin(), it));
  
  return thresholds;
}

vector<double> computeThresholds(const vector<double>& controls, const vector<double>& cases) {
  vector<double> predictor;
  predictor.insert(predictor.end(), cases.begin(), cases.end());
  predictor.insert(predictor.end(), controls.begin(), controls.end());
  return computeThresholds(predictor);
}

pair<vector<double>, vector<double>> computeSeSp(const vector<double>& thresholds, const vector<double>& controls, const vector<double>& cases) {
  vector<double> se(thresholds.size()), sp(thresholds.size());
//  std::pair<vector<double>, vector<double>> sesp;
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


