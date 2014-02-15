#include <Rcpp.h>
#include <vector>
#include <utility>
#include "auc.h"
#include "rocUtils.h"

using std::vector;
using std::pair;
using namespace Rcpp;

double aucCC(const std::vector<double>& controls, const std::vector<double>& cases) {
    // Compute SE/SP of sample
    vector<double> thresholds = computeThresholds(controls, cases);
    pair<vector<double>, vector<double>> sesp = computeSeSp(thresholds, controls, cases);
    return computeAuc(sesp);
}

double aucCC(const std::vector<double>& controls, const std::vector<double>& cases, 
             const std::vector<size_t>& controlsIdx, const std::vector<size_t>& casesIdx) {
    // Compute SE/SP of sample
    vector<double> thresholds = computeThresholds(controls, cases, controlsIdx, casesIdx);
    pair<vector<double>, vector<double>> sesp = computeSeSp(thresholds, controls, cases, controlsIdx, casesIdx);
    return computeAuc(sesp);
}

double computeAuc(const pair<vector<double>, vector<double>>& sesp) {
  const vector<double> se = sesp.first;
  const vector<double> sp = sesp.second;
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
