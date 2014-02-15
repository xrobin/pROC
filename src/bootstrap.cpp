#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "rocUtils.h"
#include "auc.h"

using std::vector;
using std::pair;
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<double> bootstrapAucStratified(const size_t bootN, const std::vector<double> controls, const std::vector<double> cases) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);

  size_t controlsSize(controls.size()),
         casesSize(cases.size());
  vector<size_t> controlsIdx(controlsSize),
                 casesIdx(casesSize);

  for (size_t i = 0; i < bootN; i++) {
    // Select random sample
    setRandomIdx(controlsSize, controlsIdx);
    setRandomIdx(casesSize, casesIdx);
  
    // Compute AUC
    aucs.push_back(aucCC(controls, cases, controlsIdx, casesIdx));
  }
  
  return aucs;
}


// [[Rcpp::export]]
std::vector<double> bootstrapAucNonStratified(const size_t bootN, const std::vector<double> controls, const std::vector<double> cases) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);
  
  vector<size_t> controlsIdx, casesIdx;

  for (size_t i = 0; i < bootN; i++) {
    // Select random sample
    setRandomNonStratifiedSample(controls.size(), cases.size(), controlsIdx, casesIdx);
  
    // Compute AUC
    aucs.push_back(aucCC(controls, cases, controlsIdx, casesIdx));
  }
  
  return aucs;
}

// // [[Rcpp::export]]
//std::vector<double> bootstrapPaucStratified(const size_t bootN, std::vector<double> controls, std::vector<double> cases) {
