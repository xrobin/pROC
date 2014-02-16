#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "Random.h"
#include "auc.h"

using std::vector;
using std::pair;


// [[Rcpp::export]]
std::vector<double> bootstrapAucStratified(const size_t bootN, const std::vector<double> controls, const std::vector<double> cases, const Rcpp::List& aucParamsList) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);
  
  // Get proper AUC params
  AucParams aucParams(aucParamsList);

  size_t controlsSize(controls.size()),
         casesSize(cases.size());
  vector<size_t> controlsIdx(controlsSize),
                 casesIdx(casesSize);

  for (size_t i = 0; i < bootN; i++) {
    // Select random sample
    setRandomIdx(controlsSize, controlsIdx);
    setRandomIdx(casesSize, casesIdx);
  
    // Compute AUC
    aucs.push_back(aucCC(controls, cases, controlsIdx, casesIdx, aucParams));
  }
  
  return aucs;
}


// [[Rcpp::export]]
std::vector<double> bootstrapAucNonStratified(const size_t bootN, const std::vector<double> controls, const std::vector<double> cases, const Rcpp::List& aucParamsList) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);
  
  // Get proper AUC params
  AucParams aucParams(aucParamsList);
  
  vector<size_t> controlsIdx, casesIdx;

  for (size_t i = 0; i < bootN; i++) {
    // Select random sample
    setRandomNonStratifiedSample(controls.size(), cases.size(), controlsIdx, casesIdx);
  
    // Compute AUC
    aucs.push_back(aucCC(controls, cases, controlsIdx, casesIdx, aucParams));
  }
  
  return aucs;
}

// // [[Rcpp::export]]
//std::vector<double> bootstrapPaucStratified(const size_t bootN, std::vector<double> controls, std::vector<double> cases) {
