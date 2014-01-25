#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "rocUtils.h"
#include "auc.h"

using std::vector;
using std::pair;
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<double> bootstrapAucStratified(const std::vector<double>& controls, const std::vector<double>& cases, const size_t bootN) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);
  
  vector<double> sampleControls(controls.size()), 
                      sampleCases(cases.size());

  for (size_t i = 0; i < bootN; i++) {
    // Select random sample
    setRandomSample(controls, sampleControls);
    setRandomSample(cases, sampleCases);
  
    // Compute AUC
    aucs.push_back(aucCC(sampleControls, sampleCases));
  }
  
  return aucs;
}