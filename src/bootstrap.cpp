#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include "rocUtils.h"
#include "auc.h"

using std::vector;
using std::pair;
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<double> bootstrapAucStratified(const size_t bootN, std::vector<double> controls, std::vector<double> cases, const std::string& direction) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);
  
  if (direction == ">") {
    for (size_t i = 0; i < controls.size(); ++i) {
      controls[i] = -controls[i];
    }
    for (size_t i = 0; i < cases.size(); ++i) {
      cases[i] = -cases[i];
    }
  }
  
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