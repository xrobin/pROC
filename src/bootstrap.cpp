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
#include <algorithm>
#include "Random.h"
#include "auc.h"
#include "RcppConversions.h"

using std::vector;
using std::pair;
using Rcpp::NumericVector;


// [[Rcpp::export]]
std::vector<double> bootstrapAucStratified(const size_t bootN, const NumericVector controls, const NumericVector cases, const Rcpp::List& aucParamsList) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);
  
  // Get proper AUC params
  AucParams aucParams = as <AucParams>(aucParamsList);

  int controlsSize(controls.size()),
         casesSize(cases.size());
  vector<int> controlsIdx(controlsSize),
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
std::vector<double> bootstrapAucNonStratified(const int bootN, const NumericVector controls, const NumericVector cases, const Rcpp::List& aucParamsList) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);
  
  // Get proper AUC params
  AucParams aucParams = as <AucParams>(aucParamsList);
  
  vector<int> controlsIdx, casesIdx;

  for (int i = 0; i < bootN; i++) {
    // Select random sample
    setRandomNonStratifiedSample(controls.size(), cases.size(), controlsIdx, casesIdx);
  
    // Compute AUC
    aucs.push_back(aucCC(controls, cases, controlsIdx, casesIdx, aucParams));
  }
  
  return aucs;
}

// // [[Rcpp::export]]
//std::vector<double> bootstrapPaucStratified(const size_t bootN, std::vector<double> controls, std::vector<double> cases) {
