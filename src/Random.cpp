#include <vector>
#include <algorithm>
#include <Rcpp.h>
using std::vector;
using Rcpp::NumericVector;


void setRandomIdx(const size_t size, vector<size_t>& idxVector) {
  NumericVector idx = Rcpp::runif(size, 0, size);
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
  NumericVector idx = Rcpp::runif(totalSize, 0, totalSize);
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