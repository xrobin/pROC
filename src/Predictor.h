#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>


/** A Predictor behaves like a concatenated vector of cases and controls
 * The first indices represents the controls, the last ones the controls
 */

class Predictor {
  const Rcpp::NumericVector& controls;
  const Rcpp::NumericVector& cases;
  const size_t nControls, nCases, nTotal;
  public:
    Predictor(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases): 
              controls(someControls), cases(someCases),
              nControls(someControls.size()), nCases(someCases.size()), nTotal(nControls + nCases) {}
    double operator[] (const size_t anIdx) const {
      return anIdx < nControls ? controls[anIdx] : cases[anIdx - nControls];
    }
    
    std::vector<size_t> getOrder() const;
    std::vector<size_t> getOrder(std::string) const;

    bool isControl(const size_t anIdx) const {
      return anIdx < nControls;
    }

    bool isCase(const size_t anIdx) const {
      return anIdx >= nControls;
    }
    
    bool isValid(const size_t anIdx) const { // Is there a predictor at this index? 
      return anIdx < nTotal;
    }
    
    
};


/** PredictorComparator: A sorter for two vectors, specifically cases & controls 
 * PredictorReverseComparator: reverse sort!

 * Example usage: 
 * vector<size_t> index(controls.size() + cases.size());
 * std::iota(index.begin(), index.end(), 0);
 * std::sort(index.begin(), index.end(), ControlCasesComparator(controls, cases));
 */
class PredictorComparator{
   const Predictor& predictor;
 public:
   PredictorComparator(const Predictor& somePredictor) : predictor(somePredictor) {}
   bool operator() (size_t i, size_t j) const {
     return predictor[i] < predictor[j];
   }
};

class PredictorReverseComparator{
   const Predictor& predictor;
 public:
   PredictorReverseComparator(const Predictor& somePredictor) : predictor(somePredictor) {}
   bool operator() (size_t i, size_t j) const {
     return predictor[j] < predictor[i];
   }
};