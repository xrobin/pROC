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