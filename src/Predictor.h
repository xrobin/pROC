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
#pragma once 

#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept> // std::out_of_range

/** A Predictor behaves like a concatenated vector of cases and controls
 * The first indices represents the controls, the last ones the controls
 */

class Predictor {
  public:
  	const Rcpp::NumericVector& controls;
  	const Rcpp::NumericVector& cases;
  	const int nControls, nCases, nTotal;
    Predictor(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases): 
              controls(someControls), cases(someCases),
              nControls(someControls.size()), nCases(someCases.size()), nTotal(nControls + nCases) {}
    double operator[] (const int anIdx) const {
      return anIdx < nControls ? controls[anIdx] : cases[anIdx - nControls];
    }

    std::vector<int> getOrder(std::string direction = ">") const;
    
    bool isControl(const int anIdx) const {
      return anIdx < nControls;
    }

    bool isCase(const int anIdx) const {
      return anIdx >= nControls;
    }
    
    bool isValid(const int anIdx) const { // Is there a predictor at this index? 
      return anIdx < nTotal;
    }
    
    /*double at(const size_t anIdx) const {
    	if (isValid(anIdx)) return this[anIdx];
    	throw std::out_of_range("Out of range!");
    }*/
};


/** ResampledPredictor: derived class of Predictor, that takes additional resampling indices */
class ResampledPredictor: public Predictor {
  const std::vector<int>& controlsIdx;
  const std::vector<int>& casesIdx;
  public:
    ResampledPredictor(const Predictor& somePredictor, const std::vector<int>& someControlsIdx, const std::vector<int>& someCasesIdx):
              Predictor(somePredictor), controlsIdx(someControlsIdx), casesIdx(someCasesIdx) {}

    double operator[] (const int anIdx) const {
      return anIdx < nControls ? Predictor::operator[] (controlsIdx[anIdx]) : Predictor::operator[] (casesIdx[anIdx - nControls] + nControls);
    }

    std::vector<int> getOrder(std::string direction = ">") const;
};


/** PredictorComparator: A sorter for two vectors, specifically cases & controls 
 * PredictorReverseComparator: reverse sort!

 * Example usage: 
 * vector<int> index(controls.size() + cases.size());
 * std::iota(index.begin(), index.end(), 0);
 * std::sort(index.begin(), index.end(), ControlCasesComparator(controls, cases));
 */
template <class P> class PredictorComparator{
   const P& predictor;
 public:
   PredictorComparator(const P& somePredictor) : predictor(somePredictor) {}
   bool operator() (int i, int j) const {
     return predictor[i] < predictor[j];
   }
};

template <class P> class PredictorReverseComparator{
   const P& predictor;
 public:
   PredictorReverseComparator(const P& somePredictor) : predictor(somePredictor) {}
   bool operator() (int i, int j) const {
     return predictor[j] < predictor[i];
   }
};


/** getOrder 
 * Get the order (indices) of a Predictor or ResampledPredictor
 */

template <class P> std::vector<int> getPredictorOrder(const P& predictor, const std::string direction = ">") {
    std::vector<int> index(predictor.nTotal);
    std::iota(index.begin(), index.end(), 0);
    if (direction == ">") {
      std::sort(index.begin(), index.end(), PredictorComparator<P>(predictor));
    }
    else {
      std::sort(index.begin(), index.end(), PredictorReverseComparator<P>(predictor));
    }
    
    return index;
}