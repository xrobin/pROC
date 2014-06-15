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
#include <string>
#include <vector>

#include <pROC/auc.h>
#include <pROC/Predictor.h>
#include <pROC/rocUtils.h> // computeThresholds

namespace pROC {

	template <typename PredictorType = Predictor>
	class ROC {
		PredictorType predictor;
		std::vector<double> sensitivity, specificity, thresholds;
		const std::string direction;
		
		void computeSeSp();
		double fullAUC() {return computeFullAuc(sensitivity, specificity);}
		double partialAUC(double from, double to, std::string focus, bool correct) {
			return computePartialAuc(sensitivity, specificity, from, to, focus, correct);
		}
		double partialAUC(double from, double to, bool focusOnSp, bool correct) {
			return computePartialAuc(sensitivity, specificity, from, to, focusOnSp, correct);
		}
		double partialAUC(AucParams params) {
			return computePartialAuc(sensitivity, specificity, params);
		}
		
		public:
			/** Constructor with controls, cases and a direction.
			 * Computes sensitivity and specificity
			 */
			ROC(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, const std::string& aDirection):
				predictor(someControls, someCases), sensitivity(), specificity(), 
				thresholds(computeThresholds(predictor, aDirection)), direction(aDirection)
				{
					computeSeSp();
			}
			/** Constructor with sensitivity and specificity pre-computed.
			 * Warning: no attempt is made to check the correctness of the data
			 */
			ROC(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, const std::string& aDirection,
			    std::vector<double> aSensitivity, std::vector<double> aSpecificity):
				predictor(someControls, someCases),sensitivity(aSensitivity), specificity(aSpecificity),
				thresholds(computeThresholds(predictor, aDirection)), direction(aDirection) {}
	
			/** Constructor from a Predictor.
			 * Computes sensitivity and specificity
			 */
			ROC(const PredictorType& aPredictor, const std::string& aDirection):
				predictor(aPredictor), sensitivity(), specificity(), thresholds(computeThresholds(aPredictor, aDirection)), 
				direction(aDirection) {
					computeSeSp();
			}
				
			double auc(bool partial = false, double from = 0.9, double to = 1, std::string& focus = "specificity", bool correct = false) {
				if (partial) return partialAUC(from, to, focus, correct);
				else return fullAUC();
			}
			double auc(bool partial = false, double from = 0.9, double to = 1, bool focusOnSp = true, bool correct = false) {
				if (partial) return partialAUC(from, to, focusOnSp, correct);
				else return fullAUC();
			}
			double auc(AucParams params) {
				if (params.partial) return partialAUC(params);
				else return fullAUC();	
			}
					
			/** returns a ROC curve with a resampled predictor. Call as:
			 * ROC<ResampledPredictorStratified> myResampledROC = ROC.getResampled<ResampledPredictorStratified>()
			 * ROC<ResampledPredictorNonStratified> myResampledROC = ROC.getResampled<ResampledPredictorNonStratified>()
			 */
			template <class ResampledPredictorType>
			ROC<ResampledPredictorType> getResampled() const {
				ResampledPredictorType aResampledPredictor(predictor);
				ROC<ResampledPredictorType> aResampledROC(aResampledPredictor, direction);
				return aResampledROC;
			}
			/** Resamples the predictor and computes the new SE/SP
			 * Works only if the predictor is a ResampledPredictor.
			 */
			void resample();
			
			std::vector<double> getSensitivity() const {return sensitivity;}
			std::vector<double> getSpecificity() const {return specificity;}
			std::vector<double> getThresholds() const {return thresholds;}
			char getDirection() const {return direction;}
			PredictorType getPredictor() const {return predictor;}
	};
}

#include <pROC/ROC.tpl>
