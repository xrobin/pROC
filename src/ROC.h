#pragma once

#include <Rcpp.h>
#include <vector>
#include <string>

#include "Predictor.h"
#include "rocUtils.h" // computeThresholds
#include "auc.h"

template <typename PredictorType>
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
			predictor(someControls, someCases), thresholds(computeThresholds(predictor)), direction(aDirection) {
				computeSeSp();
		}
		/** Constructor with sensitivity and specificity pre-computed.
		 * Warning: no attempt is made to check the correctness of the data
		 */
		ROC(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, const std::string& aDirection,
		    std::vector<double> aSensitivity, std::vector<double> aSpecificity):
			predictor(someControls, someCases),sensitivity(aSensitivity), specificity(aSpecificity),
			thresholds(computeThresholds(predictor)), direction(aDirection) {}

		/** Constructor from a Predictor.
		 * Computes sensitivity and specificity
		 */
		ROC(const PredictorType& aPredictor, const std::string& aDirection):
			predictor(aPredictor), thresholds(computeThresholds(aPredictor)), direction(aDirection) {
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
		ROC<ResampledPredictorType> getResampled() {
			ResampledPredictorType aResampledPredictor(predictor);
			ROC<ResampledPredictorType> aResampledROC(aResampledPredictor, direction);
			return aResampledROC;
		}
		/** Resamples the predictor and computes the new SE/SP
		 * Works only if the predictor is a ResampledPredictor.
		 */
		void resample() {
			predictor.resample();
			computeSeSp();
		}
		
		std::vector<double> getSensitivity() const {return sensitivity;}
		std::vector<double> getSpecificity() const {return specificity;}
		std::vector<double> getThresholds() const {return thresholds;}
		char getDirection() const {return direction;}
		PredictorType getPredictor() const {return predictor;}
};

#include "computeSeSp.h"
