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
	char direction;
	
	void computeSeSp();
	double fullAUC() {return computeFullAuc(sensitivity, specificity);}
	double partialAUC(double from, double to, std::string focus, bool correct) {
		return computePartialAuc(sensitivity, specificity, from, to, focus, correct);
	}
	double partialAUC(double from, double to, bool focusOnSp, bool correct) {
		return computePartialAuc(sensitivity, specificity, from, to, focusOnSp, correct);
	}
	
	public:
		ROC(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, char aDirection):
			predictor(someControls, someCases), thresholds(computeThresholds(predictor)), direction(aDirection) {
				computeSeSp();
			}
		ROC(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, char aDirection,
		    std::vector<double> aSensitivity, std::vector<double> aSpecificity):
			predictor(someControls, someCases),sensitivity(aSensitivity), specificity(aSpecificity),
			thresholds(computeThresholds(predictor)), direction(aDirection) {}
			
		double auc(bool partial = false, double from = 0.9, double to = 1, std::string focus = "specificity", bool correct = false) {
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
		
		std::vector<double> getSensitivity() const {return sensitivity;}
		std::vector<double> getSpecificity() const {return sensitivity;}
		std::vector<double> getThresholds() const {return thresholds;}
};

