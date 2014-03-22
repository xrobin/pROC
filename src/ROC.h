#pragma once

#include <Rcpp.h>
#include <vector>
#include <string>

#include "Predictor.h"
#include "rocUtils.h" // computeThresholds

template <typename PredictorType>
class ROC {
	PredictorType predictor;
	std::vector<double> sensitivity, specificity, thresholds;
	char direction;
	
	void computeSeSp();
	double fullAUC();
	double partialAUC(double from = 0.9, double to = 1, std::string focus = "specificity", bool correct = false);
	
	public:
		ROC(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, char aDirection): 
			predictor(someControls, someCases), thresholds(computeThresholds(predictor)), direction(aDirection) {
				computeSeSp();
			}
		double auc(bool partial = false, double from = 0.9, double to = 1, std::string focus = "specificity", bool correct = false) {
			if (partial) return partialAUC(from, to, focus, correct);
			else return fullAUC();
		}
		
		std::vector<double> getSensitivity() {return sensitivity;}
		std::vector<double> getSpecificity() {return sensitivity;}
		std::vector<double> getThresholds() {return thresholds;}
};

