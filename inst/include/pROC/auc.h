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

#include <string>
#include <vector>
#include <pROC/ROC.h>

namespace pROC {
	struct AucParams {
		bool partial, focusOnSp, correct;
		double from, to;
		
		/** Default constructor: full AUC */
		AucParams(): partial(false), focusOnSp(true), correct(false), from(.9), to(1.0) {}
		explicit AucParams(double aFrom, double aTo, std::string aFocus, bool aCorrect):
		                   partial(true), focusOnSp(aFocus == "specificity"), correct(aCorrect), from(aFrom), to(aTo) {}
		explicit AucParams(double aFrom, double aTo, bool aFocusOnSp, bool aCorrect):
		                   partial(true), focusOnSp(aFocusOnSp), correct(aCorrect), from(aFrom), to(aTo) {}
	};
	
	double computeAuc(const std::vector<double>& se, const std::vector<double>& sp, const AucParams& aucParams = AucParams());
	double computeAuc(const std::vector<double>& se, const std::vector<double>& sp, 
	                  const bool partial = false, const double from = 0.9, const double to = 1, const bool focusOnSp = true, const bool correct = false);
	double computeAuc(const std::vector<double>& se, const std::vector<double>& sp, 
	                  const bool partial = false, const double from = 0.9, const double to = 1, const std::string focus = "specificity", const bool correct = false);
	
	template<class Predictor>
	double computeAuc(const pROC::ROC<Predictor>&, const AucParams& aucParams = AucParams());
	template<class Predictor>
	double computeAuc(const pROC::ROC<Predictor>&, 
	                  const bool partial = false, const double from = 0.9, const double to = 1, const bool focusOnSp = true, const bool correct = false);
	template<class Predictor>
	double computeAuc(const pROC::ROC<Predictor>&, 
	                  const bool partial = false, const double from = 0.9, const double to = 1, const std::string focus = "specificity", const bool correct = false);

	double computeFullAuc(const std::vector<double>& se, const std::vector<double>& sp);
	double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, const AucParams& aucParams);
	double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, 
	                         const double from = 0.9, const double to = 1, const std::string focus = "specificity", const bool correct = false);
	double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, 
	                         const double from = 0.9, const double to = 1, const bool focusOnSp = true, const bool correct = false);

	template<class Predictor>
	double computeFullAuc(const pROC::ROC<Predictor>&);
	template<class Predictor>
	double computePartialAuc(const pROC::ROC<Predictor>&, const AucParams& aucParams);
	template<class Predictor>
	double computePartialAuc(const pROC::ROC<Predictor>&, 
	                         const double from = 0.9, const double to = 1, const std::string focus = "specificity", const bool correct = false);
	template<class Predictor>
	double computePartialAuc(const pROC::ROC<Predictor>&, 
	                         const double from = 0.9, const double to = 1, const bool focusOnSp = true, const bool correct = false);
	                         
	/** Computes the corrected pAUC from from to to */
	double correctPartialAuc(const double pAUC, const double from, const double to);
}

#include <pROC/auc_tpl.h>
