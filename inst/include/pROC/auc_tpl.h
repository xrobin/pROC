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
#include <pROC/auc.h>
#include <pROC/ROC.h>

template<class Predictor>
double pROC::computeAuc(const pROC::ROC<Predictor>& roc, const pROC::AucParams& aucParams) {
	return computeAuc(roc.getSensitivity(), roc.getSpecificity(), aucParams);
}

template<class Predictor>
double pROC::computeAuc(const pROC::ROC<Predictor>& roc, 
                  bool partial, double from, double to, bool focusOnSp, bool correct) {
	return computeAuc(roc.getSensitivity(), roc.getSpecificity(), partial, from, to, focusOnSp, correct);
}

template<class Predictor>
double pROC::computeAuc(const pROC::ROC<Predictor>& roc, 
                  bool partial, double from, double to, std::string focus, bool correct) {
	return computeAuc(roc.getSensitivity(), roc.getSpecificity(), partial, from, to, focus, correct);
}

template<class Predictor>
double pROC::computeFullAuc(const pROC::ROC<Predictor>& roc) {
	return pROC::computeFullAuc(roc.getSensitivity(), roc.getSpecificity());
}

template<class Predictor>
double pROC::computePartialAuc(const pROC::ROC<Predictor>& roc, const pROC::AucParams& aucParams) {
	return pROC::computePartialAuc(roc.getSensitivity(), roc.getSpecificity(), aucParams);
}

double pROC::computePartialAuc(const pROC::ROC<Predictor>& roc, 
                         double from, double to, std::string focus, bool correct) {
	return computePartialAuc(roc.getSensitivity(), roc.getSpecificity(), from, to, focus, correct);
}

template<class Predictor>
double pROC::computePartialAuc(const pROC::ROC<Predictor>& roc, 
                         double from, double to, bool focusOnSp, bool correct) {
	return computePartialAuc(roc.getSensitivity(), roc.getSpecificity(), from, to, focusOnSp, correct);
}
