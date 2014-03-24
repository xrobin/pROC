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

#include <vector>
#include <string>

struct AucParams {
  bool partial, focusOnSp, correct;
  double from, to;

  /** Default constructor: full AUC */
  AucParams(): partial(false) {} 
  explicit AucParams(double aFrom, double aTo, std::string aFocus, bool aCorrect):
	partial(true), focusOnSp(aFocus == "specificity"), correct(aCorrect), from(aFrom), to(aTo) {}
  explicit AucParams(double aFrom, double aTo, bool aFocusOnSp, bool aCorrect):
	partial(true), focusOnSp(aFocusOnSp), correct(aCorrect), from(aFrom), to(aTo) {}
};

double computeAuc(const std::pair<std::vector<double>, std::vector<double>>&, const AucParams& aucParams);

double computeFullAuc(const std::vector<double>& se, const std::vector<double>& sp);
double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, const AucParams& aucParams);
double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, 
	double from = 0.9, double to = 1, std::string focus = "specificity", bool correct = false);
double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, 
	double from = 0.9, double to = 1, bool focusOnSp = true, bool correct = false);
