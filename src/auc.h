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

struct AucParams {
  bool partial, focusOnSp, correct;
  double from, to;

  AucParams(): partial(false) {}
  explicit AucParams(const Rcpp::List& l); // non trivial constructor, defined in .cpp file
  explicit AucParams(double aFrom = 0.9, double aTo = 1, std::string aFocus = "specificity", bool aCorrect = false):
	partial(true), focusOnSp(aFocus == "specificity"), correct(aCorrect), from(aFrom), to(aTo) {}
  explicit AucParams(double aFrom = 0.9, double aTo = 1, bool aFocusOnSp = true, bool aCorrect = false):
	partial(true), focusOnSp(aFocusOnSp), correct(aCorrect), from(aFrom), to(aTo) {}
};


/*double aucCC(const Rcpp::NumericVector& controls, const Rcpp::NumericVector& cases, const AucParams& aucParams);
double aucCC(const Rcpp::NumericVector& controls, const Rcpp::NumericVector& cases, 
             const std::vector<int>& controlsIdx, const std::vector<int>& casesIdx, 
             const AucParams& aucParams);
*/
double computeAuc(const std::pair<std::vector<double>, std::vector<double>>&, const AucParams& aucParams);


double computeFullAuc(const std::vector<double>& se, const std::vector<double>& sp);
double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, const AucParams& aucParams);
double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, 
	double from = 0.9, double to = 1, std::string focus = "specificity", bool correct = false);
double computePartialAuc(const std::vector<double>& se, const std::vector<double>& sp, 
	double from = 0.9, double to = 1, bool focusOnSp = true, bool correct = false);
