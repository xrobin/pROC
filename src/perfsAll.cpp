/* pROC: Tools Receiver operating characteristic (ROC curves) with
   (partial) area under the curve, confidence intervals and comparison. 
   Copyright (C) 2014 Xavier Robin

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
using namespace Rcpp;

// [[Rcpp::export]]
List rocUtilsPerfsAllC(NumericVector thresholds, NumericVector controls, NumericVector cases, std::string direction) {
  NumericVector se(thresholds.size());
  NumericVector sp(thresholds.size());
  long tp, tn;
  long i; // iterator over cases & controls
  
  if (direction == ">") {
    for (long t = 0; t < thresholds.size(); t++) {
      if (t % 100 == 0) Rcpp::checkUserInterrupt();
      double threshold = thresholds(t);
        tp = 0;
        for (i = 0; i < cases.size(); i++) {
          if (cases(i) <= threshold) {
            tp++;
        }
      }
      se(t) = (double)tp / cases.size();
      tn = 0;
      for (i = 0; i < controls.size(); i++) {
        if (controls(i) > threshold) {
          tn++;
        }
      }
      sp(t) = (double)tn / controls.size();
    }
  }
 else {
    for (long t = 0; t < thresholds.size(); t++) {
      if (t % 100 == 0) Rcpp::checkUserInterrupt();
      double threshold = thresholds(t);
      tp = 0;
      for (i = 0; i < cases.size(); i++) {
        if (cases(i) >= threshold) {
          tp++;
        }
      }
      se(t) = (double)tp / cases.size();
      long tn = 0;
      for (i = 0; i < controls.size(); i++) {
        if (controls(i) < threshold) {
          tn++;
        }
      }
      sp(t) = (double)tn / controls.size();
    }
  }
  List ret;
  ret["se"] = se;
  ret["sp"] = sp;
  return(ret);
}
