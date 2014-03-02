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

double MWkernel(double x, double y) {
  if (y == x) {return(0.5);}
	else if (y < x) {return(1);}
	else {return(0);}
}


// [[Rcpp::export]]
List delongPlacementsCpp(List roc) {  
  NumericVector cases = roc["cases"];
	NumericVector controls = roc["controls"];
	std::size_t xsize = cases.size();
	std::size_t ysize = controls.size();
	double sum = 0;
	double MW;
	std::vector<double> X(xsize, 0.0);
	std::vector<double> Y(ysize, 0.0);

	for (std::size_t i = 0; i < xsize; i++) {
	  for (std::size_t j = 0; j < ysize; j++) {
      MW = MWkernel(cases[i], controls[j]);
      X[i] += MW ;
	    Y[j] += MW;
	    sum += MW;
	  }
	}
  for (std::size_t i = 0; i < xsize; i++) {
    X[i] /= ysize;
  }
	for (std::size_t j = 0; j < ysize; j++) {
    Y[j] /= xsize;
	}
	List ret;
	ret["theta"] = sum / xsize / ysize;
	ret["X"] = X;
	ret["Y"] = Y;
	return(ret);
}
