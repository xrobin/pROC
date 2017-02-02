/* pROC: Tools Receiver operating characteristic (ROC curves) with
   (partial) area under the curve, confidence intervals and comparison. 
   Copyright (C) 2014 Xavier Robin

   * modified by Stefan Siegert <stefan_siegert@gmx.de>, January 2016

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
List delongPlacementsCpp_old(List roc) {  
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

bool _cmp(std::pair<int, double> l, std::pair<int, double> r) {
  return l.second < r.second;
}

// [[Rcpp::export]]
List delongPlacementsCpp(List roc) {  

  int i, j, k, m, n, mdupl, ndupl, L;

  std::vector<double> cases = roc["cases"];
  std::vector<double> controls = roc["controls"];
  m = cases.size();
  n = controls.size();
  L = m + n;

  // concatenate cases and controls into a vector of L pairs of the form
  // (index, value), also save class labels (1 for cases, 0 for controls)
  std::vector< std::pair<int, double> > Z;
  std::vector< bool > labels;
  for (i = 0; i < m; i++) {
    Z.push_back(std::pair<int, double>(i, cases.at(i)));
    labels.push_back(true);
  }
  for (j = 0; j < n; j++) {
    Z.push_back(std::pair<int, double>(m+j, controls.at(j)));
    labels.push_back(false);
  }

  // sort Z from smallest to largest value, so Z holds the order indices and
  // order statistics of all classifiers
  std::sort(Z.begin(), Z.end(), _cmp);

  // the following calculates the "Delong-placements" X and Y in a single pass
  // over the vector Z, instead of having to double loop over all pairs of
  // (X_i, Y_j)
  std::vector< double > XY(L, 0.0);  // vector to hold the unnormalised X and Y values
  std::vector< int > X_inds, Y_inds; // temporary vectors to save indices of duplicates
  m = n = i = 0;                     // initialisation
  while (i < L) {
    X_inds.clear();
    Y_inds.clear();
    mdupl = ndupl = 0;
    while(1) {
      j = Z.at(i).first;
      if (labels.at(j)) {
        mdupl++;
        X_inds.push_back(j);
      } else {
        ndupl++;
        Y_inds.push_back(j);
      }
      if (i == L-1) {
        break;
      }
      if (Z.at(i).second != Z.at(i+1).second) {
        break;
      }
      i++;
    }
    for (k = 0; k < mdupl; k++) {
      XY.at(X_inds.at(k)) = n + ndupl/2.0;
    }
    for (k = 0; k < ndupl; k++) {
      XY.at(Y_inds.at(k)) = m + mdupl/2.0;
    }
    n += ndupl;
    m += mdupl;
    i++;
  }

  double sum = 0.0;
  std::vector<double> X, Y;

  for (i = 0; i < L; i++) {
    if (labels.at(i)) {
      sum += XY.at(i);
      X.push_back(XY.at(i) / n);
    } else {
      Y.push_back(1.0 - XY.at(i) / m);
    }
  }

  List ret;
  ret["theta"] = sum / (m*n);
  ret["X"] = X;
  ret["Y"] = Y;
  return(ret);
}
